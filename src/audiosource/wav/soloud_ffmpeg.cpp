/*
SoLoud audio engine - ffmpeg interface
Copyright (c) 2025 William Horvath

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

   1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.

   2. Altered source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.

   3. This notice may not be removed or altered from any source
   distribution.
*/

#if defined(WITH_FFMPEG) && __has_include(<libavcodec/avcodec.h>) && ((defined(_WIN32) || defined(_WIN64)) || defined(__linux__))

#include <algorithm>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "soloud_ffmpeg.h"
#include "soloud_ffmpeg_load.h"
#include "soloud_file.h"

using namespace SoLoud::FFmpeg::FFmpegLoader::FFmpegFuncs;

namespace SoLoud::FFmpeg
{

static const unsigned int IO_BUFFER_SIZE = 65536;

struct FFmpegIOContext
{
	File *file;
	int pos;
};

struct FFmpegDecoder
{
	AVFormatContext *formatContext;
	AVCodecContext *codecContext;
	SwrContext *swrContext;
	AVFrame *frame;
	AVPacket *packet;
	FFmpegIOContext *ioContext;
	AVIOContext *avioContext;
	unsigned char *ioBuffer;

	int audioStreamIndex;
	unsigned int channels;
	unsigned int sampleRate;
	unsigned long long totalFrames;
	unsigned long long currentFrame;
	unsigned int encoderDelaySamples;
	AVChannelLayout channelLayout;

	// frame buffering for partial reads
	float *frameBuffer;
	unsigned int frameBufferSize;
	unsigned int frameBufferOffset;
	unsigned int frameBufferValidSamples;

	// pts-based position tracking
	long long lastFramePts;
	bool validPts;
	unsigned long long bufferedFramePosition;
	bool bufferedFramePositionValid;

	bool endOfStream;
	bool seekable;
};

static int ffmpeg_read_func(void *opaque, unsigned char *buf, int buf_size)
{
	FFmpegIOContext *ctx = (FFmpegIOContext *)opaque;
	if (!ctx || !ctx->file || !buf || buf_size <= 0)
		return AVERROR_EOF;

	int remaining = ctx->file->length() - ctx->pos;
	if (remaining <= 0)
		return AVERROR_EOF;

	int toRead = (buf_size < remaining) ? buf_size : remaining;
	ctx->file->seek(ctx->pos);
	int read = ctx->file->read(buf, toRead);
	if (read > 0)
		ctx->pos += read;

	return read > 0 ? read : AVERROR_EOF;
}

static int64_t ffmpeg_seek_func(void *opaque, int64_t offset, int whence)
{
	FFmpegIOContext *ctx = (FFmpegIOContext *)opaque;
	if (!ctx || !ctx->file)
		return -1;

	switch (whence)
	{
	case SEEK_SET:
		ctx->pos = (int)offset;
		break;
	case SEEK_CUR:
		ctx->pos += (int)offset;
		break;
	case SEEK_END:
		ctx->pos = ctx->file->length() + (int)offset;
		break;
	case AVSEEK_SIZE:
		return ctx->file->length();
	}

	ctx->pos = std::max(0, std::min<int>(ctx->pos, ctx->file->length()));
	return ctx->pos;
}

static unsigned long long ptsToFrames(FFmpegDecoder *decoder, int64_t pts)
{
	if (pts == AV_NOPTS_VALUE)
		return decoder->currentFrame;

	AVStream *stream = decoder->formatContext->streams[decoder->audioStreamIndex];

	// convert pts from stream time_base to seconds, then to sample frames
	double timeSeconds = (double)pts * stream->time_base.num / stream->time_base.den;
	unsigned long long framePosition = (unsigned long long)(timeSeconds * decoder->sampleRate);

	// subtract encoder delay to match our adjusted frame counting
	if (framePosition >= decoder->encoderDelaySamples)
		framePosition -= decoder->encoderDelaySamples;

	return framePosition;
}

static bool decodeNextFrame(FFmpegDecoder *decoder)
{
	if (decoder->endOfStream)
		return false;

	while (true)
	{
		int ret = av_read_frame(decoder->formatContext, decoder->packet);
		if (ret < 0)
		{
			// flush decoder at end of stream
			avcodec_send_packet(decoder->codecContext, nullptr);
			if (avcodec_receive_frame(decoder->codecContext, decoder->frame) >= 0)
			{
				// store frame pts if available
				if (decoder->frame->pts != AV_NOPTS_VALUE)
				{
					decoder->lastFramePts = decoder->frame->pts;
					decoder->validPts = true;
				}

				return decoder->frame->nb_samples > 0;
			}
			decoder->endOfStream = true;
			return false;
		}

		if (decoder->packet->stream_index != decoder->audioStreamIndex)
		{
			av_packet_unref(decoder->packet);
			continue;
		}

		if (avcodec_send_packet(decoder->codecContext, decoder->packet) >= 0)
		{
			av_packet_unref(decoder->packet);

			if (avcodec_receive_frame(decoder->codecContext, decoder->frame) >= 0)
			{
				// store frame pts if available
				if (decoder->frame->pts != AV_NOPTS_VALUE)
				{
					decoder->lastFramePts = decoder->frame->pts;
					decoder->validPts = true;
				}

				return decoder->frame->nb_samples > 0;
			}
		}
		else
		{
			av_packet_unref(decoder->packet);
		}
	}
}

static void convertFrameToBuffer(FFmpegDecoder *decoder)
{
	if (!decoder->frame || decoder->frame->nb_samples <= 0)
		return;

	unsigned int samples = decoder->frame->nb_samples;
	unsigned int channels = decoder->channels;

	// ensure frame buffer is large enough
	unsigned int requiredSize = samples * channels;
	if (decoder->frameBufferSize < requiredSize)
	{
		delete[] decoder->frameBuffer;
		decoder->frameBuffer = new float[requiredSize];
		decoder->frameBufferSize = requiredSize;
	}

	// update position tracking when we start consuming a new frame
	if (decoder->validPts && decoder->lastFramePts != AV_NOPTS_VALUE)
	{
		decoder->bufferedFramePosition = ptsToFrames(decoder, decoder->lastFramePts);
		decoder->bufferedFramePositionValid = true;

		// sync manual tracking if it's invalid (after seek)
		if (decoder->currentFrame == 0 && decoder->bufferedFramePosition > 0)
			decoder->currentFrame = decoder->bufferedFramePosition;
	}

	if (decoder->codecContext->sample_fmt == AV_SAMPLE_FMT_FLTP)
	{
		// direct copy from planar float data
		for (unsigned int ch = 0; ch < channels; ch++)
		{
			float *src = (float *)decoder->frame->data[ch];
			float *dst = decoder->frameBuffer + ch * samples;
			memcpy(dst, src, samples * sizeof(float));
		}
	}
	else
	{
		// convert via resampler
		const uint8_t *inputData[MAX_CHANNELS];
		if (av_sample_fmt_is_planar(decoder->codecContext->sample_fmt))
		{
			for (unsigned int ch = 0; ch < channels; ch++)
			{
				inputData[ch] = decoder->frame->data[ch];
			}
		}
		else
		{
			inputData[0] = decoder->frame->data[0];
		}

		uint8_t *outputPtrs[MAX_CHANNELS];
		for (unsigned int ch = 0; ch < channels; ch++)
		{
			outputPtrs[ch] = (uint8_t *)(decoder->frameBuffer + ch * samples);
		}

		int converted = swr_convert(decoder->swrContext, outputPtrs, samples, inputData, samples);
		if (converted > 0)
			samples = converted;
		else
			samples = 0;
	}

	decoder->frameBufferValidSamples = samples;
	decoder->frameBufferOffset = 0;
}

FFmpegDecoder *open(File *aFile, bool oneshot)
{
	if (!FFmpegLoader::isAvailable() || !aFile)
		return nullptr;

	auto *decoder = new FFmpegDecoder();
	memset((void *)decoder, 0, sizeof(FFmpegDecoder));

	AVStream *stream = nullptr;
	AVCodecParameters *codecpar = nullptr;
	const AVCodec *codec = nullptr;
	AVDictionary *format_opts = nullptr;
	double duration = 0.0;
	unsigned long long totalFrames = 0;
	unsigned int encoderDelaySamples = 0;

	// setup io context
	decoder->ioContext = new FFmpegIOContext();
	decoder->ioContext->file = aFile;
	decoder->ioContext->pos = 0;

	decoder->ioBuffer = (unsigned char *)av_malloc(IO_BUFFER_SIZE);
	if (!decoder->ioBuffer)
		goto cleanup;

	decoder->avioContext = avio_alloc_context(decoder->ioBuffer, IO_BUFFER_SIZE, 0, decoder->ioContext, ffmpeg_read_func, nullptr, ffmpeg_seek_func);
	if (!decoder->avioContext)
		goto cleanup;

	decoder->formatContext = avformat_alloc_context();
	if (!decoder->formatContext)
		goto cleanup;

	decoder->formatContext->pb = decoder->avioContext;

	if (!oneshot)
	{
		av_dict_set(&format_opts, "flags", "low_delay", 0);
		av_dict_set(&format_opts, "seek2any", "1", 0);
		av_dict_set(&format_opts, "accurate_seek", "1", 0);
	}

	if (avformat_open_input(&decoder->formatContext, nullptr, nullptr, &format_opts) < 0)
	{
		av_dict_free(&format_opts);
		goto cleanup;
	}
	av_dict_free(&format_opts);

	if (avformat_find_stream_info(decoder->formatContext, nullptr) < 0)
		goto cleanup;

	// find audio stream
	decoder->audioStreamIndex = -1;
	for (unsigned int i = 0; i < decoder->formatContext->nb_streams; i++)
	{
		if (decoder->formatContext->streams[i]->codecpar->codec_type == AVMEDIA_TYPE_AUDIO)
		{
			decoder->audioStreamIndex = i;
			break;
		}
	}

	if (decoder->audioStreamIndex == -1)
		goto cleanup;

	// setup decoder
	stream = decoder->formatContext->streams[decoder->audioStreamIndex];
	codecpar = stream->codecpar;
	codec = avcodec_find_decoder(codecpar->codec_id);
	if (!codec)
		goto cleanup;

	decoder->codecContext = avcodec_alloc_context3(codec);
	if (!decoder->codecContext)
		goto cleanup;

	if (avcodec_parameters_to_context(decoder->codecContext, codecpar) < 0)
		goto cleanup;

	decoder->codecContext->pkt_timebase = stream->time_base;
	decoder->codecContext->thread_count = 1; // we're already a thread, it's not worth spawning yet more threads and keeping track of that synchronization as well
	decoder->codecContext->thread_type = 0;

	if (codec->capabilities & AV_CODEC_CAP_VARIABLE_FRAME_SIZE)
		decoder->codecContext->flags2 |= AV_CODEC_FLAG2_CHUNKS;

	if (avcodec_open2(decoder->codecContext, codec, nullptr) < 0)
		goto cleanup;

	// setup resampler if needed
	if (codecpar->format != AV_SAMPLE_FMT_FLTP)
	{
		decoder->swrContext = swr_alloc();
		if (!decoder->swrContext)
			goto cleanup;

		if (av_channel_layout_copy(&decoder->channelLayout, &codecpar->ch_layout) < 0)
			goto cleanup;

		if (decoder->channelLayout.nb_channels == 0)
		{
			av_channel_layout_uninit(&decoder->channelLayout);
			av_channel_layout_default(&decoder->channelLayout, codecpar->ch_layout.nb_channels);
		}

		av_opt_set_chlayout(decoder->swrContext, "in_chlayout", &decoder->channelLayout, 0);
		av_opt_set_chlayout(decoder->swrContext, "out_chlayout", &decoder->channelLayout, 0);
		av_opt_set_int(decoder->swrContext, "in_sample_rate", codecpar->sample_rate, 0);
		av_opt_set_int(decoder->swrContext, "out_sample_rate", codecpar->sample_rate, 0);
		av_opt_set_sample_fmt(decoder->swrContext, "in_sample_fmt", decoder->codecContext->sample_fmt, 0);
		av_opt_set_sample_fmt(decoder->swrContext, "out_sample_fmt", AV_SAMPLE_FMT_FLTP, 0);

		if (swr_init(decoder->swrContext) < 0)
			goto cleanup;
	}
	else
	{
		if (av_channel_layout_copy(&decoder->channelLayout, &codecpar->ch_layout) < 0)
			goto cleanup;

		if (decoder->channelLayout.nb_channels == 0)
		{
			av_channel_layout_uninit(&decoder->channelLayout);
			av_channel_layout_default(&decoder->channelLayout, codecpar->ch_layout.nb_channels);
		}
	}

	decoder->frame = av_frame_alloc();
	decoder->packet = av_packet_alloc();
	if (!decoder->frame || !decoder->packet)
		goto cleanup;

	decoder->channels = decoder->channelLayout.nb_channels;
	decoder->sampleRate = codecpar->sample_rate;
	decoder->currentFrame = 0;
	decoder->endOfStream = false;
	decoder->seekable = (decoder->formatContext->pb && decoder->formatContext->pb->seekable);

	// get encoder delay from the stream we just opened
	if (stream->start_time != AV_NOPTS_VALUE && stream->start_time > 0)
	{
		double startTimeSeconds = (double)stream->start_time * stream->time_base.num / stream->time_base.den;
		encoderDelaySamples = (unsigned int)floor(startTimeSeconds * decoder->sampleRate);
	}
	decoder->encoderDelaySamples = encoderDelaySamples;

	// calculate total frames - subtract encoder delay for BASS-like behavior
	if (stream->duration > 0 && stream->time_base.den > 0)
		duration = (double)stream->duration * stream->time_base.num / stream->time_base.den;
	else if (decoder->formatContext->duration > 0)
		duration = (double)decoder->formatContext->duration / AV_TIME_BASE;
	else
		duration = 300.0; // fallback

	totalFrames = (unsigned int)floor(duration * decoder->sampleRate);
	decoder->totalFrames = totalFrames > encoderDelaySamples ? totalFrames - encoderDelaySamples : totalFrames;

	decoder->frameBuffer = nullptr;
	decoder->frameBufferSize = 0;
	decoder->frameBufferOffset = 0;
	decoder->frameBufferValidSamples = 0;

	// initialize pts tracking
	decoder->lastFramePts = AV_NOPTS_VALUE;
	decoder->validPts = false;
	decoder->bufferedFramePosition = 0;
	decoder->bufferedFramePositionValid = false;

	// establish proper starting position accounting for encoder delay
	if (!seekToFrame(decoder, 0))
	{
#ifdef _DEBUG
		SoLoud::logStdout("ffmpeg: failed to establish initial position\n");
#endif
		goto cleanup;
	}

	return decoder;

cleanup:
	close(decoder);
	return nullptr;
}

void close(FFmpegDecoder *decoder)
{
	if (!decoder)
		return;

	// cleanup format context
	if (decoder->formatContext)
	{
		avformat_close_input(&decoder->formatContext);
		decoder->formatContext = nullptr;
	}

	// cleanup codec context
	if (decoder->codecContext)
	{
		avcodec_free_context(&decoder->codecContext);
		decoder->codecContext = nullptr;
	}

	// cleanup resampler
	if (decoder->swrContext)
	{
		swr_free(&decoder->swrContext);
		decoder->swrContext = nullptr;
	}

	if (decoder->frame)
	{
		av_frame_free(&decoder->frame);
		decoder->frame = nullptr;
	}

	if (decoder->packet)
	{
		av_packet_free(&decoder->packet);
		decoder->packet = nullptr;
	}

	// cleanup channel layout
	av_channel_layout_uninit(&decoder->channelLayout);

	// cleanup io context
	unsigned char *bufferToFree = nullptr;
	if (decoder->avioContext)
	{
		bufferToFree = decoder->avioContext->buffer;
		avio_context_free(&decoder->avioContext);
		if (bufferToFree)
			av_free(bufferToFree);
	}
	else if (decoder->ioBuffer)
	{
		av_free(decoder->ioBuffer);
	}

	delete decoder->ioContext;
	delete[] decoder->frameBuffer;
	delete decoder;
}

bool seekToFrame(FFmpegDecoder *decoder, unsigned long long frameIndex)
{
	if (!decoder || !decoder->seekable)
	{
#ifdef _DEBUG
		SoLoud::logStdout("ffmpeg: seek rejected - decoder %s, seekable: %s\n", decoder ? "valid" : "null", decoder && decoder->seekable ? "yes" : "no");
#endif
		return false;
	}

	unsigned long long currentPos = getCurrentFrame(decoder);
	if (frameIndex == currentPos)
	{
#ifdef _DEBUG
		SoLoud::logStdout("ffmpeg: seek skipped - already at frame %llu\n", frameIndex);
#endif
		return true;
	}

	AVStream *stream = decoder->formatContext->streams[decoder->audioStreamIndex];
	unsigned long long adjustedFrameIndex = frameIndex + decoder->encoderDelaySamples;
	long long targetTimestamp = av_rescale_q((int64_t)adjustedFrameIndex, {1, (int)decoder->sampleRate}, stream->time_base);

#ifdef _DEBUG
	SoLoud::logStdout("ffmpeg: seeking from frame=%llu (%.3fs) to frame %llu (%.3fs), adjusted=%llu, timestamp=%lld\n", currentPos, (double)currentPos / decoder->sampleRate,
	       frameIndex, (double)frameIndex / decoder->sampleRate, adjustedFrameIndex, targetTimestamp);
#endif

	// seek to keyframe at or before target
	int ret = avformat_seek_file(decoder->formatContext, decoder->audioStreamIndex, INT64_MIN, targetTimestamp, targetTimestamp, AVSEEK_FLAG_BACKWARD);
	if (ret < 0)
	{
#ifdef _DEBUG
		SoLoud::logStdout("ffmpeg: seek failed with error %d\n", ret);
#endif
		return false;
	}

	// reset decoder state
	avcodec_flush_buffers(decoder->codecContext);

	// flush resampler buffers
	if (decoder->swrContext)
	{
		long long delay = swr_get_delay(decoder->swrContext, decoder->sampleRate);
		if (delay > 0)
		{
#ifdef _DEBUG
			SoLoud::logStdout("ffmpeg: flushing %lld samples from resampler\n", delay);
#endif
			swr_drop_output(decoder->swrContext, (int)delay);
		}
	}

	// invalidate frame buffer
	decoder->frameBufferOffset = 0;
	decoder->frameBufferValidSamples = 0;
	decoder->endOfStream = false;

	// reset position tracking
	decoder->currentFrame = 0;
	decoder->lastFramePts = AV_NOPTS_VALUE;
	decoder->validPts = false;
	decoder->bufferedFramePosition = 0;
	decoder->bufferedFramePositionValid = false;

	// decode first frame to establish where we actually landed
	if (!decodeNextFrame(decoder))
	{
#ifdef _DEBUG
		SoLoud::logStdout("ffmpeg: failed to decode frame after seek\n");
#endif
		return false;
	}
	convertFrameToBuffer(decoder); // frame is now in buffer, ready to be consumed
	av_frame_unref(decoder->frame);

	unsigned long long keyframePos = getCurrentFrame(decoder);

#ifdef _DEBUG
	SoLoud::logStdout("ffmpeg: landed at keyframe %llu (%.3fs), target is frame %llu (%.3fs)\n", keyframePos, (double)keyframePos / decoder->sampleRate, frameIndex,
	       (double)frameIndex / decoder->sampleRate);
#endif

	// if we need to advance further, read forward (discarding audio)
	if (keyframePos < frameIndex)
	{
		unsigned long long framesToSkip = frameIndex - keyframePos;

#ifdef _DEBUG
		SoLoud::logStdout("ffmpeg: skipping forward %llu frames to reach target\n", framesToSkip);
#endif

		// temporary buffer for discarded audio
		float *tempBuffer = new float[framesToSkip * decoder->channels];
		unsigned long long framesRead = readFrames(decoder, framesToSkip, tempBuffer);
		delete[] tempBuffer;

		if (framesRead != framesToSkip)
		{
#ifdef _DEBUG
			SoLoud::logStdout("ffmpeg: skip-forward incomplete - skipped %llu/%llu frames\n", framesRead, framesToSkip);
#endif
		}
	}
	else if (keyframePos == frameIndex)
	{
		// exactly where we want to be, frame is ready in buffer
#ifdef _DEBUG
		SoLoud::logStdout("ffmpeg: landed exactly on target frame\n");
#endif
	}
	else
	{
		// this shouldn't happen with AVSEEK_FLAG_BACKWARD
#ifdef _DEBUG
		SoLoud::logStdout("ffmpeg: warning - landed past target by %llu frames\n", keyframePos - frameIndex);
#endif
	}

#ifdef _DEBUG
	SoLoud::logStdout("ffmpeg: seek completed - target frame %llu (%.3fs), actual frame %llu (%.3fs)\n", frameIndex, (double)frameIndex / decoder->sampleRate,
	       getCurrentFrame(decoder), (double)getCurrentFrame(decoder) / decoder->sampleRate);
#endif

	return true;
}

unsigned long long readFrames(FFmpegDecoder *decoder, unsigned long long framesToRead, float *buffer)
{
	if (!decoder || !buffer || framesToRead == 0)
		return 0;

	unsigned long long framesRead = 0;
	unsigned int channels = decoder->channels;

	while (framesRead < framesToRead && !decoder->endOfStream)
	{
		// use buffered frame data if available
		if (decoder->frameBufferValidSamples > decoder->frameBufferOffset)
		{
			unsigned long long available = decoder->frameBufferValidSamples - decoder->frameBufferOffset;
			unsigned long long toCopy = std::min(available, framesToRead - framesRead);

			if (toCopy > 0)
			{
				// copy frame data in planar format
				for (unsigned int ch = 0; ch < channels; ch++)
				{
					float *src = decoder->frameBuffer + ch * decoder->frameBufferValidSamples + decoder->frameBufferOffset;
					float *dst = buffer + ch * framesToRead + framesRead;
					memcpy(dst, src, toCopy * sizeof(float));
				}

				decoder->frameBufferOffset += (unsigned int)toCopy;
				framesRead += toCopy;
				decoder->currentFrame += toCopy;
			}
		}
		else
		{
			// need to decode next frame
			if (!decodeNextFrame(decoder))
				break;

			convertFrameToBuffer(decoder);
			av_frame_unref(decoder->frame);
		}
	}

	return framesRead;
}

result loadToMemory(File *aFile, float **aData, unsigned int *aSampleCount, unsigned int *aChannels, float *aSampleRate)
{
	if (!aFile || !aData || !aSampleCount || !aChannels || !aSampleRate)
		return FILE_LOAD_FAILED;

	constexpr bool oneshot = true;
	FFmpegDecoder *decoder = open(aFile, oneshot);
	if (!decoder)
		return FILE_LOAD_FAILED;

	*aChannels = getChannels(decoder);
	*aSampleRate = (float)getSampleRate(decoder);
	unsigned long long totalFrames = getTotalFrameCount(decoder);

	if (totalFrames == 0)
	{
		close(decoder);
		return FILE_LOAD_FAILED;
	}

	*aSampleCount = (unsigned int)totalFrames;
	*aData = new float[totalFrames * (*aChannels)];

	unsigned long long framesRead = readFrames(decoder, totalFrames, *aData);
	close(decoder);

	if (framesRead != totalFrames)
	{
		delete[] *aData;
		*aData = nullptr;
		return FILE_LOAD_FAILED;
	}

	return SO_NO_ERROR;
}

unsigned int getChannels(FFmpegDecoder *decoder)
{
	return decoder ? decoder->channels : 0;
}

unsigned int getSampleRate(FFmpegDecoder *decoder)
{
	return decoder ? decoder->sampleRate : 0;
}

unsigned long long getTotalFrameCount(FFmpegDecoder *decoder)
{
	return decoder ? decoder->totalFrames : 0;
}

unsigned long long getCurrentFrame(FFmpegDecoder *decoder)
{
	if (!decoder)
		return 0;

	// prefer buffered frame position if we have valid timestamps
	if (decoder->bufferedFramePositionValid)
	{
		// return position of buffered frame + how much we've consumed from it
		return decoder->bufferedFramePosition + decoder->frameBufferOffset;
	}

	// fallback to manual tracking
	return decoder->currentFrame;
}

} // namespace SoLoud::FFmpeg
#endif
