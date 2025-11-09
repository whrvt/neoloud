/*
SoLoud audio engine
Copyright (c) 2013-2018 Jari Komppa

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

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "dr_flac.h"
#include "dr_mp3.h"
#include "dr_wav.h"
#include "stb_vorbis.h"

#include "soloud_mpg123.h"

#include "soloud_ffmpeg.h"
#include "soloud_ffmpeg_load.h"

#include "soloud.h"
#include "soloud_file.h"
#include "soloud_wavstream.h"

namespace SoLoud
{
namespace // static
{

// All dr_* libraries use the same seek origin values: 0=SET, 1=CUR, 2=END
// These are duplicated to avoid UBSan warnings, even though the underlying types are the same
size_t drwav_read_func(void *pUserData, void *pBufferOut, size_t bytesToRead)
{
	File *fp = (File *)pUserData;
	return fp->read((unsigned char *)pBufferOut, (unsigned int)bytesToRead);
}

drwav_bool32 drwav_seek_func(void *pUserData, int offset, drwav_seek_origin origin)
{
	File *fp = (File *)pUserData;

	if (origin == 1) // CUR
		offset += fp->pos();
	else if (origin == 2) // END
		offset += fp->length();

	fp->seek(offset);
	return 1;
}

drwav_bool32 drwav_tell_func(void *pUserData, drwav_int64 *pCursor)
{
	File *fp = (File *)pUserData;
	*pCursor = (drwav_int64)fp->pos();
	return 1;
}

size_t drmp3_read_func(void *pUserData, void *pBufferOut, size_t bytesToRead)
{
	File *fp = (File *)pUserData;
	return fp->read((unsigned char *)pBufferOut, (unsigned int)bytesToRead);
}

drmp3_bool32 drmp3_seek_func(void *pUserData, int offset, drmp3_seek_origin origin)
{
	File *fp = (File *)pUserData;

	if (origin == 1) // CUR
		offset += fp->pos();
	else if (origin == 2) // END
		offset += fp->length();

	fp->seek(offset);
	return 1;
}

drmp3_bool32 drmp3_tell_func(void *pUserData, drmp3_int64 *pCursor)
{
	File *fp = (File *)pUserData;
	*pCursor = (drmp3_int64)fp->pos();
	return 1;
}

size_t drflac_read_func(void *pUserData, void *pBufferOut, size_t bytesToRead)
{
	File *fp = (File *)pUserData;
	return fp->read((unsigned char *)pBufferOut, (unsigned int)bytesToRead);
}

drflac_bool32 drflac_seek_func(void *pUserData, int offset, drflac_seek_origin origin)
{
	File *fp = (File *)pUserData;

	if (origin == 1) // CUR
		offset += fp->pos();
	else if (origin == 2) // END
		offset += fp->length();

	fp->seek(offset);
	return 1;
}

drflac_bool32 drflac_tell_func(void *pUserData, drflac_int64 *pCursor)
{
	File *fp = (File *)pUserData;
	*pCursor = (drflac_int64)fp->pos();
	return 1;
}

} // namespace

WavStreamInstance::WavStreamInstance(WavStream *aParent)
    : mCodec()
{
	mOggFrameSize = 0;
	mParent = aParent;
	mOffset = 0;
	mFile = nullptr;

	if (aParent->mMemFile)
	{
		MemoryFile *mf = new MemoryFile();
		mFile = mf;
		mf->openMem(aParent->mMemFile->getMemPtr(), aParent->mMemFile->length(), false, false);
	}
	else if (aParent->mFilename)
	{
		DiskFile *df = new DiskFile;
		mFile = df;
		df->open(aParent->mFilename);
	}
	else if (aParent->mStreamFile)
	{
		mFile = aParent->mStreamFile;
		mFile->seek(0); // stb_vorbis assumes file offset to be at start of ogg
	}
	else
	{
		return;
	}

	if (mFile)
	{
		if (mParent->mFiletype == WAVSTREAM_WAV)
		{
			mCodec.mWav = new drwav;
			if (!drwav_init(mCodec.mWav, drwav_read_func, drwav_seek_func, drwav_tell_func, (void *)mFile, nullptr))
			{
				delete mCodec.mWav;
				mCodec.mWav = nullptr;
				if (mFile != mParent->mStreamFile)
					delete mFile;
				mFile = nullptr;
			}
		}
		else if (mParent->mFiletype == WAVSTREAM_OGG)
		{
			int e;

			mCodec.mOgg = stb_vorbis_open_file((Soloud_Filehack *)mFile, 0, &e, nullptr);

			if (!mCodec.mOgg)
			{
				if (mFile != mParent->mStreamFile)
					delete mFile;
				mFile = nullptr;
			}
			mOggFrameSize = 0;
			mOggFrameOffset = 0;
			mOggOutputs = nullptr;
		}
		else if (mParent->mFiletype == WAVSTREAM_FLAC)
		{
			mCodec.mFlac = drflac_open(drflac_read_func, drflac_seek_func, drflac_tell_func, (void *)mFile, nullptr);
			if (!mCodec.mFlac)
			{
				if (mFile != mParent->mStreamFile)
					delete mFile;
				mFile = nullptr;
			}
		}
		else if (mParent->mFiletype == WAVSTREAM_MPG123)
		{
			mCodec.mMpg123 = MPG123::open(mFile);
			if (!mCodec.mMpg123)
			{
				if (mFile != mParent->mStreamFile)
					delete mFile;
				mFile = nullptr;
			}
		}
		else if (mParent->mFiletype == WAVSTREAM_DRMP3)
		{
			mCodec.mDrmp3 = new drmp3;
			if (!drmp3_init(mCodec.mDrmp3, drmp3_read_func, drmp3_seek_func, drmp3_tell_func, nullptr, (void *)mFile, nullptr))
			{
				delete mCodec.mDrmp3;
				mCodec.mDrmp3 = nullptr;
				if (mFile != mParent->mStreamFile)
					delete mFile;
				mFile = nullptr;
			}
			else if (mParent->mMp3SeekPointCount > 0 && mParent->mMp3SeekPoints != nullptr)
			{
				drmp3_bind_seek_table(mCodec.mDrmp3, mParent->mMp3SeekPointCount, mParent->mMp3SeekPoints);
			}
		}
		else if (mParent->mFiletype == WAVSTREAM_FFMPEG)
		{
			mCodec.mFfmpeg = FFmpeg::open(mFile);
			if (!mCodec.mFfmpeg)
			{
				if (mFile != mParent->mStreamFile)
					delete mFile;
				mFile = nullptr;
			}
		}
		else
		{
			if (mFile != mParent->mStreamFile)
				delete mFile;
			mFile = nullptr;
			return;
		}
	}
}

WavStreamInstance::~WavStreamInstance()
{
	switch (mParent->mFiletype)
	{
	case WAVSTREAM_OGG:
		if (mCodec.mOgg)
		{
			stb_vorbis_close(mCodec.mOgg);
			mCodec.mOgg = nullptr;
		}
		break;
	case WAVSTREAM_FLAC:
		if (mCodec.mFlac)
		{
			drflac_close(mCodec.mFlac);
			mCodec.mFlac = nullptr;
		}
		break;
	case WAVSTREAM_MPG123:
		if (mCodec.mMpg123)
		{
			MPG123::close(mCodec.mMpg123);
			mCodec.mMpg123 = nullptr;
		}
		break;
	case WAVSTREAM_DRMP3:
		if (mCodec.mDrmp3)
		{
			drmp3_uninit(mCodec.mDrmp3);
			delete mCodec.mDrmp3;
			mCodec.mDrmp3 = nullptr;
		}
		break;
	case WAVSTREAM_WAV:
		if (mCodec.mWav)
		{
			drwav_uninit(mCodec.mWav);
			delete mCodec.mWav;
			mCodec.mWav = nullptr;
		}
		break;
	case WAVSTREAM_FFMPEG:
		if (mCodec.mFfmpeg)
		{
			FFmpeg::close(mCodec.mFfmpeg);
			mCodec.mFfmpeg = nullptr;
		}
		break;
	}
	if (mFile != mParent->mStreamFile)
	{
		delete mFile;
		mFile = nullptr;
	}
}

static unsigned int fillFromOggFrames(stb_vorbis *vorbis, unsigned int &frameSize, unsigned int &frameOffset, float **&frameOutputs, float *buffer,
                                      unsigned int samplesToRead, unsigned int bufferSize, unsigned int channels)
{
	unsigned int totalSamples = 0;

	while (totalSamples < samplesToRead)
	{
		// check if we need a new frame
		if (frameOffset >= frameSize)
		{
			frameSize = stb_vorbis_get_frame_float(vorbis, nullptr, &frameOutputs);
			frameOffset = 0;

			if (frameSize == 0)
				break; // end of stream
		}

		// calculate how many samples we can copy from current frame
		unsigned int samplesInFrame = frameSize - frameOffset;
		unsigned int samplesToCopy = samplesToRead - totalSamples;
		if (samplesToCopy > samplesInFrame)
			samplesToCopy = samplesInFrame;

		// copy samples with planar layout
		for (unsigned int s = 0; s < samplesToCopy; s++)
		{
			for (unsigned int ch = 0; ch < channels; ch++)
			{
				buffer[ch * bufferSize + totalSamples + s] = frameOutputs[ch][frameOffset + s];
			}
		}

		totalSamples += samplesToCopy;
		frameOffset += samplesToCopy;
	}

	return totalSamples;
}

unsigned int WavStreamInstance::getAudio(float *aBuffer, unsigned int aSamplesToRead, unsigned int /* aBufferSize */)
{
	unsigned int offset = 0;
	float tmp[SAMPLE_GRANULARITY * MAX_CHANNELS];
	if (aBuffer == nullptr || mFile == nullptr)
		return 0;
	switch (mParent->mFiletype)
	{
	case WAVSTREAM_FLAC: {
		unsigned int i, j, k;

		for (i = 0; i < aSamplesToRead; i += SAMPLE_GRANULARITY)
		{
			unsigned int blockSize = (aSamplesToRead - i) > SAMPLE_GRANULARITY ? SAMPLE_GRANULARITY : aSamplesToRead - i;
			offset += (unsigned int)drflac_read_pcm_frames_f32(mCodec.mFlac, blockSize, tmp);

			for (j = 0; j < blockSize; j++)
			{
				for (k = 0; k < mChannels; k++)
				{
					aBuffer[k * aSamplesToRead + i + j] = tmp[j * mCodec.mFlac->channels + k];
				}
			}
		}
		mOffset += offset;
		return offset;
	}
	break;
	case WAVSTREAM_MPG123: {
		unsigned int i, j, k;

		for (i = 0; i < aSamplesToRead; i += SAMPLE_GRANULARITY)
		{
			unsigned int blockSize = (aSamplesToRead - i) > SAMPLE_GRANULARITY ? SAMPLE_GRANULARITY : aSamplesToRead - i;
			size_t framesRead = MPG123::readFrames(mCodec.mMpg123, blockSize, tmp);
			offset += (unsigned int)framesRead;

			for (j = 0; j < framesRead; j++)
			{
				for (k = 0; k < mChannels; k++)
				{
					aBuffer[k * aSamplesToRead + i + j] = tmp[j * MPG123::getChannels(mCodec.mMpg123) + k];
				}
			}

			if (framesRead < blockSize)
				break;
		}
		mOffset += offset;
		return offset;
	}
	break;
	case WAVSTREAM_DRMP3: {
		unsigned int i, j, k;

		for (i = 0; i < aSamplesToRead; i += SAMPLE_GRANULARITY)
		{
			unsigned int blockSize = (aSamplesToRead - i) > SAMPLE_GRANULARITY ? SAMPLE_GRANULARITY : aSamplesToRead - i;
			offset += (unsigned int)drmp3_read_pcm_frames_f32(mCodec.mDrmp3, blockSize, tmp);

			for (j = 0; j < blockSize; j++)
			{
				for (k = 0; k < mChannels; k++)
				{
					aBuffer[k * aSamplesToRead + i + j] = tmp[j * mCodec.mDrmp3->channels + k];
				}
			}
		}
		mOffset += offset;
		return offset;
	}
	break;
	case WAVSTREAM_OGG: {
		unsigned int i;

		for (i = 0; i < aSamplesToRead; i += SAMPLE_GRANULARITY)
		{
			unsigned int blockSize = (aSamplesToRead - i) > SAMPLE_GRANULARITY ? SAMPLE_GRANULARITY : aSamplesToRead - i;
			unsigned int samplesRead =
			    fillFromOggFrames(mCodec.mOgg, mOggFrameSize, mOggFrameOffset, mOggOutputs, aBuffer + i, blockSize, aSamplesToRead, mChannels);
			offset += samplesRead;

			if (samplesRead < blockSize)
				break; // end of stream
		}

		mOffset += offset;
		return offset;
	}
	break;
	case WAVSTREAM_WAV: {
		unsigned int i, j, k;

		for (i = 0; i < aSamplesToRead; i += SAMPLE_GRANULARITY)
		{
			unsigned int blockSize = (aSamplesToRead - i) > SAMPLE_GRANULARITY ? SAMPLE_GRANULARITY : aSamplesToRead - i;
			offset += (unsigned int)drwav_read_pcm_frames_f32(mCodec.mWav, blockSize, tmp);

			for (j = 0; j < blockSize; j++)
			{
				for (k = 0; k < mChannels; k++)
				{
					aBuffer[k * aSamplesToRead + i + j] = tmp[j * mCodec.mWav->channels + k];
				}
			}
		}
		mOffset += offset;
		return offset;
	}
	break;
	case WAVSTREAM_FFMPEG: {
		if (mCodec.mFfmpeg)
		{
			offset = (unsigned int)FFmpeg::readFrames(mCodec.mFfmpeg, aSamplesToRead, aBuffer);
			mOffset += offset;
		}
		return offset;
	}
	break;
	}
	return aSamplesToRead;
}

result WavStreamInstance::seek(double aSeconds, float *mScratch, unsigned int mScratchSize)
{
	if (mParent->mFiletype == WAVSTREAM_DRMP3 && mCodec.mDrmp3)
	{
		drmp3_uint64 targetFrame = (drmp3_uint64)floor(aSeconds * mCodec.mDrmp3->sampleRate);
		if (drmp3_seek_to_pcm_frame(mCodec.mDrmp3, targetFrame))
		{
			// Since the position that we just sought to might not be *exactly*
			// the position we asked for, we're re-calculating the position just
			// for the sake of correctness.
			mOffset = (unsigned int)mCodec.mDrmp3->currentPCMFrame;
			double newPosition = static_cast<double>(mOffset) / mBaseSamplerate;
			mStreamPosition = newPosition;
			return SO_NO_ERROR;
		}
		return UNKNOWN_ERROR;
	}
	else if (mParent->mFiletype == WAVSTREAM_MPG123 && mCodec.mMpg123)
	{
		off_t targetFrame = (off_t)floor(aSeconds * mBaseSamplerate);
		off_t seekResult = MPG123::seekToFrame(mCodec.mMpg123, targetFrame);
		if (seekResult >= 0)
		{
			mOffset = (unsigned int)seekResult; // libmpg123 helpfully returns the actual position seeked to
			double newPosition = static_cast<double>(mOffset) / mBaseSamplerate;
			mStreamPosition = newPosition;
			return SO_NO_ERROR;
		}
		return UNKNOWN_ERROR;
	}
	else if (mParent->mFiletype == WAVSTREAM_OGG && mCodec.mOgg)
	{
		unsigned int pos = (unsigned int)floor(mBaseSamplerate * aSeconds);

		if (stb_vorbis_seek(mCodec.mOgg, pos) == 1)
		{
			// only reset frame state after successful seek
			mOggFrameSize = 0;
			mOggFrameOffset = 0;

			// get actual position after seek (may not be exact)
			mOffset = stb_vorbis_get_sample_offset(mCodec.mOgg);
			double newPosition = static_cast<double>(mOffset) / mBaseSamplerate;
			mStreamPosition = newPosition;
			return SO_NO_ERROR;
		}
		return UNKNOWN_ERROR;
	}
	else if (mParent->mFiletype == WAVSTREAM_WAV && mCodec.mWav)
	{
		drwav_uint64 targetFrame = (drwav_uint64)floor(aSeconds * mCodec.mWav->sampleRate);

		if (drwav_seek_to_pcm_frame(mCodec.mWav, targetFrame))
		{
			/* Same as above */
			mOffset = (unsigned int)mCodec.mWav->readCursorInPCMFrames;
			double newPosition = static_cast<double>(mOffset) / mBaseSamplerate;
			mStreamPosition = newPosition;
			return SO_NO_ERROR;
		}
		return UNKNOWN_ERROR;
	}
	else if (mParent->mFiletype == WAVSTREAM_FLAC && mCodec.mFlac)
	{
		drflac_uint64 targetFrame = (drflac_uint64)floor(aSeconds * mCodec.mFlac->sampleRate);

		if (drflac_seek_to_pcm_frame(mCodec.mFlac, targetFrame))
		{
			/* Same as above */
			mOffset = (unsigned int)mCodec.mFlac->currentPCMFrame;
			double newPosition = static_cast<double>(mOffset) / mBaseSamplerate;
			mStreamPosition = newPosition;
			return SO_NO_ERROR;
		}
		return UNKNOWN_ERROR;
	}
	else if (mParent->mFiletype == WAVSTREAM_FFMPEG && mCodec.mFfmpeg)
	{
		unsigned long long targetFrame = (unsigned long long)floor(aSeconds * mBaseSamplerate);

		if (FFmpeg::seekToFrame(mCodec.mFfmpeg, targetFrame))
		{
			mOffset = (unsigned int)FFmpeg::getCurrentFrame(mCodec.mFfmpeg);
			mStreamPosition = static_cast<double>(mOffset) / mBaseSamplerate;
			return SO_NO_ERROR;
		}
		return UNKNOWN_ERROR;
	}
	return AudioSourceInstance::seek(aSeconds, mScratch, mScratchSize);
}

result WavStreamInstance::rewind()
{
	switch (mParent->mFiletype)
	{
	case WAVSTREAM_OGG:
		if (mCodec.mOgg)
		{
			stb_vorbis_seek_start(mCodec.mOgg);
		}
		break;
	case WAVSTREAM_FLAC:
		if (mCodec.mFlac)
		{
			drflac_seek_to_pcm_frame(mCodec.mFlac, 0);
		}
		break;
	case WAVSTREAM_MPG123:
		if (mCodec.mMpg123)
		{
			MPG123::seekToFrame(mCodec.mMpg123, 0);
		}
		break;
	case WAVSTREAM_DRMP3:
		if (mCodec.mDrmp3)
		{
			drmp3_seek_to_pcm_frame(mCodec.mDrmp3, 0);
		}
		break;
	case WAVSTREAM_WAV:
		if (mCodec.mWav)
		{
			drwav_seek_to_pcm_frame(mCodec.mWav, 0);
		}
		break;
	case WAVSTREAM_FFMPEG:
		if (mCodec.mFfmpeg)
		{
			FFmpeg::seekToFrame(mCodec.mFfmpeg, 0);
		}
		break;
	}
	mOffset = 0;
	mStreamPosition = 0.0f;
	return 0;
}

bool WavStreamInstance::hasEnded()
{
	// check if we've reached the estimated sample count
	if (mOffset >= mParent->mSampleCount)
	{
		return true;
	}

	// for codecs that have reliable end-of-stream indicators, check those too
	if ((mParent->mFiletype == WAVSTREAM_DRMP3 && mCodec.mDrmp3 && mCodec.mDrmp3->atEnd) ||
	    (mParent->mFiletype == WAVSTREAM_MPG123 && mCodec.mMpg123 && MPG123::isAtEnd(mCodec.mMpg123)))
	{
		return true;
	}

	// for OGG: if we have no current frame and are at/past estimated end, consider ended
	// this handles cases where stb_vorbis_stream_length_in_samples was inaccurate
	if (mParent->mFiletype == WAVSTREAM_OGG && mOggFrameSize == 0 && mOffset >= mParent->mSampleCount)
	{
		return true;
	}

	return false;
}

WavStream::WavStream(bool preferFFmpeg)
{
	mFilename = nullptr;
	mSampleCount = 0;
	mFiletype = WAVSTREAM_WAV;
	mMemFile = nullptr;
	mStreamFile = nullptr;
	mMp3SeekPoints = nullptr;
	mMp3SeekPointCount = 0;
	mPreferFFmpeg = preferFFmpeg;
}

WavStream::~WavStream()
{
	stop();
	delete[] mFilename;
	delete mMemFile;
	delete[] mMp3SeekPoints;
}

#define MAKEDWORD(a, b, c, d) (((d) << 24) | ((c) << 16) | ((b) << 8) | (a))

result WavStream::loadwav(File *fp)
{
	fp->seek(0);
	drwav decoder;

	if (!drwav_init(&decoder, drwav_read_func, drwav_seek_func, drwav_tell_func, (void *)fp, nullptr))
		return FILE_LOAD_FAILED;

	mChannels = decoder.channels;
	if (mChannels > MAX_CHANNELS)
	{
		mChannels = MAX_CHANNELS;
	}

	mBaseSamplerate = (float)decoder.sampleRate;
	mSampleCount = (unsigned int)decoder.totalPCMFrameCount;
	mFiletype = WAVSTREAM_WAV;
	drwav_uninit(&decoder);

	return SO_NO_ERROR;
}

result WavStream::loadogg(File *fp)
{
	fp->seek(0);
	int e;
	stb_vorbis *v;
	v = stb_vorbis_open_file((Soloud_Filehack *)fp, 0, &e, nullptr);
	if (v == nullptr)
		return FILE_LOAD_FAILED;

	stb_vorbis_comment comment = stb_vorbis_get_comment(v);
	// skip broken "avc[...]"-encoded vorbis files, ffmpeg handles it correctly, stb_vorbis doesn't
	for (int i = 0; i < comment.comment_list_length; i++)
	{
		const char *currcom = nullptr;
		// no comment, allow this file
		if (!(comment.comment_list && (currcom = comment.comment_list[i]) && (*currcom != '\0')))
		{
			break;
		}

		if (strncmp("encoder=", currcom, (sizeof("encoder=") - 1)) == 0)
		{
			currcom += (sizeof("encoder=") - 1);
			if (*currcom != '\0' && (strncmp("avc", currcom, (sizeof("avc") - 1)) == 0))
			{
				return FILE_LOAD_FAILED;
			}
			else
			{
				// we got "encoder=" something other than avc[...], so allow this file
				break;
			}
		}
	}

	stb_vorbis_info info = stb_vorbis_get_info(v);
	mChannels = info.channels;
	if (info.channels > MAX_CHANNELS)
	{
		mChannels = MAX_CHANNELS;
	}
	mBaseSamplerate = (float)info.sample_rate;
	unsigned int samples = stb_vorbis_stream_length_in_samples(v);
	stb_vorbis_close(v);
	mFiletype = WAVSTREAM_OGG;

	// NOTE: stb_vorbis_stream_length_in_samples may be inaccurate for some files
	// we handle this gracefully in hasEnded() and getAudio()
	mSampleCount = samples;

	return 0;
}

result WavStream::loadflac(File *fp)
{
	fp->seek(0);
	drflac *decoder = drflac_open(drflac_read_func, drflac_seek_func, drflac_tell_func, (void *)fp, nullptr);

	if (decoder == nullptr)
		return FILE_LOAD_FAILED;

	mChannels = decoder->channels;
	if (mChannels > MAX_CHANNELS)
	{
		mChannels = MAX_CHANNELS;
	}

	mBaseSamplerate = (float)decoder->sampleRate;
	mSampleCount = (unsigned int)decoder->totalPCMFrameCount;
	mFiletype = WAVSTREAM_FLAC;
	drflac_close(decoder);

	return SO_NO_ERROR;
}

result WavStream::loadmpg123(File *fp)
{
	fp->seek(0);
	MPG123::MPG123Decoder *decoder = MPG123::open(fp);

	if (!decoder)
		return FILE_LOAD_FAILED;

	mChannels = MPG123::getChannels(decoder);
	if (mChannels > MAX_CHANNELS)
	{
		mChannels = MAX_CHANNELS;
	}

	mBaseSamplerate = (float)MPG123::getSampleRate(decoder);
	off_t totalFrames = MPG123::getTotalFrameCount(decoder);

	if (totalFrames <= 0)
	{
		MPG123::close(decoder);
		return FILE_LOAD_FAILED;
	}

	mSampleCount = (unsigned int)totalFrames;
	mFiletype = WAVSTREAM_MPG123;

	MPG123::close(decoder);
	return SO_NO_ERROR;
}

result WavStream::loaddrmp3(File *fp)
{
	fp->seek(0);
	drmp3 decoder;
	if (!drmp3_init(&decoder, drmp3_read_func, drmp3_seek_func, drmp3_tell_func, nullptr, (void *)fp, nullptr))
		return FILE_LOAD_FAILED;

	mChannels = decoder.channels;
	if (mChannels > MAX_CHANNELS)
	{
		mChannels = MAX_CHANNELS;
	}

	drmp3_uint64 samples = drmp3_get_pcm_frame_count(&decoder);

	if (!samples)
	{
		drmp3_uninit(&decoder);
		return FILE_LOAD_FAILED;
	}

	// validate by trying to decode a couple (2) frames, so that we can actually report an error and fall back to ffmpeg if it fails
	float temp_buffer[2304 * 2];
	drmp3_seek_to_pcm_frame(&decoder, 0);

	int successful_frames = 0;
	for (int i = 0; i < 2; i++)
	{
		if (drmp3_read_pcm_frames_f32(&decoder, 1152, &temp_buffer[0]) > 0)
			successful_frames++;
	}

	if (successful_frames < 2)
	{
		drmp3_uninit(&decoder);
		return FILE_LOAD_FAILED;
	}

	mBaseSamplerate = (float)decoder.sampleRate;
	mSampleCount = (unsigned int)samples;
	mFiletype = WAVSTREAM_DRMP3;

	delete[] mMp3SeekPoints;
	mMp3SeekPoints = nullptr;
	mMp3SeekPointCount = 0;

	double fileLengthSeconds = static_cast<double>(mSampleCount) / mBaseSamplerate;
	double seekPointIntervalSeconds = 5.0; // 1 seek point every 5 seconds

	// cap seek points to avoid memory blowup
	mMp3SeekPointCount = (drmp3_uint32)((fileLengthSeconds / seekPointIntervalSeconds) + 1);
	mMp3SeekPointCount = mMp3SeekPointCount < 16 ? 16 : mMp3SeekPointCount;
	mMp3SeekPointCount = mMp3SeekPointCount > 1000 ? 1000 : mMp3SeekPointCount;

	mMp3SeekPoints = new drmp3_seek_point[mMp3SeekPointCount];
	if (!drmp3_calculate_seek_points(&decoder, &mMp3SeekPointCount, mMp3SeekPoints))
	{
		delete[] mMp3SeekPoints;
		mMp3SeekPoints = nullptr;
		mMp3SeekPointCount = 0;
	}

	drmp3_uninit(&decoder);

	return SO_NO_ERROR;
}

result WavStream::loadffmpeg(File *fp)
{
	if (mSoloud)
		mSoloud->lockAudioMutex_internal();
	result retval = (!FFmpeg::FFmpegLoader::init() || !FFmpeg::FFmpegLoader::isAvailable()) ? FILE_LOAD_FAILED : SO_NO_ERROR;
	if (mSoloud)
		mSoloud->unlockAudioMutex_internal();
	if (retval != SO_NO_ERROR)
	{
#ifdef _DEBUG
		printf("debug: failed to load ffmpeg %s\n", FFmpeg::FFmpegLoader::getErrorDetails().c_str());
#endif
		return retval;
	}

	fp->seek(0);
	FFmpeg::FFmpegDecoder *decoder = FFmpeg::open(fp);
	if (!decoder)
		return FILE_LOAD_FAILED;

	mChannels = FFmpeg::getChannels(decoder);
	if (mChannels > MAX_CHANNELS)
	{
		mChannels = MAX_CHANNELS;
	}
	mBaseSamplerate = (float)FFmpeg::getSampleRate(decoder);
	mSampleCount = (unsigned int)FFmpeg::getTotalFrameCount(decoder);
	mFiletype = WAVSTREAM_FFMPEG;

	FFmpeg::close(decoder);
	return SO_NO_ERROR;
}

result WavStream::load(const char *aFilename)
{
	delete[] mFilename;
	delete mMemFile;
	delete[] mMp3SeekPoints;
	mMp3SeekPoints = nullptr;
	mMp3SeekPointCount = 0;
	mMemFile = nullptr;
	mFilename = nullptr;
	mSampleCount = 0;
	DiskFile fp;
	result res = fp.open(aFilename);
	if (res != SO_NO_ERROR)
		return res;

	int len = (int)strlen(aFilename);
	mFilename = new char[len + 1];
	memcpy(mFilename, aFilename, len);
	mFilename[len] = 0;

	res = parse(&fp);

	if (res != SO_NO_ERROR)
	{
		delete[] mFilename;
		mFilename = nullptr;
		return res;
	}

	return res;
}

result WavStream::loadMem(const unsigned char *aData, unsigned int aDataLen, bool aCopy, bool aTakeOwnership)
{
	delete[] mFilename;
	delete mMemFile;
	delete[] mMp3SeekPoints;
	mMp3SeekPoints = nullptr;
	mMp3SeekPointCount = 0;
	mStreamFile = nullptr;
	mMemFile = nullptr;
	mFilename = nullptr;
	mSampleCount = 0;

	if (aData == nullptr || aDataLen == 0)
		return INVALID_PARAMETER;

	MemoryFile *mf = new MemoryFile();
	int res = mf->openMem(aData, aDataLen, aCopy, aTakeOwnership);
	if (res != SO_NO_ERROR)
	{
		delete mf;
		return res;
	}

	res = parse(mf);

	if (res != SO_NO_ERROR)
	{
		delete mf;
		return res;
	}

	mMemFile = mf;

	return 0;
}

result WavStream::loadToMem(const char *aFilename)
{
	DiskFile df;
	int res = df.open(aFilename);
	if (res == SO_NO_ERROR)
	{
		res = loadFileToMem(&df);
	}
	return res;
}

result WavStream::loadFile(File *aFile)
{
	delete[] mFilename;
	delete mMemFile;
	delete[] mMp3SeekPoints;
	mMp3SeekPoints = nullptr;
	mMp3SeekPointCount = 0;
	mStreamFile = nullptr;
	mMemFile = nullptr;
	mFilename = nullptr;
	mSampleCount = 0;

	int res = parse(aFile);

	if (res != SO_NO_ERROR)
	{
		return res;
	}

	mStreamFile = aFile;

	return 0;
}

result WavStream::loadFileToMem(File *aFile)
{
	delete[] mFilename;
	delete mMemFile;
	delete[] mMp3SeekPoints;
	mMp3SeekPoints = nullptr;
	mMp3SeekPointCount = 0;
	mStreamFile = nullptr;
	mMemFile = nullptr;
	mFilename = nullptr;
	mSampleCount = 0;

	MemoryFile *mf = new MemoryFile();
	int res = mf->openFileToMem(aFile);
	if (res != SO_NO_ERROR)
	{
		delete mf;
		return res;
	}

	res = parse(mf);

	if (res != SO_NO_ERROR)
	{
		delete mf;
		return res;
	}

	mMemFile = mf;

	return res;
}

result WavStream::parse(File *aFile)
{
	if (mPreferFFmpeg && loadffmpeg(aFile) == SO_NO_ERROR)
	{
		return SO_NO_ERROR;
	}

	result res = FILE_LOAD_FAILED;
	int tag = aFile->read32();

	if (tag == MAKEDWORD('O', 'g', 'g', 'S'))
	{
		res = loadogg(aFile);
	}
	else if (tag == MAKEDWORD('R', 'I', 'F', 'F'))
	{
		res = loadwav(aFile);
	}
	else if (tag == MAKEDWORD('f', 'L', 'a', 'C'))
	{
		res = loadflac(aFile);
	}

	if (res != SO_NO_ERROR)
	{
		aFile->seek(0);
		if (loadmpg123(aFile) == SO_NO_ERROR)
		{
			res = SO_NO_ERROR;
		}
		else if (loaddrmp3(aFile) == SO_NO_ERROR)
		{
			res = SO_NO_ERROR;
		}
		else if (!mPreferFFmpeg && loadffmpeg(aFile) == SO_NO_ERROR)
		{
			res = SO_NO_ERROR;
		}
		else
		{
			res = FILE_LOAD_FAILED;
		}
	}
	return res;
}

AudioSourceInstance *WavStream::createInstance()
{
	return new WavStreamInstance(this);
}

double WavStream::getLength()
{
	if (mBaseSamplerate == 0)
		return 0;
	return static_cast<double>(mSampleCount) / mBaseSamplerate;
}
}; // namespace SoLoud
