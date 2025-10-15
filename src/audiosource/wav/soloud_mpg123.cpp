/*
SoLoud audio engine
Copyright (c) 2013-2018 Jari Komppa
Copyright (c) 2025 William Horvath (mpg123 interface)

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

#ifdef WITH_LIBMPG123

#include "soloud_file.h"

#include "soloud_mpg123.h"

#include <cstring>
#include <mpg123.h>

namespace SoLoud::MPG123
{

namespace // static
{

mpg123_ssize_t readCallback(void *handle, void *buf, size_t count)
{
	File *file = static_cast<File *>(handle);
	return static_cast<mpg123_ssize_t>(file->read(static_cast<unsigned char *>(buf), static_cast<unsigned int>(count)));
}

off_t seekCallback(void *handle, off_t offset, int whence)
{
	File *file = static_cast<File *>(handle);

	switch (whence)
	{
	case SEEK_SET:
		file->seek(static_cast<int>(offset));
		break;
	case SEEK_CUR:
		file->seek(static_cast<int>(file->pos() + offset));
		break;
	case SEEK_END:
		file->seek(static_cast<int>(file->length() + offset));
		break;
	default:
		return -1;
	}

	return static_cast<off_t>(file->pos());
}
} // namespace

MPG123Decoder *open(File *aFile)
{
	if (!aFile)
		return nullptr;

	auto *decoder = new MPG123Decoder();
	memset(decoder, 0, sizeof(MPG123Decoder));

	decoder->file = aFile;
	decoder->handle = mpg123_new(nullptr, nullptr);

	if (!decoder->handle)
	{
		delete decoder;
		return nullptr;
	}

	// configure some opts
	// 64kb leniency, allow some slightly messed up mp3s but don't make it unlimited (otherwise we'll accept complete garbage as valid mp3 data)
	mpg123_param(decoder->handle, MPG123_RESYNC_LIMIT, 65536, 0);
	mpg123_param(decoder->handle, MPG123_REMOVE_FLAGS, MPG123_AUTO_RESAMPLE, 0); // soloud handles this already
	mpg123_param(decoder->handle, MPG123_ADD_FLAGS, MPG123_GAPLESS, 0);          // this is default, but still we want it
	mpg123_param(decoder->handle, MPG123_ADD_FLAGS, MPG123_NO_FRANKENSTEIN, 0);  // no stitched-together mp3s
	mpg123_param(decoder->handle, MPG123_ADD_FLAGS, MPG123_SEEKBUFFER, 0);       // more accurate seeking
	mpg123_param(decoder->handle, MPG123_ADD_FLAGS, MPG123_SKIP_ID3V2, 0);       // we don't care
	mpg123_param(decoder->handle, MPG123_ADD_FLAGS, MPG123_FORCE_FLOAT, 0);      // important! it makes life a lot easier if we always get float32 from libmpg123
	mpg123_param(decoder->handle, MPG123_PREFRAMES, 1, 0);                       // layer 3 prefill for seeking
#ifndef _DEBUG
	mpg123_param(decoder->handle, MPG123_ADD_FLAGS, MPG123_QUIET, 0);
#endif

	// set up custom I/O
	if (mpg123_replace_reader_handle(decoder->handle, readCallback, seekCallback, nullptr) != MPG123_OK)
	{
		mpg123_delete(decoder->handle);
		delete decoder;
		return nullptr;
	}

	aFile->seek(0);

	// open with custom handle
	if (mpg123_open_handle(decoder->handle, aFile) != MPG123_OK)
	{
		mpg123_delete(decoder->handle);
		delete decoder;
		return nullptr;
	}

	// get format information
	long rate = 0;
	int channels = 0, encoding = 0;
	if (mpg123_getformat(decoder->handle, &rate, &channels, &encoding) != MPG123_OK)
	{
		mpg123_close(decoder->handle);
		mpg123_delete(decoder->handle);
		delete decoder;
		return nullptr;
	}

	decoder->rate = rate;
	decoder->channels = channels;

	// get total frame count
	decoder->totalFrames = mpg123_length(decoder->handle);
	bool didFullScan = false;

	if (decoder->totalFrames < MPG123_OK)
	{
		if (mpg123_scan(decoder->handle) == MPG123_OK)
		{
			decoder->totalFrames = mpg123_length(decoder->handle);
			didFullScan = true;
		}
		else
		{
			decoder->totalFrames = 0;
		}
	}

	// probe end of file to detect trailing metadata that affects length
	// if it keeps increasing after 3 tries just give up and do a full mpg123_scan
	if (!didFullScan && decoder->totalFrames > 0)
	{
		int samplesPerMpegFrame = mpg123_spf(decoder->handle);
		bool needsFullScan = false;

		// try up to 3 probe iterations
		for (int probeIteration = 0; probeIteration < 3; probeIteration++)
		{
			off_t seekTarget = decoder->totalFrames - (samplesPerMpegFrame * 3LL);

			if (seekTarget <= 0 || mpg123_seek(decoder->handle, seekTarget, SEEK_SET) < 0)
				break;

			off_t probeBeginLength = decoder->totalFrames;
			unsigned char *audio = nullptr;
			size_t bytes = 0;
			int64_t frameNum = 0;

			// decode up to 20 frames
			for (int i = 0; i < 20; i++)
			{
				int result = mpg123_decode_frame64(decoder->handle, &frameNum, &audio, &bytes);

				if (result == MPG123_DONE || result == MPG123_ERR)
					break;
			}

			off_t probeEndLength = mpg123_length(decoder->handle);
			decoder->totalFrames = probeEndLength;

			// check if length changed
			if (probeBeginLength != probeEndLength)
			{
				needsFullScan = true;
			}
			else
			{
				// length stable, we're done
				needsFullScan = false;
				break;
			}
		}

		// if length kept changing, do full scan
		if (needsFullScan)
		{
			if (mpg123_scan(decoder->handle) == MPG123_OK)
			{
				decoder->totalFrames = mpg123_length(decoder->handle);
			}
		}

		// seek back to start after probing
		mpg123_seek(decoder->handle, 0, SEEK_SET);
	}

	// validate that this is actually MPEG audio that mpg123 can decode
	struct mpg123_frameinfo2 frameInfo{};
	if (decoder->totalFrames <= 0 || mpg123_info2(decoder->handle, &frameInfo) != MPG123_OK ||
	    (frameInfo.layer < 1 || frameInfo.layer > 3 || frameInfo.version < 0 || frameInfo.version > 2 || frameInfo.rate <= 0 || frameInfo.bitrate <= 0))
	{
		// frame count couldn't be determined, or no frame info, or layer isn't between 1 and 3
		mpg123_close(decoder->handle);
		mpg123_delete(decoder->handle);
		delete decoder;
		return nullptr;
	}

	return decoder;
}

void close(MPG123Decoder *aDecoder)
{
	if (!aDecoder)
		return;

	if (aDecoder->handle)
	{
		mpg123_close(aDecoder->handle);
		mpg123_delete(aDecoder->handle);
	}

	delete[] aDecoder->tempBuffer;
	delete aDecoder;
}

int getChannels(MPG123Decoder *aDecoder)
{
	return aDecoder ? aDecoder->channels : 0;
}

int getSampleRate(MPG123Decoder *aDecoder)
{
	return aDecoder ? static_cast<int>(aDecoder->rate) : 0;
}

off_t getTotalFrameCount(MPG123Decoder *aDecoder)
{
	return aDecoder ? aDecoder->totalFrames : 0;
}

size_t readFrames(MPG123Decoder *aDecoder, size_t aFrameCount, float *aBuffer)
{
	if (!aDecoder || !aDecoder->handle || !aBuffer)
		return 0;

	size_t done = 0;
	size_t requiredBytes = aFrameCount * aDecoder->channels * sizeof(float);

	// reuse buffer, resize only when needed
	if (aDecoder->tempBufferSize < requiredBytes)
	{
		delete[] aDecoder->tempBuffer;
		aDecoder->tempBuffer = new unsigned char[requiredBytes];
		aDecoder->tempBufferSize = requiredBytes;
	}

	// read raw data from mpg123 (always float32 due to MPG123_FORCE_FLOAT)
	int result = mpg123_read(aDecoder->handle, aDecoder->tempBuffer, requiredBytes, &done);

	if (result == MPG123_NEW_FORMAT)
	{
		// format changed mid-stream, update decoder info
		long rate = 0;
		int channels = 0, encoding = 0;
		if (mpg123_getformat(aDecoder->handle, &rate, &channels, &encoding) == MPG123_OK)
		{
			aDecoder->rate = rate;
			aDecoder->channels = channels;
		}
		return 0; // caller should try again
	}

	if (result == MPG123_DONE)
		aDecoder->ended = true;
	else
		aDecoder->ended = false;

	if (result != MPG123_OK && result != MPG123_DONE)
		return 0;

	// we're forcing float32, so just copy
	size_t samplesRead = done / sizeof(float);
	size_t framesRead = samplesRead / aDecoder->channels;
	memcpy(aBuffer, aDecoder->tempBuffer, done);

	return framesRead;
}

off_t seekToFrame(MPG123Decoder *aDecoder, off_t aFrame)
{
	if (!aDecoder || !aDecoder->handle)
		return -1;

	off_t result = mpg123_seek(aDecoder->handle, aFrame, SEEK_SET);

	if (result == MPG123_DONE)
		aDecoder->ended = true;
	else
		aDecoder->ended = false;

	return result >= 0 ? result : aDecoder->totalFrames;
}

off_t getCurrentFrame(MPG123Decoder *aDecoder)
{
	if (!aDecoder || !aDecoder->handle)
		return 0;

	off_t result = mpg123_tell(aDecoder->handle);

	if (result == MPG123_DONE)
		aDecoder->ended = true;
	else
		aDecoder->ended = false;

	return result >= 0 ? result : aDecoder->totalFrames;
}

bool isAtEnd(MPG123Decoder *aDecoder)
{
	if (!aDecoder || !aDecoder->handle)
		return true;

	return aDecoder->ended;
}
} // namespace SoLoud::MPG123

#endif // WITH_LIBMPG123
