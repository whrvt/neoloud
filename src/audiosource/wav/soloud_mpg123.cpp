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

#include "soloud.h"
#include "soloud_file.h"

#include "soloud_mpg123.h"

#include <cstring>
#include <mpg123.h>
#include <mutex>

namespace SoLoud::MPG123
{
static std::once_flag gMpg123InitFlag;
static std::mutex gMpg123DeinitMutex;
static bool gMpg123Initialized = false;

bool init()
{
	bool success = true;
	std::call_once(gMpg123InitFlag, [&success]() {
		if (mpg123_init() != MPG123_OK)
		{
			success = false;
			return;
		}
		gMpg123Initialized = true;
	});

	return success && gMpg123Initialized;
}

void deinit()
{
	std::lock_guard<std::mutex> lock(gMpg123DeinitMutex);
	if (gMpg123Initialized)
	{
		mpg123_exit();
		gMpg123Initialized = false;
	}
}

ssize_t readCallback(void *handle, void *buf, size_t count)
{
	File *file = static_cast<File *>(handle);
	return static_cast<ssize_t>(file->read(static_cast<unsigned char *>(buf), static_cast<unsigned int>(count)));
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

MPG123Decoder *open(File *aFile)
{
	if (!init())
		return nullptr;

	if (!aFile)
		return nullptr;

	MPG123Decoder *decoder = new MPG123Decoder();
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
	long rate;
	int channels, encoding;
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
	if (decoder->totalFrames == MPG123_ERR)
	{
		// try to get frame length info before expensive scan
		off_t frameLength = mpg123_framelength(decoder->handle);
		if (frameLength > 0)
		{
			decoder->totalFrames = frameLength;
		}
		else
		{
			// fallback: scan the file to get accurate length
			off_t pos = mpg123_tell(decoder->handle);
			if (mpg123_scan(decoder->handle) == MPG123_OK)
			{
				decoder->totalFrames = mpg123_length(decoder->handle);
				mpg123_seek(decoder->handle, pos, SEEK_SET);
			}
			else
			{
				decoder->totalFrames = 0;
			}
		}
	}

	// validate that this is actually MPEG audio that mpg123 can decode
	struct mpg123_frameinfo2 frameInfo;
	if (decoder->totalFrames <= 0 || mpg123_info2(decoder->handle, &frameInfo) != MPG123_OK || (frameInfo.layer < 1 || frameInfo.layer > 3))
	{
		// frame count couldn't be determined, or no frame info, or layer isn't between 1 and 3
		mpg123_close(decoder->handle);
		mpg123_delete(decoder->handle);
		delete decoder;
		return nullptr;
	}

	decoder->initialized = true;
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
		long rate;
		int channels, encoding;
		if (mpg123_getformat(aDecoder->handle, &rate, &channels, &encoding) == MPG123_OK)
		{
			aDecoder->rate = rate;
			aDecoder->channels = channels;
		}
		return 0; // caller should try again
	}

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

	return mpg123_seek(aDecoder->handle, aFrame, SEEK_SET);
}

off_t getCurrentFrame(MPG123Decoder *aDecoder)
{
	if (!aDecoder || !aDecoder->handle)
		return 0;

	return mpg123_tell(aDecoder->handle);
}
} // namespace SoLoud::MPG123
