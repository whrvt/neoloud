/*
SoLoud audio engine - mpg123 interface
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

#ifdef WITH_LIBMPG123

#include "soloud_file.h"

#include "soloud_mpg123.h"

#include <cstdlib>
#include <cstring>
#include <memory> // for unique_ptr

#include <mpg123.h>

namespace SoLoud::MPG123
{

struct MPG123Decoder
{
	struct MPG123HandleDeleter
	{
		void operator()(mpg123_handle *handle) const
		{
			if (!handle)
				return;
			mpg123_close(handle);
			mpg123_delete(handle);
		}
	};
	struct MPG123Handle : public std::unique_ptr<mpg123_handle, MPG123HandleDeleter>
	{
		operator mpg123_handle *() const { return this->get(); }
	};

	struct CFree
	{
		inline void operator()(void *p) const noexcept { free(p); }
	};
	struct TempBuffer : public std::unique_ptr<unsigned char, CFree>
	{
		operator unsigned char *() const { return this->get(); }
	};

	MPG123Handle handle;
	TempBuffer tempBuffer;
	size_t tempBufferSize;
	long rate;
	off_t totalFrames;
	int channels;
	bool ended;
};

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

	auto decoder = std::make_unique<MPG123Decoder>();
	memset((void *)decoder.get(), 0, sizeof(MPG123Decoder));

	decoder->handle.reset(mpg123_new(nullptr, nullptr));

	if (!decoder->handle)
	{
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
	// libmpg123 does not allow setting a custom log output function unfortunately
	mpg123_param(decoder->handle, MPG123_ADD_FLAGS, MPG123_QUIET, 0);
#else
	// (noisy)
	// mpg123_param(decoder->handle, MPG123_REMOVE_FLAGS, MPG123_QUIET, 0);
	// mpg123_param(decoder->handle, MPG123_VERBOSE, 3, 0);
#endif

	// set up custom I/O
	if (mpg123_replace_reader_handle(decoder->handle, readCallback, seekCallback, nullptr) != MPG123_OK)
	{
		return nullptr;
	}

	aFile->seek(0);

	// open with custom file handle
	if (mpg123_open_handle(decoder->handle, aFile) != MPG123_OK)
	{
		return nullptr;
	}

	// get format information
	long rate = 0;
	int channels = 0, encoding = 0;
	if (mpg123_getformat(decoder->handle, &rate, &channels, &encoding) != MPG123_OK)
	{
		return nullptr;
	}

	decoder->rate = rate;
	decoder->channels = channels;

	// validate that this is actually MPEG audio that mpg123 can decode
	struct mpg123_frameinfo2 frameInfo{};
	if (mpg123_info2(decoder->handle, &frameInfo) != MPG123_OK ||
	    (frameInfo.layer < 1 || frameInfo.layer > 3 || frameInfo.version < 0 || frameInfo.version > 2 || frameInfo.rate <= 0 || frameInfo.bitrate <= 0))
	{
		return nullptr;
	}

	// get total frame count
	decoder->totalFrames = mpg123_length(decoder->handle);

	// When mpg123 parses Xing/LAME/VBRI headers, it sets vbr to MPG123_VBR or MPG123_ABR.
	// If vbr == MPG123_CBR, either the file is genuine CBR, or it's VBR without proper
	// headers. In the latter case, mpg123_length estimates duration from first-frame
	// bitrate which can be wildly inaccurate. Sample a few frames to detect this case.
	bool needsFullScan = decoder->totalFrames < MPG123_OK;

	if (!needsFullScan && frameInfo.vbr == MPG123_CBR)
	{
		off_t pos = mpg123_tell(decoder->handle);
		int firstBitrate = frameInfo.bitrate;

		int i = 0;
		for (; i < 20; i++)
		{
			if (mpg123_framebyframe_next(decoder->handle) != MPG123_OK)
				break;

			struct mpg123_frameinfo2 sampleInfo{};
			if (mpg123_info2(decoder->handle, &sampleInfo) == MPG123_OK && sampleInfo.bitrate != firstBitrate)
			{
				needsFullScan = true;
				break;
			}
		}

		// Even if we don't need a full scan, update totalFrames after decoding a few since
		// the estimate can be more accurate after doing so.
		if (!needsFullScan && i == 20)
			decoder->totalFrames = mpg123_length(decoder->handle);

		mpg123_seek(decoder->handle, pos, SEEK_SET);
	}

	if (needsFullScan)
	{
		off_t pos = mpg123_tell(decoder->handle);
		if (mpg123_scan(decoder->handle) == MPG123_OK)
		{
			const off_t newLength = mpg123_length(decoder->handle);
#ifdef _DEBUG
			if (newLength != decoder->totalFrames)
			{
				SoLoud::logStdout("SoLoud::MPG123::open: did a full scan and got length %ld (was %ld)\n", newLength, decoder->totalFrames);
			}
#endif
			decoder->totalFrames = newLength;
		}
		else if (decoder->totalFrames < MPG123_OK)
		{
			decoder->totalFrames = 0;
		}
		mpg123_seek(decoder->handle, pos, SEEK_SET);
	}

	if (decoder->totalFrames <= 0)
		return nullptr;

	return decoder.release();
}

void close(MPG123Decoder *aDecoder)
{
	if (!aDecoder)
		return;

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
		unsigned char *tmpAlloc = nullptr;
		unsigned char *old = aDecoder->tempBuffer.release();

		if (!(tmpAlloc = static_cast<unsigned char *>(realloc(old, requiredBytes))))
		{
			aDecoder->tempBuffer.reset(old);
			return 0; // we are in trouble...
		}

		aDecoder->tempBuffer.reset(tmpAlloc);
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
