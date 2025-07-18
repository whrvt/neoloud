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
#include "soloud_wav.h"

namespace SoLoud
{
WavInstance::WavInstance(Wav *aParent)
{
	mParent = aParent;
	mOffset = 0;
}

unsigned int WavInstance::getAudio(float *aBuffer, unsigned int aSamplesToRead, unsigned int aBufferSize)
{
	if (aBuffer == nullptr || mParent->mData == nullptr)
		return 0;

	unsigned int dataleft = mParent->mSampleCount - mOffset;
	unsigned int copylen = dataleft;
	if (copylen > aSamplesToRead)
		copylen = aSamplesToRead;

	unsigned int i;
	for (i = 0; i < mChannels; i++)
	{
		memcpy(aBuffer + i * aBufferSize, mParent->mData + mOffset + i * mParent->mSampleCount, sizeof(float) * copylen);
	}

	mOffset += copylen;
	return copylen;
}

result WavInstance::rewind()
{
	mOffset = 0;
	mStreamPosition = 0.0f;
	return 0;
}

bool WavInstance::hasEnded()
{
	if (!(mFlags & AudioSourceInstance::LOOPING) && mOffset >= mParent->mSampleCount)
	{
		return true;
	}
	return false;
}

Wav::Wav(bool preferFFmpeg)
{
	mData = nullptr;
	mSampleCount = 0;
	mPreferFFmpeg = preferFFmpeg;
}

Wav::~Wav()
{
	stop();
	delete[] mData;
}

#define MAKEDWORD(a, b, c, d) (((d) << 24) | ((c) << 16) | ((b) << 8) | (a))

result Wav::loadwav(MemoryFile *aReader)
{
	drwav decoder;

	if (!drwav_init_memory(&decoder, aReader->getMemPtr(), aReader->length(), nullptr))
	{
		return FILE_LOAD_FAILED;
	}

	drwav_uint64 samples = decoder.totalPCMFrameCount;

	if (!samples)
	{
		drwav_uninit(&decoder);
		return FILE_LOAD_FAILED;
	}

	mChannels = (unsigned int)decoder.channels;
	if (mChannels > MAX_CHANNELS)
		mChannels = MAX_CHANNELS;

	mBaseSamplerate = (float)decoder.sampleRate;
	mSampleCount = (unsigned int)samples;
	mData = new float[(size_t)(mSampleCount * mChannels)];
	memset(mData, 0, (size_t)(mSampleCount * mChannels) * sizeof(float));

	unsigned int i, j, k;
	for (i = 0; i < mSampleCount; i += SAMPLE_GRANULARITY)
	{
		float tmp[SAMPLE_GRANULARITY * MAX_CHANNELS];
		unsigned int blockSize = (mSampleCount - i) > SAMPLE_GRANULARITY ? SAMPLE_GRANULARITY : mSampleCount - i;
		drwav_read_pcm_frames_f32(&decoder, blockSize, tmp);
		for (j = 0; j < blockSize; j++)
		{
			for (k = 0; k < mChannels; k++)
			{
				mData[k * mSampleCount + i + j] = tmp[j * mChannels + k];
			}
		}
	}
	drwav_uninit(&decoder);

	return SO_NO_ERROR;
}

result Wav::loadogg(MemoryFile *aReader)
{
	int e = 0;
	stb_vorbis *vorbis = nullptr;
	vorbis = stb_vorbis_open_memory(aReader->getMemPtr(), aReader->length(), &e, nullptr);

	if (nullptr == vorbis)
	{
		return FILE_LOAD_FAILED;
	}

	stb_vorbis_info info = stb_vorbis_get_info(vorbis);
	unsigned int samples = stb_vorbis_stream_length_in_samples(vorbis);

	mChannels = (unsigned int)info.channels;
	if (mChannels > MAX_CHANNELS)
		mChannels = MAX_CHANNELS;

	mBaseSamplerate = (float)info.sample_rate;
	mSampleCount = samples;
	mData = new float[(size_t)(mSampleCount * mChannels)];
	memset(mData, 0, (size_t)(mSampleCount * mChannels) * sizeof(float));

	samples = 0;
	while (1)
	{
		float **outputs;
		int n = stb_vorbis_get_frame_float(vorbis, nullptr, &outputs);
		if (n == 0)
		{
			break;
		}

		unsigned int ch;
		for (ch = 0; ch < mChannels; ch++)
			memcpy(mData + samples + mSampleCount * ch, outputs[ch], sizeof(float) * n);

		samples += n;
	}
	stb_vorbis_close(vorbis);

	return 0;
}

result Wav::loadmpg123(MemoryFile *aReader)
{
	MPG123::MPG123Decoder *decoder = MPG123::open(aReader);

	if (!decoder)
		return FILE_LOAD_FAILED;

	int channels = MPG123::getChannels(decoder);
	int sampleRate = MPG123::getSampleRate(decoder);
	off_t totalFrames = MPG123::getTotalFrameCount(decoder);

	if (channels <= 0 || sampleRate <= 0 || totalFrames <= 0)
	{
		MPG123::close(decoder);
		return FILE_LOAD_FAILED;
	}

	mChannels = (unsigned int)channels;
	if (mChannels > MAX_CHANNELS)
		mChannels = MAX_CHANNELS;

	mBaseSamplerate = (float)sampleRate;
	mSampleCount = (unsigned int)totalFrames;
	mData = new float[(size_t)(mSampleCount * mChannels)];
	memset(mData, 0, (size_t)(mSampleCount * mChannels) * sizeof(float));

	unsigned int i, j, k;
	for (i = 0; i < mSampleCount; i += SAMPLE_GRANULARITY)
	{
		float tmp[SAMPLE_GRANULARITY * MAX_CHANNELS];
		unsigned int blockSize = (mSampleCount - i) > SAMPLE_GRANULARITY ? SAMPLE_GRANULARITY : mSampleCount - i;
		size_t framesRead = MPG123::readFrames(decoder, blockSize, tmp);

		if (framesRead == 0)
			break;

		for (j = 0; j < framesRead; j++)
		{
			for (k = 0; k < mChannels; k++)
			{
				mData[k * mSampleCount + i + j] = tmp[j * mChannels + k];
			}
		}

		if (framesRead < blockSize)
			break;
	}

	MPG123::close(decoder);
	return SO_NO_ERROR;
}

result Wav::loaddrmp3(MemoryFile *aReader)
{
	drmp3 decoder;

	if (!drmp3_init_memory(&decoder, aReader->getMemPtr(), aReader->length(), nullptr))
	{
		return FILE_LOAD_FAILED;
	}

	drmp3_uint64 samples = drmp3_get_pcm_frame_count(&decoder);

	if (!samples)
	{
		drmp3_uninit(&decoder);
		return FILE_LOAD_FAILED;
	}

	mChannels = (unsigned int)decoder.channels;
	if (mChannels > MAX_CHANNELS)
		mChannels = MAX_CHANNELS;

	mBaseSamplerate = (float)decoder.sampleRate;
	mSampleCount = (unsigned int)samples;
	mData = new float[(size_t)(mSampleCount * mChannels)];
	memset(mData, 0, (size_t)(mSampleCount * mChannels) * sizeof(float));

	drmp3_seek_to_pcm_frame(&decoder, 0);

	unsigned int i, j, k;
	for (i = 0; i < mSampleCount; i += SAMPLE_GRANULARITY)
	{
		float tmp[SAMPLE_GRANULARITY * MAX_CHANNELS];
		unsigned int blockSize = (mSampleCount - i) > SAMPLE_GRANULARITY ? SAMPLE_GRANULARITY : mSampleCount - i;
		drmp3_read_pcm_frames_f32(&decoder, blockSize, tmp);
		for (j = 0; j < blockSize; j++)
		{
			for (k = 0; k < mChannels; k++)
			{
				mData[k * mSampleCount + i + j] = tmp[j * mChannels + k];
			}
		}
	}
	drmp3_uninit(&decoder);

	return SO_NO_ERROR;
}

result Wav::loadflac(MemoryFile *aReader)
{
	drflac *decoder = drflac_open_memory(aReader->mDataPtr, aReader->mDataLength, nullptr);

	if (!decoder)
	{
		return FILE_LOAD_FAILED;
	}

	drflac_uint64 samples = decoder->totalPCMFrameCount;

	if (!samples)
	{
		drflac_close(decoder);
		return FILE_LOAD_FAILED;
	}

	mChannels = (unsigned int)decoder->channels;
	if (mChannels > MAX_CHANNELS)
		mChannels = MAX_CHANNELS;

	mBaseSamplerate = (float)decoder->sampleRate;
	mSampleCount = (unsigned int)samples;
	mData = new float[(size_t)(mSampleCount * mChannels)];
	memset(mData, 0, (size_t)(mSampleCount * mChannels) * sizeof(float));

	drflac_seek_to_pcm_frame(decoder, 0);

	unsigned int i, j, k;
	for (i = 0; i < mSampleCount; i += SAMPLE_GRANULARITY)
	{
		float tmp[SAMPLE_GRANULARITY * MAX_CHANNELS];
		unsigned int blockSize = (mSampleCount - i) > SAMPLE_GRANULARITY ? SAMPLE_GRANULARITY : mSampleCount - i;
		drflac_read_pcm_frames_f32(decoder, blockSize, tmp);
		for (j = 0; j < blockSize; j++)
		{
			for (k = 0; k < mChannels; k++)
			{
				mData[k * mSampleCount + i + j] = tmp[j * mChannels + k];
			}
		}
	}
	drflac_close(decoder);

	return SO_NO_ERROR;
}

result Wav::loadffmpeg(MemoryFile *aReader)
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

	float *data = nullptr;
	unsigned int channels = 2, sampleCount = 0;
	float sampleRate = 44100.0f;

	result res = FFmpeg::loadToMemory(aReader, &data, &sampleCount, &channels, &sampleRate);
	if (res != SO_NO_ERROR || sampleCount == 0)
	{
		if (data)
			delete[] data;
		return FILE_LOAD_FAILED;
	}

	mData = data;
	mSampleCount = sampleCount;
	mChannels = channels;
	mBaseSamplerate = sampleRate;

	return SO_NO_ERROR;
}

result Wav::testAndLoadFile(MemoryFile *aReader)
{
	delete[] mData;
	mData = nullptr;
	mSampleCount = 0;
	mChannels = 1;

	if (mPreferFFmpeg && loadffmpeg(aReader) == SO_NO_ERROR)
	{
		return SO_NO_ERROR;
	}

	result res = FILE_LOAD_FAILED;

	int tag = aReader->read32();
	if (tag == MAKEDWORD('O', 'g', 'g', 'S'))
	{
		res = loadogg(aReader);
	}
	else if (tag == MAKEDWORD('R', 'I', 'F', 'F'))
	{
		res = loadwav(aReader);
	}
	else if (tag == MAKEDWORD('f', 'L', 'a', 'C'))
	{
		res = loadflac(aReader);
	}

	if (res != SO_NO_ERROR)
	{
		aReader->seek(0);
		if (loadmpg123(aReader) == SO_NO_ERROR)
		{
			res = SO_NO_ERROR;
		}
		else if (loaddrmp3(aReader) == SO_NO_ERROR)
		{
			res = SO_NO_ERROR;
		}
		else if (!mPreferFFmpeg && loadffmpeg(aReader) == SO_NO_ERROR)
		{
			res = SO_NO_ERROR;
		}
	}
	return res;
}

result Wav::load(const char *aFilename)
{
	if (aFilename == nullptr)
		return INVALID_PARAMETER;
	stop();
	DiskFile dr;
	int res = dr.open(aFilename);
	if (res == SO_NO_ERROR)
		return loadFile(&dr);
	return res;
}

result Wav::loadMem(const unsigned char *aMem, unsigned int aLength, bool aCopy, bool aTakeOwnership)
{
	if (aMem == nullptr || aLength == 0)
		return INVALID_PARAMETER;
	stop();

	MemoryFile dr;
	dr.openMem(aMem, aLength, aCopy, aTakeOwnership);
	return testAndLoadFile(&dr);
}

result Wav::loadFile(File *aFile)
{
	if (!aFile)
		return INVALID_PARAMETER;
	stop();

	MemoryFile mr;
	result res = mr.openFileToMem(aFile);

	if (res != SO_NO_ERROR)
	{
		return res;
	}
	return testAndLoadFile(&mr);
}

AudioSourceInstance *Wav::createInstance()
{
	return new WavInstance(this);
}

double Wav::getLength()
{
	if (mBaseSamplerate == 0)
		return 0;
	return mSampleCount / mBaseSamplerate;
}

result Wav::loadRawWave8(unsigned char *aMem, unsigned int aLength, float aSamplerate, unsigned int aChannels)
{
	if (aMem == nullptr || aLength == 0 || aSamplerate <= 0 || aChannels < 1)
		return INVALID_PARAMETER;
	stop();
	delete[] mData;
	mData = new float[aLength];
	mSampleCount = aLength / aChannels;
	mChannels = aChannels;
	mBaseSamplerate = aSamplerate;
	unsigned int i;
	for (i = 0; i < aLength; i++)
		mData[i] = ((signed)aMem[i] - 128) / (float)0x80;
	return SO_NO_ERROR;
}

result Wav::loadRawWave16(short *aMem, unsigned int aLength, float aSamplerate, unsigned int aChannels)
{
	if (aMem == nullptr || aLength == 0 || aSamplerate <= 0 || aChannels < 1)
		return INVALID_PARAMETER;
	stop();
	delete[] mData;
	mData = new float[aLength];
	mSampleCount = aLength / aChannels;
	mChannels = aChannels;
	mBaseSamplerate = aSamplerate;
	unsigned int i;
	for (i = 0; i < aLength; i++)
		mData[i] = ((signed short)aMem[i]) / (float)0x8000;
	return SO_NO_ERROR;
}

result Wav::loadRawWave(float *aMem, unsigned int aLength, float aSamplerate, unsigned int aChannels, bool aCopy, bool aTakeOwndership)
{
	if (aMem == nullptr || aLength == 0 || aSamplerate <= 0 || aChannels < 1)
		return INVALID_PARAMETER;
	stop();
	delete[] mData;
	if (aCopy == true || aTakeOwndership == false)
	{
		mData = new float[aLength];
		memcpy(mData, aMem, sizeof(float) * aLength);
	}
	else
	{
		mData = aMem;
	}
	mSampleCount = aLength / aChannels;
	mChannels = aChannels;
	mBaseSamplerate = aSamplerate;
	return SO_NO_ERROR;
}
}; // namespace SoLoud
