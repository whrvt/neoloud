/*
SoLoud audio engine
Copyright (c) 2013-2020 Jari Komppa

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

#include "soloud_fft.h"
#include "soloud_internal.h"
#include "soloud_thread.h"
#include <float.h> // _controlfp
#include <math.h>  // sin
#include <stdlib.h>
#include <string.h>

#ifdef SOLOUD_SSE_INTRINSICS
#include <xmmintrin.h>

#include <cstddef>
#ifdef _M_IX86
#include <emmintrin.h>
#endif
#endif

// #define FLOATING_POINT_DEBUG

namespace SoLoud
{
AlignedFloatBuffer::AlignedFloatBuffer()
{
	mBasePtr = 0;
	mData = 0;
	mFloats = 0;
}

result AlignedFloatBuffer::init(unsigned int aFloats)
{
	delete[] mBasePtr;
	mBasePtr = 0;
	mData = 0;
	mFloats = aFloats;
#ifndef SOLOUD_SSE_INTRINSICS
	mBasePtr = new unsigned char[aFloats * sizeof(float)];
	if (mBasePtr == NULL)
		return OUT_OF_MEMORY;
	mData = (float *)mBasePtr;
#else
	mBasePtr = new unsigned char[aFloats * sizeof(float) + 16];
	if (mBasePtr == NULL)
		return OUT_OF_MEMORY;
	mData = (float *)(((size_t)mBasePtr + 15) & ~15);
#endif
	return SO_NO_ERROR;
}

void AlignedFloatBuffer::clear()
{
	memset(mData, 0, sizeof(float) * mFloats);
}

AlignedFloatBuffer::~AlignedFloatBuffer()
{
	delete[] mBasePtr;
}

TinyAlignedFloatBuffer::TinyAlignedFloatBuffer()
{
	unsigned char *basePtr = &mActualData[0];
	mData = (float *)(((size_t)basePtr + 15) & ~15);
}

Soloud::Soloud()
{
#ifdef FLOATING_POINT_DEBUG
	unsigned int u;
	u = _controlfp(0, 0);
	u = u & ~(_EM_INVALID | /*_EM_DENORMAL |*/ _EM_ZERODIVIDE | _EM_OVERFLOW /*| _EM_UNDERFLOW  | _EM_INEXACT*/);
	_controlfp(u, _MCW_EM);
#endif
	mResampler = SOLOUD_DEFAULT_RESAMPLER;
	mInsideAudioThreadMutex = false;
	mScratchSize = 0;
	mSamplerate = 0;
	mBufferSize = 0;
	mFlags = 0;
	mGlobalVolume = 0;
	mPlayIndex = 0;
	mBackendData = NULL;
	mAudioThreadMutex = NULL;
	mPostClipScaler = 0;
	mBackendCleanupFunc = NULL;
	mBackendPauseFunc = NULL;
	mBackendResumeFunc = NULL;
	mChannels = 2;
	mStreamTime = 0;
	mLastClockedTime = 0;
	mAudioSourceID = 1;
	mBackendString = 0;
	mBackendID = 0;
	mActiveVoiceDirty = true;
	mActiveVoiceCount = 0;
	int i;
	for (i = 0; i < VOICE_COUNT; i++)
		mActiveVoice[i] = 0;
	for (i = 0; i < FILTERS_PER_STREAM; i++)
	{
		mFilter[i] = NULL;
		mFilterInstance[i] = NULL;
	}
	for (i = 0; i < 256; i++)
	{
		mFFTData[i] = 0;
		mVisualizationWaveData[i] = 0;
		mWaveData[i] = 0;
	}
	for (i = 0; i < MAX_CHANNELS; i++)
	{
		mVisualizationChannelVolume[i] = 0;
	}
	for (i = 0; i < VOICE_COUNT; i++)
	{
		mVoice[i] = 0;
	}
	mVoiceGroup = 0;
	mVoiceGroupCount = 0;

	m3dPosition[0] = 0;
	m3dPosition[1] = 0;
	m3dPosition[2] = 0;
	m3dAt[0] = 0;
	m3dAt[1] = 0;
	m3dAt[2] = -1;
	m3dUp[0] = 0;
	m3dUp[1] = 1;
	m3dUp[2] = 0;
	m3dVelocity[0] = 0;
	m3dVelocity[1] = 0;
	m3dVelocity[2] = 0;
	m3dSoundSpeed = 343.3f;
	mMaxActiveVoices = 16;
	mHighestVoice = 0;
	for (i = 0; i < 3 * MAX_CHANNELS; i++)
		m3dSpeakerPosition[i] = 0;
}

Soloud::~Soloud()
{
	// let's stop all sounds before deinit, so we don't mess up our mutexes
	stopAll();
	deinit();
	unsigned int i;
	for (i = 0; i < FILTERS_PER_STREAM; i++)
	{
		delete mFilterInstance[i];
	}
	for (i = 0; i < mVoiceGroupCount; i++)
		delete[] mVoiceGroup[i];
	delete[] mVoiceGroup;
}

void Soloud::deinit()
{
	// Make sure no audio operation is currently pending
	lockAudioMutex_internal();
	unlockAudioMutex_internal();
	SOLOUD_ASSERT(!mInsideAudioThreadMutex);
	stopAll();
	if (mBackendCleanupFunc)
		mBackendCleanupFunc(this);
	mBackendCleanupFunc = 0;
	if (mAudioThreadMutex)
		Thread::destroyMutex(mAudioThreadMutex);
	mAudioThreadMutex = NULL;
}

result Soloud::init(unsigned int aFlags, unsigned int aBackend, unsigned int aSamplerate, unsigned int aBufferSize, unsigned int aChannels)
{
	if (aBackend >= BACKEND_MAX || aChannels == 3 || aChannels == 5 || aChannels == 7 || aChannels > MAX_CHANNELS)
		return INVALID_PARAMETER;

	deinit();

	mAudioThreadMutex = Thread::createMutex();

	mBackendID = 0;
	mBackendString = 0;

	int samplerate = 44100;
	int buffersize = 2048;
	int inited = 0;

	if (aSamplerate != Soloud::AUTO)
		samplerate = aSamplerate;
	else
		samplerate = Soloud::AUTO;

	if (aBufferSize != Soloud::AUTO)
		buffersize = aBufferSize;
	else
		buffersize = Soloud::AUTO;

	if (!inited && (aBackend == Soloud::MINIAUDIO || aBackend == Soloud::AUTO))
	{
		int ret = miniaudio_init(this, aFlags, samplerate, buffersize, aChannels);
		if (ret == 0)
		{
			inited = 1;
			mBackendID = Soloud::MINIAUDIO;
		}

		if (ret != 0 && aBackend != Soloud::AUTO)
			return ret;
	}

#if defined(WITH_SDL3)
	if (!inited && (aBackend == Soloud::SDL3 || aBackend == Soloud::AUTO))
	{
		int ret = sdl3_init(this, aFlags, samplerate, buffersize, aChannels);
		if (ret == 0)
		{
			inited = 1;
			mBackendID = Soloud::SDL3;
		}

		if (ret != 0 && aBackend != Soloud::AUTO)
			return ret;
	}
#endif

	if (!inited && (aBackend == Soloud::NOSOUND || aBackend == Soloud::AUTO))
	{
		if (aBufferSize == Soloud::AUTO)
			buffersize = 2048;

		int ret = nosound_init(this, aFlags, samplerate, buffersize, aChannels);
		if (ret == 0)
		{
			inited = 1;
			mBackendID = Soloud::NOSOUND;
		}

		if (ret != 0 && aBackend != Soloud::AUTO)
			return ret;
	}

	if (!inited && (aBackend == Soloud::NULLDRIVER))
	{
		if (aBufferSize == Soloud::AUTO)
			buffersize = 2048;

		int ret = null_init(this, aFlags, samplerate, buffersize, aChannels);
		if (ret == 0)
		{
			inited = 1;
			mBackendID = Soloud::NULLDRIVER;
		}

		if (ret != 0)
			return ret;
	}

	if (!inited && aBackend != Soloud::AUTO)
		return NOT_IMPLEMENTED;
	if (!inited)
		return UNKNOWN_ERROR;
	return 0;
}

result Soloud::pause()
{
	if (mBackendPauseFunc)
		return mBackendPauseFunc(this);

	return NOT_IMPLEMENTED;
}

result Soloud::resume()
{
	if (mBackendResumeFunc)
		return mBackendResumeFunc(this);

	return NOT_IMPLEMENTED;
}

void Soloud::postinit_internal(unsigned int aSamplerate, unsigned int aBufferSize, unsigned int aFlags, unsigned int aChannels)
{
	mGlobalVolume = 1;
	mChannels = aChannels;
	mSamplerate = aSamplerate;
	mBufferSize = aBufferSize;
	mScratchSize = (aBufferSize + 15) & (~0xf); // round to the next div by 16
	if (mScratchSize < SAMPLE_GRANULARITY * 2)
		mScratchSize = SAMPLE_GRANULARITY * 2;
	if (mScratchSize < 4096)
		mScratchSize = 4096;
	mScratch.init(mScratchSize * MAX_CHANNELS);
	mOutputScratch.init(mScratchSize * MAX_CHANNELS);

	mFlags = aFlags;
	mPostClipScaler = 0.95f;
	switch (mChannels)
	{
	case 1:
		m3dSpeakerPosition[0 * 3 + 0] = 0;
		m3dSpeakerPosition[0 * 3 + 1] = 0;
		m3dSpeakerPosition[0 * 3 + 2] = 1;
		break;
	case 2:
		m3dSpeakerPosition[0 * 3 + 0] = 2;
		m3dSpeakerPosition[0 * 3 + 1] = 0;
		m3dSpeakerPosition[0 * 3 + 2] = 1;
		m3dSpeakerPosition[1 * 3 + 0] = -2;
		m3dSpeakerPosition[1 * 3 + 1] = 0;
		m3dSpeakerPosition[1 * 3 + 2] = 1;
		break;
	case 4:
		m3dSpeakerPosition[0 * 3 + 0] = 2;
		m3dSpeakerPosition[0 * 3 + 1] = 0;
		m3dSpeakerPosition[0 * 3 + 2] = 1;
		m3dSpeakerPosition[1 * 3 + 0] = -2;
		m3dSpeakerPosition[1 * 3 + 1] = 0;
		m3dSpeakerPosition[1 * 3 + 2] = 1;
		// I suppose technically the second pair should be straight left & right,
		// but I prefer moving them a bit back to mirror the front speakers.
		m3dSpeakerPosition[2 * 3 + 0] = 2;
		m3dSpeakerPosition[2 * 3 + 1] = 0;
		m3dSpeakerPosition[2 * 3 + 2] = -1;
		m3dSpeakerPosition[3 * 3 + 0] = -2;
		m3dSpeakerPosition[3 * 3 + 1] = 0;
		m3dSpeakerPosition[3 * 3 + 2] = -1;
		break;
	case 6:
		m3dSpeakerPosition[0 * 3 + 0] = 2;
		m3dSpeakerPosition[0 * 3 + 1] = 0;
		m3dSpeakerPosition[0 * 3 + 2] = 1;
		m3dSpeakerPosition[1 * 3 + 0] = -2;
		m3dSpeakerPosition[1 * 3 + 1] = 0;
		m3dSpeakerPosition[1 * 3 + 2] = 1;

		// center and subwoofer.
		m3dSpeakerPosition[2 * 3 + 0] = 0;
		m3dSpeakerPosition[2 * 3 + 1] = 0;
		m3dSpeakerPosition[2 * 3 + 2] = 1;
		// Sub should be "mix of everything". We'll handle it as a special case and make it a null vector.
		m3dSpeakerPosition[3 * 3 + 0] = 0;
		m3dSpeakerPosition[3 * 3 + 1] = 0;
		m3dSpeakerPosition[3 * 3 + 2] = 0;

		// I suppose technically the second pair should be straight left & right,
		// but I prefer moving them a bit back to mirror the front speakers.
		m3dSpeakerPosition[4 * 3 + 0] = 2;
		m3dSpeakerPosition[4 * 3 + 1] = 0;
		m3dSpeakerPosition[4 * 3 + 2] = -1;
		m3dSpeakerPosition[5 * 3 + 0] = -2;
		m3dSpeakerPosition[5 * 3 + 1] = 0;
		m3dSpeakerPosition[5 * 3 + 2] = -1;
		break;
	case 8:
		m3dSpeakerPosition[0 * 3 + 0] = 2;
		m3dSpeakerPosition[0 * 3 + 1] = 0;
		m3dSpeakerPosition[0 * 3 + 2] = 1;
		m3dSpeakerPosition[1 * 3 + 0] = -2;
		m3dSpeakerPosition[1 * 3 + 1] = 0;
		m3dSpeakerPosition[1 * 3 + 2] = 1;

		// center and subwoofer.
		m3dSpeakerPosition[2 * 3 + 0] = 0;
		m3dSpeakerPosition[2 * 3 + 1] = 0;
		m3dSpeakerPosition[2 * 3 + 2] = 1;
		// Sub should be "mix of everything". We'll handle it as a special case and make it a null vector.
		m3dSpeakerPosition[3 * 3 + 0] = 0;
		m3dSpeakerPosition[3 * 3 + 1] = 0;
		m3dSpeakerPosition[3 * 3 + 2] = 0;

		// side
		m3dSpeakerPosition[4 * 3 + 0] = 2;
		m3dSpeakerPosition[4 * 3 + 1] = 0;
		m3dSpeakerPosition[4 * 3 + 2] = 0;
		m3dSpeakerPosition[5 * 3 + 0] = -2;
		m3dSpeakerPosition[5 * 3 + 1] = 0;
		m3dSpeakerPosition[5 * 3 + 2] = 0;

		// back
		m3dSpeakerPosition[6 * 3 + 0] = 2;
		m3dSpeakerPosition[6 * 3 + 1] = 0;
		m3dSpeakerPosition[6 * 3 + 2] = -1;
		m3dSpeakerPosition[7 * 3 + 0] = -2;
		m3dSpeakerPosition[7 * 3 + 1] = 0;
		m3dSpeakerPosition[7 * 3 + 2] = -1;
		break;
	}
}

const char *Soloud::getErrorString(result aErrorCode) const
{
	switch (aErrorCode)
	{
	case SO_NO_ERROR:
		return "No error";
	case INVALID_PARAMETER:
		return "Some parameter is invalid";
	case FILE_NOT_FOUND:
		return "File not found";
	case FILE_LOAD_FAILED:
		return "File found, but could not be loaded";
	case DLL_NOT_FOUND:
		return "DLL not found, or wrong DLL";
	case OUT_OF_MEMORY:
		return "Out of memory";
	case NOT_IMPLEMENTED:
		return "Feature not implemented";
		/*case UNKNOWN_ERROR: return "Other error";*/
	}
	return "Other error";
}

float *Soloud::getWave()
{
	int i;
	lockAudioMutex_internal();
	for (i = 0; i < 256; i++)
		mWaveData[i] = mVisualizationWaveData[i];
	unlockAudioMutex_internal();
	return mWaveData;
}

float Soloud::getApproximateVolume(unsigned int aChannel)
{
	if (aChannel > mChannels)
		return 0;
	float vol = 0;
	lockAudioMutex_internal();
	vol = mVisualizationChannelVolume[aChannel];
	unlockAudioMutex_internal();
	return vol;
}

float *Soloud::calcFFT()
{
	lockAudioMutex_internal();
	float temp[1024];
	int i;
	for (i = 0; i < 256; i++)
	{
		temp[i * 2] = mVisualizationWaveData[i];
		temp[i * 2 + 1] = 0;
		temp[i + 512] = 0;
		temp[i + 768] = 0;
	}
	unlockAudioMutex_internal();

	SoLoud::FFT::fft1024(temp);

	for (i = 0; i < 256; i++)
	{
		float real = temp[i * 2];
		float imag = temp[i * 2 + 1];
		mFFTData[i] = (float)sqrt(real * real + imag * imag);
	}

	return mFFTData;
}

#if defined(SOLOUD_SSE_INTRINSICS)
void Soloud::clip_internal(AlignedFloatBuffer &aBuffer, AlignedFloatBuffer &aDestBuffer, unsigned int aSamples, float aVolume0, float aVolume1)
{
	float vd = (aVolume1 - aVolume0) / aSamples;
	float v = aVolume0;
	unsigned int i, j, c, d;
	unsigned int samplequads = (aSamples + 3) / 4; // rounded up

	// Clip
	if (mFlags & CLIP_ROUNDOFF)
	{
		float nb = -1.65f;
		__m128 negbound = _mm_load_ps1(&nb);
		float pb = 1.65f;
		__m128 posbound = _mm_load_ps1(&pb);
		float ls = 0.87f;
		__m128 linearscale = _mm_load_ps1(&ls);
		float cs = -0.1f;
		__m128 cubicscale = _mm_load_ps1(&cs);
		float nw = -0.9862875f;
		__m128 negwall = _mm_load_ps1(&nw);
		float pw = 0.9862875f;
		__m128 poswall = _mm_load_ps1(&pw);
		__m128 postscale = _mm_load_ps1(&mPostClipScaler);
		TinyAlignedFloatBuffer volumes;
		volumes.mData[0] = v;
		volumes.mData[1] = v + vd;
		volumes.mData[2] = v + vd + vd;
		volumes.mData[3] = v + vd + vd + vd;
		vd *= 4;
		__m128 vdelta = _mm_load_ps1(&vd);
		c = 0;
		d = 0;
		for (j = 0; j < mChannels; j++)
		{
			__m128 vol = _mm_load_ps(volumes.mData);

			for (i = 0; i < samplequads; i++)
			{
				// float f1 = origdata[c] * v;	c++; v += vd;
				__m128 f = _mm_load_ps(&aBuffer.mData[c]);
				c += 4;
				f = _mm_mul_ps(f, vol);
				vol = _mm_add_ps(vol, vdelta);

				// float u1 = (f1 > -1.65f);
				__m128 u = _mm_cmpgt_ps(f, negbound);

				// float o1 = (f1 < 1.65f);
				__m128 o = _mm_cmplt_ps(f, posbound);

				// f1 = (0.87f * f1 - 0.1f * f1 * f1 * f1) * u1 * o1;
				__m128 lin = _mm_mul_ps(f, linearscale);
				__m128 cubic = _mm_mul_ps(f, f);
				cubic = _mm_mul_ps(cubic, f);
				cubic = _mm_mul_ps(cubic, cubicscale);
				f = _mm_add_ps(cubic, lin);

				// f1 = f1 * u1 + !u1 * -0.9862875f;
				__m128 lowmask = _mm_andnot_ps(u, negwall);
				__m128 ilowmask = _mm_and_ps(u, f);
				f = _mm_add_ps(lowmask, ilowmask);

				// f1 = f1 * o1 + !o1 * 0.9862875f;
				__m128 himask = _mm_andnot_ps(o, poswall);
				__m128 ihimask = _mm_and_ps(o, f);
				f = _mm_add_ps(himask, ihimask);

				// outdata[d] = f1 * postclip; d++;
				f = _mm_mul_ps(f, postscale);
				_mm_store_ps(&aDestBuffer.mData[d], f);
				d += 4;
			}
		}
	}
	else
	{
		float nb = -1.0f;
		__m128 negbound = _mm_load_ps1(&nb);
		float pb = 1.0f;
		__m128 posbound = _mm_load_ps1(&pb);
		__m128 postscale = _mm_load_ps1(&mPostClipScaler);
		TinyAlignedFloatBuffer volumes;
		volumes.mData[0] = v;
		volumes.mData[1] = v + vd;
		volumes.mData[2] = v + vd + vd;
		volumes.mData[3] = v + vd + vd + vd;
		vd *= 4;
		__m128 vdelta = _mm_load_ps1(&vd);
		c = 0;
		d = 0;
		for (j = 0; j < mChannels; j++)
		{
			__m128 vol = _mm_load_ps(volumes.mData);
			for (i = 0; i < samplequads; i++)
			{
				// float f1 = aBuffer.mData[c] * v; c++; v += vd;
				__m128 f = _mm_load_ps(&aBuffer.mData[c]);
				c += 4;
				f = _mm_mul_ps(f, vol);
				vol = _mm_add_ps(vol, vdelta);

				// f1 = (f1 <= -1) ? -1 : (f1 >= 1) ? 1 : f1;
				f = _mm_max_ps(f, negbound);
				f = _mm_min_ps(f, posbound);

				// aDestBuffer.mData[d] = f1 * mPostClipScaler; d++;
				f = _mm_mul_ps(f, postscale);
				_mm_store_ps(&aDestBuffer.mData[d], f);
				d += 4;
			}
		}
	}
}
#else // fallback code
void Soloud::clip_internal(AlignedFloatBuffer &aBuffer, AlignedFloatBuffer &aDestBuffer, unsigned int aSamples, float aVolume0, float aVolume1)
{
	float vd = (aVolume1 - aVolume0) / aSamples;
	float v = aVolume0;
	unsigned int i, j, c, d;
	unsigned int samplequads = (aSamples + 3) / 4; // rounded up
	// Clip
	if (mFlags & CLIP_ROUNDOFF)
	{
		c = 0;
		d = 0;
		for (j = 0; j < mChannels; j++)
		{
			v = aVolume0;
			for (i = 0; i < samplequads; i++)
			{
				float f1 = aBuffer.mData[c] * v;
				c++;
				v += vd;
				float f2 = aBuffer.mData[c] * v;
				c++;
				v += vd;
				float f3 = aBuffer.mData[c] * v;
				c++;
				v += vd;
				float f4 = aBuffer.mData[c] * v;
				c++;
				v += vd;

				f1 = (f1 <= -1.65f) ? -0.9862875f : (f1 >= 1.65f) ? 0.9862875f : (0.87f * f1 - 0.1f * f1 * f1 * f1);
				f2 = (f2 <= -1.65f) ? -0.9862875f : (f2 >= 1.65f) ? 0.9862875f : (0.87f * f2 - 0.1f * f2 * f2 * f2);
				f3 = (f3 <= -1.65f) ? -0.9862875f : (f3 >= 1.65f) ? 0.9862875f : (0.87f * f3 - 0.1f * f3 * f3 * f3);
				f4 = (f4 <= -1.65f) ? -0.9862875f : (f4 >= 1.65f) ? 0.9862875f : (0.87f * f4 - 0.1f * f4 * f4 * f4);

				aDestBuffer.mData[d] = f1 * mPostClipScaler;
				d++;
				aDestBuffer.mData[d] = f2 * mPostClipScaler;
				d++;
				aDestBuffer.mData[d] = f3 * mPostClipScaler;
				d++;
				aDestBuffer.mData[d] = f4 * mPostClipScaler;
				d++;
			}
		}
	}
	else
	{
		c = 0;
		d = 0;
		for (j = 0; j < mChannels; j++)
		{
			v = aVolume0;
			for (i = 0; i < samplequads; i++)
			{
				float f1 = aBuffer.mData[c] * v;
				c++;
				v += vd;
				float f2 = aBuffer.mData[c] * v;
				c++;
				v += vd;
				float f3 = aBuffer.mData[c] * v;
				c++;
				v += vd;
				float f4 = aBuffer.mData[c] * v;
				c++;
				v += vd;

				f1 = (f1 <= -1) ? -1 : (f1 >= 1) ? 1 : f1;
				f2 = (f2 <= -1) ? -1 : (f2 >= 1) ? 1 : f2;
				f3 = (f3 <= -1) ? -1 : (f3 >= 1) ? 1 : f3;
				f4 = (f4 <= -1) ? -1 : (f4 >= 1) ? 1 : f4;

				aDestBuffer.mData[d] = f1 * mPostClipScaler;
				d++;
				aDestBuffer.mData[d] = f2 * mPostClipScaler;
				d++;
				aDestBuffer.mData[d] = f3 * mPostClipScaler;
				d++;
				aDestBuffer.mData[d] = f4 * mPostClipScaler;
				d++;
			}
		}
	}
}
#endif

void panAndExpand(AudioSourceInstance *aVoice, float *aBuffer, unsigned int aSamplesToRead, unsigned int aBufferSize, float *aScratch, unsigned int aChannels)
{
#ifdef SOLOUD_SSE_INTRINSICS
	SOLOUD_ASSERT(((size_t)aBuffer & 0xf) == 0);
	SOLOUD_ASSERT(((size_t)aScratch & 0xf) == 0);
	SOLOUD_ASSERT(((size_t)aBufferSize & 0xf) == 0);
#endif
	float pan[MAX_CHANNELS];  // current speaker volume
	float pand[MAX_CHANNELS]; // destination speaker volume
	float pani[MAX_CHANNELS]; // speaker volume increment per sample
	unsigned int j, k;
	for (k = 0; k < aChannels; k++)
	{
		pan[k] = aVoice->mCurrentChannelVolume[k];
		pand[k] = aVoice->mChannelVolume[k] * aVoice->mOverallVolume;
		pani[k] = (pand[k] - pan[k]) / aSamplesToRead; // TODO: this is a bit inconsistent.. but it's a hack to begin with
	}

	int ofs = 0;
	switch (aChannels)
	{
	case 1: // Target is mono. Sum everything. (1->1, 2->1, 4->1, 6->1, 8->1)
		for (j = 0, ofs = 0; j < aVoice->mChannels; j++, ofs += aBufferSize)
		{
			pan[0] = aVoice->mCurrentChannelVolume[0];
			for (k = 0; k < aSamplesToRead; k++)
			{
				pan[0] += pani[0];
				aBuffer[k] += aScratch[ofs + k] * pan[0];
			}
		}
		break;
	case 2:
		switch (aVoice->mChannels)
		{
		case 8: // 8->2, just sum lefties and righties, add a bit of center and sub?
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				float s1 = aScratch[j];
				float s2 = aScratch[aBufferSize + j];
				float s3 = aScratch[aBufferSize * 2 + j];
				float s4 = aScratch[aBufferSize * 3 + j];
				float s5 = aScratch[aBufferSize * 4 + j];
				float s6 = aScratch[aBufferSize * 5 + j];
				float s7 = aScratch[aBufferSize * 6 + j];
				float s8 = aScratch[aBufferSize * 7 + j];
				aBuffer[j + 0] += 0.2f * (s1 + s3 + s4 + s5 + s7) * pan[0];
				aBuffer[j + aBufferSize] += 0.2f * (s2 + s3 + s4 + s6 + s8) * pan[1];
			}
			break;
		case 6: // 6->2, just sum lefties and righties, add a bit of center and sub?
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				float s1 = aScratch[j];
				float s2 = aScratch[aBufferSize + j];
				float s3 = aScratch[aBufferSize * 2 + j];
				float s4 = aScratch[aBufferSize * 3 + j];
				float s5 = aScratch[aBufferSize * 4 + j];
				float s6 = aScratch[aBufferSize * 5 + j];
				aBuffer[j + 0] += 0.3f * (s1 + s3 + s4 + s5) * pan[0];
				aBuffer[j + aBufferSize] += 0.3f * (s2 + s3 + s4 + s6) * pan[1];
			}
			break;
		case 4: // 4->2, just sum lefties and righties
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				float s1 = aScratch[j];
				float s2 = aScratch[aBufferSize + j];
				float s3 = aScratch[aBufferSize * 2 + j];
				float s4 = aScratch[aBufferSize * 3 + j];
				aBuffer[j + 0] += 0.5f * (s1 + s3) * pan[0];
				aBuffer[j + aBufferSize] += 0.5f * (s2 + s4) * pan[1];
			}
			break;
		case 2: // 2->2
#if defined(SOLOUD_SSE_INTRINSICS)
		{
			int c = 0;
			// if ((aBufferSize & 3) == 0)
			{
				unsigned int samplequads = aSamplesToRead / 4; // rounded down
				TinyAlignedFloatBuffer pan0;
				pan0.mData[0] = pan[0] + pani[0];
				pan0.mData[1] = pan[0] + pani[0] * 2;
				pan0.mData[2] = pan[0] + pani[0] * 3;
				pan0.mData[3] = pan[0] + pani[0] * 4;
				TinyAlignedFloatBuffer pan1;
				pan1.mData[0] = pan[1] + pani[1];
				pan1.mData[1] = pan[1] + pani[1] * 2;
				pan1.mData[2] = pan[1] + pani[1] * 3;
				pan1.mData[3] = pan[1] + pani[1] * 4;
				pani[0] *= 4;
				pani[1] *= 4;
				__m128 pan0delta = _mm_load_ps1(&pani[0]);
				__m128 pan1delta = _mm_load_ps1(&pani[1]);
				__m128 p0 = _mm_load_ps(pan0.mData);
				__m128 p1 = _mm_load_ps(pan1.mData);

				for (j = 0; j < samplequads; j++)
				{
					__m128 f0 = _mm_load_ps(aScratch + c);
					__m128 c0 = _mm_mul_ps(f0, p0);
					__m128 f1 = _mm_load_ps(aScratch + c + aBufferSize);
					__m128 c1 = _mm_mul_ps(f1, p1);
					__m128 o0 = _mm_load_ps(aBuffer + c);
					__m128 o1 = _mm_load_ps(aBuffer + c + aBufferSize);
					c0 = _mm_add_ps(c0, o0);
					c1 = _mm_add_ps(c1, o1);
					_mm_store_ps(aBuffer + c, c0);
					_mm_store_ps(aBuffer + c + aBufferSize, c1);
					p0 = _mm_add_ps(p0, pan0delta);
					p1 = _mm_add_ps(p1, pan1delta);
					c += 4;
				}
			}

			// If buffer size or samples to read are not divisible by 4, handle leftovers
			for (j = c; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				float s1 = aScratch[j];
				float s2 = aScratch[aBufferSize + j];
				aBuffer[j + 0] += s1 * pan[0];
				aBuffer[j + aBufferSize] += s2 * pan[1];
			}
		}
#else // fallback
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				float s1 = aScratch[j];
				float s2 = aScratch[aBufferSize + j];
				aBuffer[j + 0] += s1 * pan[0];
				aBuffer[j + aBufferSize] += s2 * pan[1];
			}
#endif
		break;
		case 1: // 1->2
#if defined(SOLOUD_SSE_INTRINSICS)
		{
			int c = 0;
			// if ((aBufferSize & 3) == 0)
			{
				unsigned int samplequads = aSamplesToRead / 4; // rounded down
				TinyAlignedFloatBuffer pan0;
				pan0.mData[0] = pan[0] + pani[0];
				pan0.mData[1] = pan[0] + pani[0] * 2;
				pan0.mData[2] = pan[0] + pani[0] * 3;
				pan0.mData[3] = pan[0] + pani[0] * 4;
				TinyAlignedFloatBuffer pan1;
				pan1.mData[0] = pan[1] + pani[1];
				pan1.mData[1] = pan[1] + pani[1] * 2;
				pan1.mData[2] = pan[1] + pani[1] * 3;
				pan1.mData[3] = pan[1] + pani[1] * 4;
				pani[0] *= 4;
				pani[1] *= 4;
				__m128 pan0delta = _mm_load_ps1(&pani[0]);
				__m128 pan1delta = _mm_load_ps1(&pani[1]);
				__m128 p0 = _mm_load_ps(pan0.mData);
				__m128 p1 = _mm_load_ps(pan1.mData);

				for (j = 0; j < samplequads; j++)
				{
					__m128 f = _mm_load_ps(aScratch + c);
					__m128 c0 = _mm_mul_ps(f, p0);
					__m128 c1 = _mm_mul_ps(f, p1);
					__m128 o0 = _mm_load_ps(aBuffer + c);
					__m128 o1 = _mm_load_ps(aBuffer + c + aBufferSize);
					c0 = _mm_add_ps(c0, o0);
					c1 = _mm_add_ps(c1, o1);
					_mm_store_ps(aBuffer + c, c0);
					_mm_store_ps(aBuffer + c + aBufferSize, c1);
					p0 = _mm_add_ps(p0, pan0delta);
					p1 = _mm_add_ps(p1, pan1delta);
					c += 4;
				}
			}
			// If buffer size or samples to read are not divisible by 4, handle leftovers
			for (j = c; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				float s = aScratch[j];
				aBuffer[j + 0] += s * pan[0];
				aBuffer[j + aBufferSize] += s * pan[1];
			}
		}
#else // fallback
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				float s = aScratch[j];
				aBuffer[j + 0] += s * pan[0];
				aBuffer[j + aBufferSize] += s * pan[1];
			}
#endif
		break;
		}
		break;
	case 4:
		switch (aVoice->mChannels)
		{
		case 8: // 8->4, add a bit of center, sub?
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				pan[2] += pani[2];
				pan[3] += pani[3];
				float s1 = aScratch[j];
				float s2 = aScratch[aBufferSize + j];
				float s3 = aScratch[aBufferSize * 2 + j];
				float s4 = aScratch[aBufferSize * 3 + j];
				float s5 = aScratch[aBufferSize * 4 + j];
				float s6 = aScratch[aBufferSize * 5 + j];
				float s7 = aScratch[aBufferSize * 6 + j];
				float s8 = aScratch[aBufferSize * 7 + j];
				float c = (s3 + s4) * 0.7f;
				aBuffer[j + 0] += s1 * pan[0] + c;
				aBuffer[j + aBufferSize] += s2 * pan[1] + c;
				aBuffer[j + aBufferSize * 2] += 0.5f * (s5 + s7) * pan[2];
				aBuffer[j + aBufferSize * 3] += 0.5f * (s6 + s8) * pan[3];
			}
			break;
		case 6: // 6->4, add a bit of center, sub?
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				pan[2] += pani[2];
				pan[3] += pani[3];
				float s1 = aScratch[j];
				float s2 = aScratch[aBufferSize + j];
				float s3 = aScratch[aBufferSize * 2 + j];
				float s4 = aScratch[aBufferSize * 3 + j];
				float s5 = aScratch[aBufferSize * 4 + j];
				float s6 = aScratch[aBufferSize * 5 + j];
				float c = (s3 + s4) * 0.7f;
				aBuffer[j + 0] += s1 * pan[0] + c;
				aBuffer[j + aBufferSize] += s2 * pan[1] + c;
				aBuffer[j + aBufferSize * 2] += s5 * pan[2];
				aBuffer[j + aBufferSize * 3] += s6 * pan[3];
			}
			break;
		case 4: // 4->4
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				pan[2] += pani[2];
				pan[3] += pani[3];
				float s1 = aScratch[j];
				float s2 = aScratch[aBufferSize + j];
				float s3 = aScratch[aBufferSize * 2 + j];
				float s4 = aScratch[aBufferSize * 3 + j];
				aBuffer[j + 0] += s1 * pan[0];
				aBuffer[j + aBufferSize] += s2 * pan[1];
				aBuffer[j + aBufferSize * 2] += s3 * pan[2];
				aBuffer[j + aBufferSize * 3] += s4 * pan[3];
			}
			break;
		case 2: // 2->4
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				pan[2] += pani[2];
				pan[3] += pani[3];
				float s1 = aScratch[j];
				float s2 = aScratch[aBufferSize + j];
				aBuffer[j + 0] += s1 * pan[0];
				aBuffer[j + aBufferSize] += s2 * pan[1];
				aBuffer[j + aBufferSize * 2] += s1 * pan[2];
				aBuffer[j + aBufferSize * 3] += s2 * pan[3];
			}
			break;
		case 1: // 1->4
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				pan[2] += pani[2];
				pan[3] += pani[3];
				float s = aScratch[j];
				aBuffer[j + 0] += s * pan[0];
				aBuffer[j + aBufferSize] += s * pan[1];
				aBuffer[j + aBufferSize * 2] += s * pan[2];
				aBuffer[j + aBufferSize * 3] += s * pan[3];
			}
			break;
		}
		break;
	case 6:
		switch (aVoice->mChannels)
		{
		case 8: // 8->6
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				pan[2] += pani[2];
				pan[3] += pani[3];
				pan[4] += pani[4];
				pan[5] += pani[5];
				float s1 = aScratch[j];
				float s2 = aScratch[aBufferSize + j];
				float s3 = aScratch[aBufferSize * 2 + j];
				float s4 = aScratch[aBufferSize * 3 + j];
				float s5 = aScratch[aBufferSize * 4 + j];
				float s6 = aScratch[aBufferSize * 5 + j];
				float s7 = aScratch[aBufferSize * 6 + j];
				float s8 = aScratch[aBufferSize * 7 + j];
				aBuffer[j + 0] += s1 * pan[0];
				aBuffer[j + aBufferSize] += s2 * pan[1];
				aBuffer[j + aBufferSize * 2] += s3 * pan[2];
				aBuffer[j + aBufferSize * 3] += s4 * pan[3];
				aBuffer[j + aBufferSize * 4] += 0.5f * (s5 + s7) * pan[4];
				aBuffer[j + aBufferSize * 5] += 0.5f * (s6 + s8) * pan[5];
			}
			break;
		case 6: // 6->6
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				pan[2] += pani[2];
				pan[3] += pani[3];
				pan[4] += pani[4];
				pan[5] += pani[5];
				float s1 = aScratch[j];
				float s2 = aScratch[aBufferSize + j];
				float s3 = aScratch[aBufferSize * 2 + j];
				float s4 = aScratch[aBufferSize * 3 + j];
				float s5 = aScratch[aBufferSize * 4 + j];
				float s6 = aScratch[aBufferSize * 5 + j];
				aBuffer[j + 0] += s1 * pan[0];
				aBuffer[j + aBufferSize] += s2 * pan[1];
				aBuffer[j + aBufferSize * 2] += s3 * pan[2];
				aBuffer[j + aBufferSize * 3] += s4 * pan[3];
				aBuffer[j + aBufferSize * 4] += s5 * pan[4];
				aBuffer[j + aBufferSize * 5] += s6 * pan[5];
			}
			break;
		case 4: // 4->6
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				pan[2] += pani[2];
				pan[3] += pani[3];
				pan[4] += pani[4];
				pan[5] += pani[5];
				float s1 = aScratch[j];
				float s2 = aScratch[aBufferSize + j];
				float s3 = aScratch[aBufferSize * 2 + j];
				float s4 = aScratch[aBufferSize * 3 + j];
				aBuffer[j + 0] += s1 * pan[0];
				aBuffer[j + aBufferSize] += s2 * pan[1];
				aBuffer[j + aBufferSize * 2] += 0.5f * (s1 + s2) * pan[2];
				aBuffer[j + aBufferSize * 3] += 0.25f * (s1 + s2 + s3 + s4) * pan[3];
				aBuffer[j + aBufferSize * 4] += s3 * pan[4];
				aBuffer[j + aBufferSize * 5] += s4 * pan[5];
			}
			break;
		case 2: // 2->6
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				pan[2] += pani[2];
				pan[3] += pani[3];
				pan[4] += pani[4];
				pan[5] += pani[5];
				float s1 = aScratch[j];
				float s2 = aScratch[aBufferSize + j];
				aBuffer[j + 0] += s1 * pan[0];
				aBuffer[j + aBufferSize] += s2 * pan[1];
				aBuffer[j + aBufferSize * 2] += 0.5f * (s1 + s2) * pan[2];
				aBuffer[j + aBufferSize * 3] += 0.5f * (s1 + s2) * pan[3];
				aBuffer[j + aBufferSize * 4] += s1 * pan[4];
				aBuffer[j + aBufferSize * 5] += s2 * pan[5];
			}
			break;
		case 1: // 1->6
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				pan[2] += pani[2];
				pan[3] += pani[3];
				pan[4] += pani[4];
				pan[5] += pani[5];
				float s = aScratch[j];
				aBuffer[j + 0] += s * pan[0];
				aBuffer[j + aBufferSize] += s * pan[1];
				aBuffer[j + aBufferSize * 2] += s * pan[2];
				aBuffer[j + aBufferSize * 3] += s * pan[3];
				aBuffer[j + aBufferSize * 4] += s * pan[4];
				aBuffer[j + aBufferSize * 5] += s * pan[5];
			}
			break;
		}
		break;
	case 8:
		switch (aVoice->mChannels)
		{
		case 8: // 8->8
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				pan[2] += pani[2];
				pan[3] += pani[3];
				pan[4] += pani[4];
				pan[5] += pani[5];
				pan[6] += pani[6];
				pan[7] += pani[7];
				float s1 = aScratch[j];
				float s2 = aScratch[aBufferSize + j];
				float s3 = aScratch[aBufferSize * 2 + j];
				float s4 = aScratch[aBufferSize * 3 + j];
				float s5 = aScratch[aBufferSize * 4 + j];
				float s6 = aScratch[aBufferSize * 5 + j];
				float s7 = aScratch[aBufferSize * 6 + j];
				float s8 = aScratch[aBufferSize * 7 + j];
				aBuffer[j + 0] += s1 * pan[0];
				aBuffer[j + aBufferSize] += s2 * pan[1];
				aBuffer[j + aBufferSize * 2] += s3 * pan[2];
				aBuffer[j + aBufferSize * 3] += s4 * pan[3];
				aBuffer[j + aBufferSize * 4] += s5 * pan[4];
				aBuffer[j + aBufferSize * 5] += s6 * pan[5];
				aBuffer[j + aBufferSize * 6] += s7 * pan[6];
				aBuffer[j + aBufferSize * 7] += s8 * pan[7];
			}
			break;
		case 6: // 6->8
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				pan[2] += pani[2];
				pan[3] += pani[3];
				pan[4] += pani[4];
				pan[5] += pani[5];
				pan[6] += pani[6];
				pan[7] += pani[7];
				float s1 = aScratch[j];
				float s2 = aScratch[aBufferSize + j];
				float s3 = aScratch[aBufferSize * 2 + j];
				float s4 = aScratch[aBufferSize * 3 + j];
				float s5 = aScratch[aBufferSize * 4 + j];
				float s6 = aScratch[aBufferSize * 5 + j];
				aBuffer[j + 0] += s1 * pan[0];
				aBuffer[j + aBufferSize] += s2 * pan[1];
				aBuffer[j + aBufferSize * 2] += s3 * pan[2];
				aBuffer[j + aBufferSize * 3] += s4 * pan[3];
				aBuffer[j + aBufferSize * 4] += 0.5f * (s5 + s1) * pan[4];
				aBuffer[j + aBufferSize * 5] += 0.5f * (s6 + s2) * pan[5];
				aBuffer[j + aBufferSize * 6] += s5 * pan[6];
				aBuffer[j + aBufferSize * 7] += s6 * pan[7];
			}
			break;
		case 4: // 4->8
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				pan[2] += pani[2];
				pan[3] += pani[3];
				pan[4] += pani[4];
				pan[5] += pani[5];
				pan[6] += pani[6];
				pan[7] += pani[7];
				float s1 = aScratch[j];
				float s2 = aScratch[aBufferSize + j];
				float s3 = aScratch[aBufferSize * 2 + j];
				float s4 = aScratch[aBufferSize * 3 + j];
				aBuffer[j + 0] += s1 * pan[0];
				aBuffer[j + aBufferSize] += s2 * pan[1];
				aBuffer[j + aBufferSize * 2] += 0.5f * (s1 + s2) * pan[2];
				aBuffer[j + aBufferSize * 3] += 0.25f * (s1 + s2 + s3 + s4) * pan[3];
				aBuffer[j + aBufferSize * 4] += 0.5f * (s1 + s3) * pan[4];
				aBuffer[j + aBufferSize * 5] += 0.5f * (s2 + s4) * pan[5];
				aBuffer[j + aBufferSize * 6] += s3 * pan[4];
				aBuffer[j + aBufferSize * 7] += s4 * pan[5];
			}
			break;
		case 2: // 2->8
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				pan[2] += pani[2];
				pan[3] += pani[3];
				pan[4] += pani[4];
				pan[5] += pani[5];
				pan[6] += pani[6];
				pan[7] += pani[7];
				float s1 = aScratch[j];
				float s2 = aScratch[aBufferSize + j];
				aBuffer[j + 0] += s1 * pan[0];
				aBuffer[j + aBufferSize] += s2 * pan[1];
				aBuffer[j + aBufferSize * 2] += 0.5f * (s1 + s2) * pan[2];
				aBuffer[j + aBufferSize * 3] += 0.5f * (s1 + s2) * pan[3];
				aBuffer[j + aBufferSize * 4] += s1 * pan[4];
				aBuffer[j + aBufferSize * 5] += s2 * pan[5];
				aBuffer[j + aBufferSize * 6] += s1 * pan[6];
				aBuffer[j + aBufferSize * 7] += s2 * pan[7];
			}
			break;
		case 1: // 1->8
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				pan[2] += pani[2];
				pan[3] += pani[3];
				pan[4] += pani[4];
				pan[5] += pani[5];
				pan[6] += pani[6];
				pan[7] += pani[7];
				float s = aScratch[j];
				aBuffer[j + 0] += s * pan[0];
				aBuffer[j + aBufferSize] += s * pan[1];
				aBuffer[j + aBufferSize * 2] += s * pan[2];
				aBuffer[j + aBufferSize * 3] += s * pan[3];
				aBuffer[j + aBufferSize * 4] += s * pan[4];
				aBuffer[j + aBufferSize * 5] += s * pan[5];
				aBuffer[j + aBufferSize * 6] += s * pan[6];
				aBuffer[j + aBufferSize * 7] += s * pan[7];
			}
			break;
		}
		break;
	}

	for (k = 0; k < aChannels; k++)
		aVoice->mCurrentChannelVolume[k] = pand[k];
}

// Helper function to ensure we have enough source data in the resample buffer
unsigned int ensureSourceData_internal(AudioSourceInstance *voice, unsigned int samplesNeeded, float *scratchBuffer, unsigned int scratchSize)
{
	const unsigned int CHUNK_SIZE = AudioSourceInstance::CHUNK_SIZE;
	const unsigned int RESAMPLE_BUFFER_SIZE = AudioSourceInstance::RESAMPLE_BUFFER_SIZE;

	// Calculate how many samples we currently have available
	unsigned int availableSamples = voice->mResampleBufferFill - voice->mResampleBufferPos;

	// If we have enough samples, return early
	if (availableSamples >= samplesNeeded)
	{
		return availableSamples;
	}

	// Compact buffer when read position advances significantly
	if (voice->mResampleBufferPos >= CHUNK_SIZE && availableSamples > 0)
	{
		for (unsigned int ch = 0; ch < voice->mChannels; ch++)
		{
			memmove(voice->mResampleBuffer[ch], voice->mResampleBuffer[ch] + voice->mResampleBufferPos, availableSamples * sizeof(float));
		}
		voice->mResampleBufferFill = availableSamples;
		voice->mResampleBufferPos = 0;
	}
	else if (availableSamples == 0)
	{
		// Buffer is empty, reset positions
		voice->mResampleBufferFill = 0;
		voice->mResampleBufferPos = 0;
	}

	// Determine how much space we have for new data
	unsigned int spaceAvailable = RESAMPLE_BUFFER_SIZE - voice->mResampleBufferFill;
	if (spaceAvailable == 0)
	{
		return availableSamples; // Buffer is full
	}

	// Read one chunk at a time to maintain low latency
	unsigned int samplesToRead = CHUNK_SIZE;
	if (samplesToRead > spaceAvailable)
		samplesToRead = spaceAvailable;
	if (samplesToRead > scratchSize / voice->mChannels)
		samplesToRead = scratchSize / voice->mChannels;

	if (samplesToRead == 0)
	{
		return availableSamples;
	}

	unsigned int samplesRead = 0;
	bool shouldTryToRead = !voice->hasEnded() || (voice->mFlags & AudioSourceInstance::LOOPING) || availableSamples < samplesNeeded;

	if (shouldTryToRead)
	{
		// Use scratch buffer for reading
		float *channelBuffer = scratchBuffer;
		unsigned int channelBufferSize = samplesToRead * voice->mChannels;

		if (channelBufferSize <= scratchSize)
		{
			samplesRead = voice->getAudio(channelBuffer, samplesToRead, samplesToRead);

			// Handle looping
			if (samplesRead < samplesToRead && (voice->mFlags & AudioSourceInstance::LOOPING))
			{
				while (samplesRead < samplesToRead && voice->seek(voice->mLoopPoint, channelBuffer + samplesRead * voice->mChannels,
				                                                  (samplesToRead - samplesRead) * voice->mChannels) == SO_NO_ERROR)
				{
					voice->mLoopCount++;
					unsigned int loopSamples =
					    voice->getAudio(channelBuffer + samplesRead * voice->mChannels, samplesToRead - samplesRead, samplesToRead - samplesRead);
					samplesRead += loopSamples;
					if (loopSamples == 0)
						break;
				}
			}

			// Apply filters to source data
			if (samplesRead > 0)
			{
				for (unsigned int j = 0; j < FILTERS_PER_STREAM; j++)
				{
					if (voice->mFilter[j])
					{
						voice->mFilter[j]->filter(channelBuffer, samplesRead, samplesToRead, voice->mChannels, voice->mSamplerate, voice->mStreamTime);
					}
				}
			}

			// Copy from channel-separated scratch buffer to resample buffers
			for (unsigned int ch = 0; ch < voice->mChannels; ch++)
			{
				float *srcChannel = channelBuffer + ch * samplesToRead;
				float *dstChannel = voice->mResampleBuffer[ch] + voice->mResampleBufferFill;
				memcpy(dstChannel, srcChannel, samplesRead * sizeof(float));

				// Clear remaining space if we read fewer samples than requested
				if (samplesRead < samplesToRead)
				{
					memset(dstChannel + samplesRead, 0, (samplesToRead - samplesRead) * sizeof(float));
				}
			}
		}
	}

	voice->mResampleBufferFill += samplesRead;
	return voice->mResampleBufferFill - voice->mResampleBufferPos;
}

// Improved resampling function that handles arbitrary ratios precisely
unsigned int Soloud::resampleVoicePrecise_internal(AudioSourceInstance *voice,
                                                   float *outputBuffer,
                                                   unsigned int outputSamples,
                                                   unsigned int outputStride,
                                                   double outputSampleRate,
                                                   unsigned int resampler,
                                                   float *scratchBuffer,
                                                   unsigned int scratchSize)
{
	if (outputSamples == 0 || !voice)
		return 0;

	// Calculate step size (how much we advance in source per output sample)
	double stepSize = voice->mSamplerate / outputSampleRate;

	// Determine lookahead requirements
	unsigned int lookaheadSamples = 0;
	switch (resampler)
	{
	case RESAMPLER_POINT:
		lookaheadSamples = 1;
		break;
	case RESAMPLER_LINEAR:
		lookaheadSamples = 2;
		break;
	case RESAMPLER_CATMULLROM:
		lookaheadSamples = 4;
		break;
	default:
		lookaheadSamples = 2;
		break;
	}

	// For chunked processing, we may not be able to produce all requested samples
	// if we run out of source data. That's okay - we'll get more on the next call.
	unsigned int samplesProduced = 0;
	unsigned int samplesToProcess = outputSamples;

	// Limit processing to one chunk worth of output to maintain low latency
	// This prevents excessive buffering for high sample rate ratios
	const unsigned int MAX_OUTPUT_CHUNK = AudioSourceInstance::CHUNK_SIZE;
	if (samplesToProcess > MAX_OUTPUT_CHUNK)
	{
		samplesToProcess = MAX_OUTPUT_CHUNK;
	}

	// Calculate how much input we need for this chunk of output
	unsigned int inputNeeded = (unsigned int)ceil(samplesToProcess * stepSize) + lookaheadSamples + 8;
	unsigned int availableInput = ensureSourceData_internal(voice, inputNeeded, scratchBuffer, scratchSize);

	if (availableInput < lookaheadSamples)
	{
		return 0; // Not enough input data available
	}

	// If no output buffer provided (tick-only mode), just advance position
	if (!outputBuffer)
	{
		// Calculate how many samples we can actually advance through
		double maxAdvance = (availableInput - lookaheadSamples) / stepSize;
		unsigned int actualSamples = (unsigned int)maxAdvance;
		if (actualSamples > samplesToProcess)
			actualSamples = samplesToProcess;

		voice->mPreciseSrcPosition += actualSamples * stepSize;

		// Update buffer position for consumed integer samples
		unsigned int integralConsumed = (unsigned int)floor(voice->mPreciseSrcPosition);
		if (integralConsumed > 0 && integralConsumed <= availableInput)
		{
			voice->mResampleBufferPos += integralConsumed;
			voice->mPreciseSrcPosition -= integralConsumed;
		}
		return actualSamples;
	}

	// Resample each channel
	for (unsigned int ch = 0; ch < voice->mChannels; ch++)
	{
		double srcPos = voice->mPreciseSrcPosition;
		float *src = voice->mResampleBuffer[ch] + voice->mResampleBufferPos;
		float *dst = outputBuffer + ch * outputStride;

		for (unsigned int i = 0; i < samplesToProcess; i++)
		{
			unsigned int intPos = (unsigned int)floor(srcPos);

			// Check bounds, stop if we don't have enough input data
			if (intPos + lookaheadSamples > availableInput)
			{
				samplesProduced = i;
				break;
			}

			double frac = srcPos - intPos;

			// Perform resampling with boundary protection
			switch (resampler)
			{
			case RESAMPLER_POINT:
				dst[i] = src[intPos];
				break;

			case RESAMPLER_LINEAR: {
				float s0 = src[intPos];
				float s1 = (intPos + 1 < availableInput) ? src[intPos + 1] : src[intPos];
				dst[i] = s0 + (s1 - s0) * (float)frac;
				break;
			}

			default:
			case RESAMPLER_CATMULLROM: {
				float p0 = (intPos >= 1) ? src[intPos - 1] : src[intPos];
				float p1 = src[intPos];
				float p2 = (intPos + 1 < availableInput) ? src[intPos + 1] : src[intPos];
				float p3 = (intPos + 2 < availableInput) ? src[intPos + 2] : src[intPos];

				float t = (float)frac;
				dst[i] = 0.5f * ((2 * p1) + (-p0 + p2) * t + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t * t + (-p0 + 3 * p1 - 3 * p2 + p3) * t * t * t);
				break;
			}
			}

			srcPos += stepSize;
		}

		// All channels should produce the same number of samples
		if (ch == 0)
		{
			samplesProduced = (samplesProduced == 0) ? samplesToProcess : samplesProduced;
		}
	}

	// Update position tracking with precise accumulation
	double totalAdvance = samplesProduced * stepSize;
	voice->mPreciseSrcPosition += totalAdvance;

	// Update buffer position for consumed integer samples
	unsigned int integralConsumed = (unsigned int)floor(voice->mPreciseSrcPosition);
	if (integralConsumed > 0 && integralConsumed <= availableInput)
	{
		voice->mResampleBufferPos += integralConsumed;
		voice->mPreciseSrcPosition -= integralConsumed; // Keep fractional part
	}

	return samplesProduced;
}

// dynamic resampling + chunked processing
void Soloud::mixBus_internal(float *aBuffer, unsigned int aSamplesToRead, unsigned int aBufferSize, float *aScratch, unsigned int aBus, float aSamplerate,
                             unsigned int aChannels, unsigned int aResampler)
{
	unsigned int i, j;

	// Clear accumulation buffer
	for (i = 0; i < aSamplesToRead; i++)
	{
		for (j = 0; j < aChannels; j++)
		{
			aBuffer[i + j * aBufferSize] = 0;
		}
	}

	// scratch space for chunked processing
	unsigned int voiceScratchSize = MAX_CHANNELS * aBufferSize;
	unsigned int tempScratchSize = 2048; // temp buffer for chunk-based reading
	unsigned int delayScratchSize = aChannels * aBufferSize;
	unsigned int totalPerVoice = voiceScratchSize + tempScratchSize + delayScratchSize;
	unsigned int maxVoicesInParallel = (mScratchSize * MAX_CHANNELS) / totalPerVoice;
	if (maxVoicesInParallel == 0)
		maxVoicesInParallel = 1;

	// Process each active voice using chunks
	for (i = 0; i < mActiveVoiceCount; i++)
	{
		AudioSourceInstance *voice = mVoice[mActiveVoice[i]];
		if (!voice || voice->mBusHandle != aBus)
			continue;

		bool isPaused = (voice->mFlags & AudioSourceInstance::PAUSED) != 0;
		bool isInaudible = (voice->mFlags & AudioSourceInstance::INAUDIBLE) != 0;
		bool mustTick = (voice->mFlags & AudioSourceInstance::INAUDIBLE_TICK) != 0;

		bool isAudible = !isPaused && !isInaudible;
		bool shouldProcess = isAudible || (!isPaused && mustTick);

		if (!shouldProcess)
			continue;

		unsigned int outputSamples = aSamplesToRead;
		unsigned int outputOffset = 0;

		// Handle delay samples
		if (voice->mDelaySamples > 0)
		{
			if (voice->mDelaySamples >= aSamplesToRead)
			{
				voice->mDelaySamples -= aSamplesToRead;
				continue;
			}
			else
			{
				outputOffset = voice->mDelaySamples;
				outputSamples = aSamplesToRead - voice->mDelaySamples;
				voice->mDelaySamples = 0;
			}
		}

		if (outputSamples == 0)
			continue;

		// Allocate scratch space for this voice
		unsigned int voiceSlot = i % maxVoicesInParallel;
		float *voiceScratch = aScratch + (voiceSlot * totalPerVoice);
		float *tempScratch = voiceScratch + voiceScratchSize;
		float *delayScratch = tempScratch + tempScratchSize;

		// Adjust temp scratch size if needed
		unsigned int actualTempScratchSize = tempScratchSize;
		if (voiceSlot == maxVoicesInParallel - 1)
		{
			unsigned int remainingSpace = (mScratchSize * MAX_CHANNELS) - (voiceSlot * totalPerVoice + voiceScratchSize + delayScratchSize);
			if (remainingSpace < tempScratchSize)
			{
				actualTempScratchSize = remainingSpace;
			}
		}

		// Process in chunks - may require multiple iterations to fill output buffer
		unsigned int totalProduced = 0;
		while (totalProduced < outputSamples)
		{
			unsigned int remainingSamples = outputSamples - totalProduced;
			unsigned int chunkOutput = remainingSamples;

			// Limit chunk size to maintain low latency
			const unsigned int MAX_CHUNK = AudioSourceInstance::CHUNK_SIZE;
			if (chunkOutput > MAX_CHUNK)
				chunkOutput = MAX_CHUNK;

			if (isAudible)
			{
				// Point to the right offset in the voice scratch buffer
				float *chunkScratch = voiceScratch;
				for (unsigned int ch = 0; ch < voice->mChannels; ch++)
				{
					// Offset each channel's scratch by totalProduced
					memset(chunkScratch + ch * aBufferSize + totalProduced, 0, chunkOutput * sizeof(float));
				}

				// Resample this chunk
				unsigned int samplesProduced = resampleVoicePrecise_internal(voice,
				                                                             chunkScratch + totalProduced, // Offset scratch buffer
				                                                             chunkOutput,
				                                                             aBufferSize,
				                                                             aSamplerate,
				                                                             aResampler,
				                                                             tempScratch,
				                                                             actualTempScratchSize);

				if (samplesProduced == 0)
				{
					// No more samples available, break out
					break;
				}

				totalProduced += samplesProduced;

				// If we got fewer samples than requested, we're done
				if (samplesProduced < chunkOutput)
				{
					break;
				}
			}
			else if (mustTick && !isPaused)
			{
				// Just advance the voice without producing audio
				unsigned int samplesAdvanced =
				    resampleVoicePrecise_internal(voice, NULL, chunkOutput, 0, aSamplerate, aResampler, tempScratch, actualTempScratchSize);

				if (samplesAdvanced == 0)
					break;

				totalProduced += samplesAdvanced;
			}
			else
			{
				break;
			}
		}

		// Mix the accumulated samples into the output buffer
		if (isAudible && totalProduced > 0)
		{
			if (outputOffset == 0)
			{
				// No delay, mix directly
				panAndExpand(voice, aBuffer, totalProduced, aBufferSize, voiceScratch, aChannels);
			}
			else
			{
				// Apply delay using delay scratch buffer
				memset(delayScratch, 0, aChannels * aBufferSize * sizeof(float));
				panAndExpand(voice, delayScratch, totalProduced, aBufferSize, voiceScratch, aChannels);

				// Mix with offset
				for (unsigned int ch = 0; ch < aChannels; ch++)
				{
					for (unsigned int s = 0; s < totalProduced; s++)
					{
						if (s + outputOffset < aSamplesToRead)
						{
							aBuffer[(s + outputOffset) + ch * aBufferSize] += delayScratch[s + ch * aBufferSize];
						}
					}
				}
			}
		}

		// Check if voice should be stopped - use small lookahead for chunked approach
		unsigned int remainingSamples = voice->mResampleBufferFill - voice->mResampleBufferPos;
		unsigned int lookaheadSamples = (aResampler == RESAMPLER_CATMULLROM) ? 4 : 2;
		bool bufferNearEmpty = remainingSamples <= lookaheadSamples;

		if (!(voice->mFlags & (AudioSourceInstance::LOOPING | AudioSourceInstance::DISABLE_AUTOSTOP)) && voice->hasEnded() && bufferNearEmpty)
		{
			stopVoice_internal(mActiveVoice[i]);
		}
	}
}

void Soloud::calcActiveVoices_internal()
{
	// TODO: consider whether we need to re-evaluate the active voices all the time.
	// It is a must when new voices are started, but otherwise we could get away
	// with postponing it sometimes..

	mActiveVoiceDirty = false;

	// Populate
	unsigned int i, candidates, mustlive;
	candidates = 0;
	mustlive = 0;
	for (i = 0; i < mHighestVoice; i++)
	{
		if (mVoice[i] &&
		    (!(mVoice[i]->mFlags & (AudioSourceInstance::INAUDIBLE | AudioSourceInstance::PAUSED)) || (mVoice[i]->mFlags & AudioSourceInstance::INAUDIBLE_TICK)))
		{
			mActiveVoice[candidates] = i;
			candidates++;
			if (mVoice[i]->mFlags & AudioSourceInstance::INAUDIBLE_TICK)
			{
				mActiveVoice[candidates - 1] = mActiveVoice[mustlive];
				mActiveVoice[mustlive] = i;
				mustlive++;
			}
		}
	}

	// Check for early out
	if (candidates <= mMaxActiveVoices)
	{
		// everything is audible, early out
		mActiveVoiceCount = candidates;
		return;
	}

	mActiveVoiceCount = mMaxActiveVoices;

	if (mustlive >= mMaxActiveVoices)
	{
		// Oopsie. Well, nothing to sort, since the "must live" voices already
		// ate all our active voice slots.
		// This is a potentially an error situation, but we have no way to report
		// error from here. And asserting could be bad, too.
		return;
	}

	// If we get this far, there's nothing to it: we'll have to sort the voices to find the most audible.

	// Iterative partial quicksort:
	int left = 0, stack[24], pos = 0, right;
	int len = candidates - mustlive;
	unsigned int *data = mActiveVoice + mustlive;
	int k = mActiveVoiceCount;
	for (;;)
	{
		for (; left + 1 < len; len++)
		{
			if (pos == 24)
				len = stack[pos = 0];
			int pivot = data[left];
			float pivotvol = mVoice[pivot]->mOverallVolume;
			stack[pos++] = len;
			for (right = left - 1;;)
			{
				do
				{
					right++;
				} while (mVoice[data[right]]->mOverallVolume > pivotvol);
				do
				{
					len--;
				} while (pivotvol > mVoice[data[len]]->mOverallVolume);
				if (right >= len)
					break;
				int temp = data[right];
				data[right] = data[len];
				data[len] = temp;
			}
		}
		if (pos == 0)
			break;
		if (left >= k)
			break;
		left = len;
		len = stack[--pos];
	}
	// TODO: should the rest of the voices be flagged INAUDIBLE?
}

void Soloud::mix_internal(unsigned int aSamples, unsigned int aStride)
{
#ifdef FLOATING_POINT_DEBUG
	// This needs to be done in the audio thread as well..
	static int done = 0;
	if (!done)
	{
		unsigned int u;
		u = _controlfp(0, 0);
		u = u & ~(_EM_INVALID | /*_EM_DENORMAL |*/ _EM_ZERODIVIDE | _EM_OVERFLOW /*| _EM_UNDERFLOW  | _EM_INEXACT*/);
		_controlfp(u, _MCW_EM);
		done = 1;
	}
#endif

#ifdef __arm__
	// flush to zero (FTZ) for ARM
	{
		static bool once = false;
		if (!once)
		{
			once = true;
			asm("vmsr fpscr,%0" ::"r"(1 << 24));
		}
	}
#endif

#ifdef _MCW_DN
	{
		static bool once = false;
		if (!once)
		{
			once = true;
			if (!(mFlags & NO_FPU_REGISTER_CHANGE))
			{
				_controlfp(_DN_FLUSH, _MCW_DN);
			}
		}
	}
#endif

#ifdef SOLOUD_SSE_INTRINSICS
	{
		static bool once = false;
		if (!once)
		{
			once = true;
			// Set denorm clear to zero (CTZ) and denorms are zero (DAZ) flags on.
			// This causes all math to consider really tiny values as zero, which
			// helps performance. I'd rather use constants from the sse headers,
			// but for some reason the DAZ value is not defined there(!)
			if (!(mFlags & NO_FPU_REGISTER_CHANGE))
			{
				_mm_setcsr(_mm_getcsr() | 0x8040);
			}
		}
	}
#endif

	time buffertime = (time)aSamples / mSamplerate;
	float globalVolume[2];
	mStreamTime += buffertime;
	mLastClockedTime = 0;

	lockAudioMutex_internal();

	// Read and update global volume inside mutex
	globalVolume[0] = mGlobalVolume;
	if (mGlobalVolumeFader.mActive)
	{
		mGlobalVolume = mGlobalVolumeFader.get(mStreamTime);
	}
	globalVolume[1] = mGlobalVolume;

	// Process faders. May change scratch size.
	int i;
	for (i = 0; i < (signed)mHighestVoice; i++)
	{
		if (mVoice[i] && !(mVoice[i]->mFlags & AudioSourceInstance::PAUSED))
		{
			float volume[2];

			mVoice[i]->mActiveFader = 0;

			if (mGlobalVolumeFader.mActive > 0)
			{
				mVoice[i]->mActiveFader = 1;
			}

			mVoice[i]->mStreamTime += buffertime;
			mVoice[i]->mStreamPosition += (double)buffertime * (double)mVoice[i]->mOverallRelativePlaySpeed;

			// TODO: this is actually unstable, because mStreamTime depends on the relative
			// play speed.
			if (mVoice[i]->mRelativePlaySpeedFader.mActive > 0)
			{
				float speed = mVoice[i]->mRelativePlaySpeedFader.get(mVoice[i]->mStreamTime);
				setVoiceRelativePlaySpeed_internal(i, speed);
			}

			volume[0] = mVoice[i]->mOverallVolume;
			if (mVoice[i]->mVolumeFader.mActive > 0)
			{
				mVoice[i]->mSetVolume = mVoice[i]->mVolumeFader.get(mVoice[i]->mStreamTime);
				mVoice[i]->mActiveFader = 1;
				updateVoiceVolume_internal(i);
				mActiveVoiceDirty = true;
			}
			volume[1] = mVoice[i]->mOverallVolume;

			if (mVoice[i]->mPanFader.mActive > 0)
			{
				float pan = mVoice[i]->mPanFader.get(mVoice[i]->mStreamTime);
				setVoicePan_internal(i, pan);
				mVoice[i]->mActiveFader = 1;
			}

			if (mVoice[i]->mPauseScheduler.mActive)
			{
				mVoice[i]->mPauseScheduler.get(mVoice[i]->mStreamTime);
				if (mVoice[i]->mPauseScheduler.mActive == -1)
				{
					mVoice[i]->mPauseScheduler.mActive = 0;
					setVoicePause_internal(i, 1);
				}
			}

			if (mVoice[i]->mStopScheduler.mActive)
			{
				mVoice[i]->mStopScheduler.get(mVoice[i]->mStreamTime);
				if (mVoice[i]->mStopScheduler.mActive == -1)
				{
					mVoice[i]->mStopScheduler.mActive = 0;
					stopVoice_internal(i);
				}
			}
		}
	}

	if (mActiveVoiceDirty)
		calcActiveVoices_internal();

	mixBus_internal(mOutputScratch.mData, aSamples, aStride, mScratch.mData, 0, (float)mSamplerate, mChannels, mResampler);

	for (i = 0; i < FILTERS_PER_STREAM; i++)
	{
		if (mFilterInstance[i])
		{
			mFilterInstance[i]->filter(mOutputScratch.mData, aSamples, aStride, mChannels, (float)mSamplerate, mStreamTime);
		}
	}

	unlockAudioMutex_internal();

	// Note: clipping channels*aStride, not channels*aSamples, so we're possibly clipping some unused data.
	// The buffers should be large enough for it, we just may do a few bytes of unneccessary work.
	clip_internal(mOutputScratch, mScratch, aStride, globalVolume[0], globalVolume[1]);

	if (mFlags & ENABLE_VISUALIZATION)
	{
		for (i = 0; i < MAX_CHANNELS; i++)
		{
			mVisualizationChannelVolume[i] = 0;
		}
		if (aSamples > 255)
		{
			for (i = 0; i < 256; i++)
			{
				int j;
				mVisualizationWaveData[i] = 0;
				for (j = 0; j < (signed)mChannels; j++)
				{
					float sample = mScratch.mData[i + j * aStride];
					float absvol = (float)fabs(sample);
					if (mVisualizationChannelVolume[j] < absvol)
						mVisualizationChannelVolume[j] = absvol;
					mVisualizationWaveData[i] += sample;
				}
			}
		}
		else
		{
			// Very unlikely failsafe branch
			for (i = 0; i < 256; i++)
			{
				int j;
				mVisualizationWaveData[i] = 0;
				for (j = 0; j < (signed)mChannels; j++)
				{
					float sample = mScratch.mData[(i % aSamples) + j * aStride];
					float absvol = (float)fabs(sample);
					if (mVisualizationChannelVolume[j] < absvol)
						mVisualizationChannelVolume[j] = absvol;
					mVisualizationWaveData[i] += sample;
				}
			}
		}
	}
}

void Soloud::mix(float *aBuffer, unsigned int aSamples)
{
	unsigned int stride = (aSamples + 15) & ~0xf;
	mix_internal(aSamples, stride);
	// 111222 -> 121212
	unsigned int i = 0, j = 0, c = 0;
	for (j = 0; j < mChannels; j++)
	{
		c = j * stride;
		for (i = j; i < aSamples * mChannels; i += mChannels)
		{
			aBuffer[i] = mScratch.mData[c];
			c++;
		}
	}
}

void Soloud::mixUnsigned8(unsigned char *aBuffer, unsigned int aSamples)
{
	unsigned int stride = (aSamples + 15) & ~0xf;
	mix_internal(aSamples, stride);
	// 111222 -> 121212, convert from [-1,1] float to [0,255] unsigned 8-bit
	unsigned int i = 0, j = 0, c = 0;
	for (j = 0; j < mChannels; j++)
	{
		c = j * stride;
		for (i = j; i < aSamples * mChannels; i += mChannels)
		{
			int sample = (int)(mScratch.mData[c] * 127.0f + 128.0f);
			if (sample < 0)
				sample = 0;
			if (sample > 255)
				sample = 255;
			aBuffer[i] = (unsigned char)sample;
			c++;
		}
	}
}

void Soloud::mixSigned16(short *aBuffer, unsigned int aSamples)
{
	unsigned int stride = (aSamples + 15) & ~0xf;
	mix_internal(aSamples, stride);
	// 111222 -> 121212
	unsigned int i = 0, j = 0, c = 0;
	for (j = 0; j < mChannels; j++)
	{
		c = j * stride;
		for (i = j; i < aSamples * mChannels; i += mChannels)
		{
			aBuffer[i] = (short)(mScratch.mData[c] * 0x7fff);
			c++;
		}
	}
}

void Soloud::mixSigned24(unsigned char *aBuffer, unsigned int aSamples)
{
	unsigned int stride = (aSamples + 15) & ~0xf;
	mix_internal(aSamples, stride);
	// 111222 -> 121212, convert from [-1,1] float to 24-bit signed (3 bytes per sample, little endian)
	unsigned int i = 0, j = 0, c = 0;
	for (j = 0; j < mChannels; j++)
	{
		c = j * stride;
		for (i = j; i < aSamples * mChannels; i += mChannels)
		{
			int sample = (int)(mScratch.mData[c] * 8388607.0f); // 0x7fffff
			if (sample < -8388608)
				sample = -8388608;
			if (sample > 8388607)
				sample = 8388607;

			// store as little endian 3-byte signed integer
			unsigned int destIdx = i * 3;
			aBuffer[destIdx] = (unsigned char)(sample & 0xff);
			aBuffer[destIdx + 1] = (unsigned char)((sample >> 8) & 0xff);
			aBuffer[destIdx + 2] = (unsigned char)((sample >> 16) & 0xff);
			c++;
		}
	}
}

void Soloud::mixSigned32(int *aBuffer, unsigned int aSamples)
{
	unsigned int stride = (aSamples + 15) & ~0xf;
	mix_internal(aSamples, stride);
	// 111222 -> 121212, convert from [-1,1] float to 32-bit signed
	unsigned int i = 0, j = 0, c = 0;
	for (j = 0; j < mChannels; j++)
	{
		c = j * stride;
		for (i = j; i < aSamples * mChannels; i += mChannels)
		{
			// use double precision for better accuracy with 32-bit range
			double sample = (double)mScratch.mData[c] * 2147483647.0;
			if (sample < -2147483648.0)
				sample = -2147483648.0;
			if (sample > 2147483647.0)
				sample = 2147483647.0;
			aBuffer[i] = (int)sample;
			c++;
		}
	}
}

void Soloud::lockAudioMutex_internal()
{
	if (mAudioThreadMutex)
	{
		Thread::lockMutex(mAudioThreadMutex);
	}
	SOLOUD_ASSERT(!mInsideAudioThreadMutex);
	mInsideAudioThreadMutex = true;
}

void Soloud::unlockAudioMutex_internal()
{
	SOLOUD_ASSERT(mInsideAudioThreadMutex);
	mInsideAudioThreadMutex = false;
	if (mAudioThreadMutex)
	{
		Thread::unlockMutex(mAudioThreadMutex);
	}
}

}; // namespace SoLoud
