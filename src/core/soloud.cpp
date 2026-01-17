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

#include "soloud_audiosource.h"
#include "soloud_cpu.h"
#include "soloud_fft.h"
#include "soloud_internal.h"
#include "soloud_mixing_internal.h"
#include "soloud_thread.h"

#include <cmath> // sin
#include <cstdio>
#include <cstdlib>
#include <cstring>

// #define FLOATING_POINT_DEBUG

#ifdef FLOATING_POINT_DEBUG
#include <cfloat> // _controlfp
#endif

namespace SoLoud
{

using namespace mixing;

std::unique_ptr<Mixer> Soloud::mMixer{nullptr};
Soloud::Soloud()
{
	// Initialize alignment masks and such.
	initCPUFeatures();
#ifdef FLOATING_POINT_DEBUG
	unsigned int u;
	u = _controlfp(0, 0);
	u = u & ~(_EM_INVALID | /*_EM_DENORMAL |*/ _EM_ZERODIVIDE | _EM_OVERFLOW /*| _EM_UNDERFLOW  | _EM_INEXACT*/);
	_controlfp(u, _MCW_EM);
#endif
	// Create the runtime-dispatched mixer.
	if (!mMixer)
	{
		mMixer.reset(Mixer::createMixer());
	}

	mResampler = SOLOUD_DEFAULT_RESAMPLER;
	mInsideAudioThreadMutex = false;
	mScratchSize = 0;
	mSamplerate = 0;
	mBufferSize = 0;
	mFlags = 0;
	mGlobalVolume = 0;
	mPlayIndex = 0;
	mBackendData = nullptr;
	mAudioThreadMutex = nullptr;
	mPostClipScaler = 0;
	mBackendCleanupFunc = nullptr;
	mBackendPauseFunc = nullptr;
	mBackendResumeFunc = nullptr;
	mEnumerateDevicesFunc = nullptr;
	mGetCurrentDeviceFunc = nullptr;
	mSetDeviceFunc = nullptr;
	mDeviceList = nullptr;
	mDeviceCount = 0;
	mChannels = 2;
	mStreamTime = 0;
	mLastClockedTime = 0;
	mAudioSourceID = 1;
	mBackendString = nullptr;
	mBackendID = 0;
	mActiveVoiceDirty = true;
	mActiveVoiceCount = 0;
	mActiveVoice = {};
	mFilter = {};
	mFilterInstance = {};
	mFFTData = {};
	mVisualizationWaveData = {};
	mWaveData = {};
	mVisualizationChannelVolume = {};
	mVoiceGroup = nullptr;
	mVoiceGroupCount = 0;
	mVoice = {};

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
	m3dSpeakerPosition = {};
}

namespace FFmpeg::FFmpegLoader
{
// don't need the entire header here
extern void cleanup();
} // namespace FFmpeg::FFmpegLoader

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

	// As the final step, unload any FFmpeg libraries we might have dynamically loaded (through Wavs/WavStreams).
	// Maybe slightly out of place, but there's no point in overcomplicating matters. We know for sure that it's safe
	// to do at this point.
	FFmpeg::FFmpegLoader::cleanup();
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
	mBackendCleanupFunc = nullptr;
	if (mAudioThreadMutex)
		Thread::destroyMutex(mAudioThreadMutex);
	mAudioThreadMutex = nullptr;
	mEnumerateDevicesFunc = nullptr;
	mGetCurrentDeviceFunc = nullptr;
	mSetDeviceFunc = nullptr;
	delete[] mDeviceList;
	mDeviceList = nullptr;
	mDeviceCount = 0;
}

result Soloud::init(unsigned int aFlags, unsigned int aBackend, unsigned int aSamplerate, unsigned int aBufferSize, unsigned int aChannels)
{
	if (aBackend >= BACKEND_MAX || aChannels == 3 || aChannels == 5 || aChannels == 7 || aChannels > MAX_CHANNELS)
		return INVALID_PARAMETER;

	deinit();

	mAudioThreadMutex = Thread::createMutex();

	mBackendID = 0;
	mBackendString = nullptr;

	int samplerate = 44100;
	int buffersize = 2048;
	int inited = 0;
	int failed = 0;

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
		failed = miniaudio_init(this, aFlags, samplerate, buffersize, aChannels);
		if (failed == 0)
		{
			inited = 1;
			mBackendID = Soloud::MINIAUDIO;
		}
	}

#if defined(WITH_SDL3)
	if (!inited && (aBackend == Soloud::SDL3 || aBackend == Soloud::AUTO || failed))
	{
		failed = sdl3_init(this, aFlags, samplerate, buffersize, aChannels);
		if (failed == 0)
		{
			inited = 1;
			mBackendID = Soloud::SDL3;
		}
	}
#endif

	if (failed != 0 && aBackend != Soloud::AUTO)
		return failed;

	if (!inited && (aBackend == Soloud::NOSOUND || aBackend == Soloud::AUTO))
	{
		if (aBufferSize == Soloud::AUTO)
			buffersize = 2048;
		if (aSamplerate == Soloud::AUTO)
			samplerate = 44100;

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
		if (aSamplerate == Soloud::AUTO)
			samplerate = 44100;

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
	mScratchSize = (aBufferSize + CPU_ALIGNMENT_MASK()) & ~CPU_ALIGNMENT_MASK();
	if (mScratchSize < SAMPLE_GRANULARITY * 4) // 4096
		mScratchSize = SAMPLE_GRANULARITY * 4;

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
	return mWaveData.data();
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
		mFFTData[i] = (float)std::sqrt(real * real + imag * imag);
	}

	return mFFTData.data();
}

// Helper function to ensure we have enough source data in the resample buffer
unsigned int Soloud::ensureSourceData_internal(AudioSourceInstance *voice, unsigned int samplesNeeded, float *scratchBuffer, unsigned int scratchSize)
{
	// Calculate how many samples we currently have available
	unsigned int availableSamples = voice->mResampleBufferFill - voice->mResampleBufferPos;

	// If we have enough samples, return early
	if (availableSamples >= samplesNeeded)
	{
		return availableSamples;
	}

	// Compact buffer when read position advances significantly
	// This prevents the buffer from growing indefinitely due to fractional positioning
	if (voice->mResampleBufferPos >= SAMPLE_GRANULARITY && availableSamples > 0)
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
	unsigned int spaceAvailable = AudioSourceInstance::RESAMPLE_BUFFER_SIZE - voice->mResampleBufferFill;
	if (spaceAvailable == 0)
	{
		return availableSamples; // Buffer is full
	}

	// Read only as much data is necessary (rounded to a power of 2, and with a reasonable minimum), up to SAMPLE_GRANULARITY, to maintain low latency
	unsigned int samplesToRead = samplesNeeded < 64U ? 64U : samplesNeeded;
	if (samplesToRead & (samplesToRead - 1))
	{
		unsigned int i = 1;
		while (i < samplesToRead)
			i *= 2;
		samplesToRead = i;
	}
	if (samplesToRead > SAMPLE_GRANULARITY)
		samplesToRead = SAMPLE_GRANULARITY;
	if (samplesToRead > spaceAvailable)
		samplesToRead = spaceAvailable;
	if (samplesToRead > scratchSize / voice->mChannels)
		samplesToRead = scratchSize / voice->mChannels;

	if (samplesToRead == 0)
		return availableSamples;

	unsigned int samplesRead = 0;
	bool shouldTryToRead = !voice->hasEnded() || (voice->mFlags & AudioSourceInstance::LOOPING) || availableSamples < samplesNeeded;

	if (shouldTryToRead)
	{
		// Use scratch buffer for reading - this buffer is channel-interleaved
		float *channelBuffer = scratchBuffer;
		unsigned int channelBufferSize = samplesToRead * voice->mChannels;

		if (channelBufferSize <= scratchSize)
		{
			// Keep aligned for optimal intrinsics perf
			unsigned int alignedBufferSize = (samplesToRead + CPU_ALIGNMENT_MASK()) & ~CPU_ALIGNMENT_MASK();
			// Make sure we don't exceed scratch space
			if (alignedBufferSize * voice->mChannels > scratchSize)
			{
				alignedBufferSize = (scratchSize / voice->mChannels) & ~CPU_ALIGNMENT_MASK();
				if (alignedBufferSize < samplesToRead)
					samplesToRead = alignedBufferSize;
			}

			samplesRead = voice->getAudio(channelBuffer, samplesToRead, alignedBufferSize);

			// Handle looping: continue reading from loop point if we reach the end
			if (samplesRead < samplesToRead && (voice->mFlags & AudioSourceInstance::LOOPING))
			{
				// crossfade ~1.5ms at 44.1kHz, long enough to eliminate clicks, but short enough to not be noticeable
				const unsigned int LOOP_CROSSFADE_SAMPLES = 64;

				unsigned int preLoopSamples = samplesRead;

				// Save the last sample per channel for crossfade continuity
				float lastEndSamples[MAX_CHANNELS];
				for (unsigned int ch = 0; ch < voice->mChannels; ch++)
				{
					lastEndSamples[ch] = (preLoopSamples > 0) ? channelBuffer[ch * alignedBufferSize + preLoopSamples - 1] : 0.0f;
				}

				while (samplesRead < samplesToRead && voice->seek(voice->mLoopPoint, channelBuffer, samplesToRead * voice->mChannels) == SO_NO_ERROR)
				{
					voice->mLoopCount++;
					unsigned int remaining = samplesToRead - samplesRead;

					// Check available space for loop temp data.
					// getAudio implementations use aSamplesToRead as the channel stride,
					// so we need remaining * channels floats of space.
					unsigned int mainBufferEnd = alignedBufferSize * voice->mChannels;
					unsigned int tempSpaceAvailable = scratchSize - mainBufferEnd;
					if (remaining * voice->mChannels > tempSpaceAvailable)
					{
						remaining = tempSpaceAvailable / voice->mChannels;
					}

					if (remaining == 0)
						break;

					float *loopTempBuffer = channelBuffer + mainBufferEnd;
					unsigned int loopSamples = voice->getAudio(loopTempBuffer, remaining, remaining);

					if (loopSamples == 0)
						break;

					// Copy loop samples to channelBuffer (using remaining as source stride
					// since that's what getAudio uses internally)
					for (size_t ch = 0; ch < voice->mChannels; ch++)
					{
						memcpy(channelBuffer + ch * alignedBufferSize + samplesRead, loopTempBuffer + ch * remaining, loopSamples * sizeof(float));
					}

					// Apply crossfade at the loop boundary (only on the first iteration after seeking)
					if (samplesRead == preLoopSamples && preLoopSamples > 0 && loopSamples > 0)
					{
						unsigned int crossfadeLen = (loopSamples < LOOP_CROSSFADE_SAMPLES) ? loopSamples : LOOP_CROSSFADE_SAMPLES;

						for (size_t ch = 0; ch < voice->mChannels; ch++)
						{
							float lastEndSample = lastEndSamples[ch];
							float *loopStart = channelBuffer + ch * alignedBufferSize + preLoopSamples;

							for (unsigned int i = 0; i < crossfadeLen; i++)
							{
								// t=0 at boundary (output = lastEndSample), t=1 at end (output = loopStart)
								float t = (float)(i + 1) / (float)(crossfadeLen + 1);
								// Fade from the last end sample into the loop samples
								loopStart[i] = lastEndSample * (1.0f - t) + loopStart[i] * t;
							}
						}
					}

					samplesRead += loopSamples;
				}
			}

			// Apply per-voice filters to source data before resampling
			if (samplesRead > 0)
			{
				for (unsigned int filterIdx = 0; filterIdx < FILTERS_PER_STREAM; filterIdx++)
				{
					if (voice->mFilter[filterIdx])
					{
						voice->mFilter[filterIdx]->filter(channelBuffer, samplesRead, alignedBufferSize, voice->mChannels, voice->mSamplerate, voice->mStreamTime);
					}
				}
			}

			// Copy from channel-interleaved scratch buffer to channel-separated resample buffers
			if (voice->mResampleBuffer != nullptr)
			{
				for (unsigned int ch = 0; ch < voice->mChannels; ch++)
				{
					float *srcChannel = channelBuffer + ch * alignedBufferSize;
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
	}

	voice->mResampleBufferFill += samplesRead;
	return voice->mResampleBufferFill - voice->mResampleBufferPos;
}

// High-precision resampling function that handles arbitrary sample rate ratios
unsigned int Soloud::resampleVoicePrecise_internal(AudioSourceInstance *voice,
                                                   float *outputBuffer,
                                                   unsigned int outputSamples,
                                                   unsigned int outputStride,
                                                   double outputSampleRate,
                                                   unsigned int resampler,
                                                   float *scratchBuffer,
                                                   unsigned int scratchSize)
{
	using namespace ResamplingConstants;

	if (outputSamples == 0 || !voice || voice->mResampleBuffer == nullptr)
		return 0;

	// Calculate step size: how much we advance in source per output sample
	double stepSize = voice->mSamplerate / outputSampleRate;

	// Determine lookahead requirements based on interpolation algorithm
	unsigned int lookaheadSamples = LINEAR_LOOKAHEAD; // Default fallback
	switch (resampler)
	{
	case RESAMPLER_POINT:
		lookaheadSamples = POINT_LOOKAHEAD;
		break;
	case RESAMPLER_LINEAR:
		lookaheadSamples = LINEAR_LOOKAHEAD;
		break;
	case RESAMPLER_CATMULLROM:
		lookaheadSamples = CATMULLROM_LOOKAHEAD;
		break;
	}

	// For chunked processing, we may not be able to produce all requested samples
	// if we run out of source data. This is normal and prevents excessive buffering.
	unsigned int samplesProduced = 0;
	unsigned int samplesToProcess = outputSamples;

	// Limit processing to one chunk worth of output to maintain low latency
	// This prevents excessive buffering for high sample rate ratios
	if (samplesToProcess > SAMPLE_GRANULARITY)
	{
		samplesToProcess = SAMPLE_GRANULARITY;
	}

	// Calculate how much input we need for this chunk of output
	unsigned int inputNeeded = (unsigned int)ceil(samplesToProcess * stepSize) + lookaheadSamples + LOOKAHEAD_SAFETY_MARGIN;
	unsigned int availableInput = ensureSourceData_internal(voice, inputNeeded, scratchBuffer, scratchSize);

	if (availableInput < lookaheadSamples)
		return 0; // Not enough input data available for interpolation

	unsigned int safeOutputCount = (unsigned int)ceil((availableInput - lookaheadSamples + 1.0 - voice->mPreciseSrcPosition) / stepSize);
	if (safeOutputCount > samplesToProcess)
		safeOutputCount = samplesToProcess;

	// If no output buffer provided (tick-only mode), just advance position
	if (!outputBuffer)
	{
		voice->mPreciseSrcPosition += safeOutputCount * stepSize;

		// Update buffer position for consumed integer samples
		unsigned int integralConsumed = (unsigned int)floor(voice->mPreciseSrcPosition);
		if (integralConsumed > 0 && integralConsumed <= availableInput)
		{
			voice->mResampleBufferPos += integralConsumed;
			voice->mPreciseSrcPosition -= integralConsumed;
		}
		return safeOutputCount;
	}

	if (safeOutputCount == 0)
		return 0;

	// Resample all channels using SIMD-optimized implementation
	// Create offset pointers for each channel (accounting for current read position)
	float *srcChannelsOffset[MAX_CHANNELS];
	for (unsigned int ch = 0; ch < voice->mChannels; ch++)
	{
		srcChannelsOffset[ch] = voice->mResampleBuffer[ch] + voice->mResampleBufferPos;
	}

	samplesProduced = safeOutputCount; // Always process this amount
	mMixer->resample_channels(srcChannelsOffset, outputBuffer, outputStride, voice->mChannels, safeOutputCount, voice->mPreciseSrcPosition, stepSize, availableInput,
	                          lookaheadSamples);

	// Update position tracking with precise accumulation
	double totalAdvance = samplesProduced * stepSize;
	voice->mPreciseSrcPosition += totalAdvance;

	// Update buffer position for consumed integer samples while preserving fractional part
	unsigned int integralConsumed = (unsigned int)floor(voice->mPreciseSrcPosition);
	if (integralConsumed > 0 && integralConsumed <= availableInput)
	{
		voice->mResampleBufferPos += integralConsumed;
		voice->mPreciseSrcPosition -= integralConsumed; // Keep fractional part for precise positioning
	}

	return samplesProduced;
}

namespace MixingConstants
{
// Scratch buffer allocation strategy
// Each voice needs space for: voice processing + temp reads + delay handling
constexpr size_t VOICE_SCRATCH_MULTIPLIER = MAX_CHANNELS;          // Voice processing buffer
constexpr size_t TEMP_READ_BUFFER_SIZE = SAMPLE_GRANULARITY * 4UL; // Temporary buffer for chunk reading

// Voice processing limits to prevent resource exhaustion
constexpr size_t MIN_VOICES_PER_BATCH = 1; // Always process at least one voice
} // namespace MixingConstants

// High-performance audio mixing with dynamic resampling and chunked processing
void Soloud::mixBus_internal(float *aBuffer, unsigned int aSamplesToRead, unsigned int aBufferSize, float *aScratch, unsigned int aBus, float aSamplerate,
                             unsigned int aChannels, unsigned int aResampler)
{
	using namespace MixingConstants;

	// Clear accumulation buffer - this is where all voices will be mixed together
	for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
	{
		for (unsigned int channel = 0; channel < aChannels; channel++)
		{
			aBuffer[sample + channel * aBufferSize] = 0;
		}
	}

	// Calculate scratch space allocation for optimal voice processing
	// Each voice needs separate buffers to avoid interference during parallel processing
	unsigned int voiceScratchSize = VOICE_SCRATCH_MULTIPLIER * aBufferSize;
	unsigned int tempScratchSize = TEMP_READ_BUFFER_SIZE;
	unsigned int delayScratchSize = aChannels * aBufferSize;
	unsigned int totalPerVoice = voiceScratchSize + tempScratchSize + delayScratchSize;

	// Align memory allocation to improve cache performance
	totalPerVoice = (totalPerVoice + CPU_MEMORY_ALIGNMENT_MASK()) & ~CPU_MEMORY_ALIGNMENT_MASK();

	// Calculate maximum number of voices we can process in parallel given memory constraints
	unsigned int maxVoicesInParallel = (mScratchSize * MAX_CHANNELS) / totalPerVoice;
	if (maxVoicesInParallel < MIN_VOICES_PER_BATCH)
		maxVoicesInParallel = MIN_VOICES_PER_BATCH;

	// Process each active voice using chunk-based approach for consistent latency
	for (unsigned int voiceIdx = 0; voiceIdx < mActiveVoiceCount; voiceIdx++)
	{
		AudioSourceInstance *voice = mVoice[mActiveVoice[voiceIdx]];
		if (!voice || voice->mBusHandle != aBus)
			continue;

		// Determine voice processing requirements
		bool isPaused = (voice->mFlags & AudioSourceInstance::PAUSED) != 0;
		bool isInaudible = (voice->mFlags & AudioSourceInstance::INAUDIBLE) != 0;
		bool mustTick = (voice->mFlags & AudioSourceInstance::INAUDIBLE_TICK) != 0;

		bool isAudible = !isPaused && !isInaudible;
		bool shouldProcess = isAudible || (!isPaused && mustTick);

		if (!shouldProcess)
			continue;

		unsigned int outputSamples = aSamplesToRead;
		unsigned int outputOffset = 0;

		// Handle delay samples: sounds can be scheduled to start later
		if (voice->mDelaySamples > 0)
		{
			if (voice->mDelaySamples >= aSamplesToRead)
			{
				// Entire buffer is still in delay period
				voice->mDelaySamples -= aSamplesToRead;
				continue;
			}
			else
			{
				// Partial delay: sound starts partway through buffer
				outputOffset = voice->mDelaySamples;
				outputSamples = aSamplesToRead - voice->mDelaySamples;
				voice->mDelaySamples = 0;
			}
		}

		if (outputSamples == 0)
			continue;

		// Allocate scratch space for this voice using round-robin allocation
		// This prevents memory fragmentation and ensures predictable performance
		unsigned int voiceSlot = voiceIdx % maxVoicesInParallel;
		float *voiceScratch = aScratch + (voiceSlot * totalPerVoice);
		float *tempScratch = voiceScratch + voiceScratchSize;
		float *delayScratch = tempScratch + tempScratchSize;

		// Adjust temporary scratch size based on available space
		unsigned int actualTempScratchSize = tempScratchSize;
		if (voiceSlot == maxVoicesInParallel - 1)
		{
			// Last slot gets any remaining space
			unsigned int remainingSpace = (mScratchSize * MAX_CHANNELS) - (voiceSlot * totalPerVoice + voiceScratchSize + delayScratchSize);
			if (remainingSpace < tempScratchSize)
			{
				actualTempScratchSize = remainingSpace;
			}
		}

		// Process in chunks to maintain consistent latency and prevent buffer overflow
		// Large sample rate differences might require multiple iterations
		unsigned int totalProduced = 0;
		while (totalProduced < outputSamples)
		{
			unsigned int remainingSamples = outputSamples - totalProduced;
			unsigned int chunkOutput = remainingSamples;

			// Limit chunk size to maintain low latency and predictable performance
			if (chunkOutput > SAMPLE_GRANULARITY)
				chunkOutput = SAMPLE_GRANULARITY;

			if (isAudible)
			{
				// Clear the chunk area in voice scratch buffer
				for (unsigned int ch = 0; ch < voice->mChannels; ch++)
				{
					memset(voiceScratch + ch * aBufferSize + totalProduced, 0, chunkOutput * sizeof(float));
				}

				// Resample this chunk from source sample rate to output sample rate
				unsigned int samplesProduced = resampleVoicePrecise_internal(voice,
				                                                             voiceScratch + totalProduced, // Offset into scratch buffer
				                                                             chunkOutput,
				                                                             aBufferSize,
				                                                             aSamplerate,
				                                                             aResampler,
				                                                             tempScratch,
				                                                             actualTempScratchSize);

				if (samplesProduced == 0)
				{
					// No more samples available from source
					break;
				}

				totalProduced += samplesProduced;

				// If we got fewer samples than requested, source is exhausted
				if (samplesProduced < chunkOutput)
				{
					break;
				}
			}
			else if (mustTick && !isPaused)
			{
				// Inaudible but must tick: advance voice position without producing audio
				// This keeps the voice synchronized for when it becomes audible again
				unsigned int samplesAdvanced =
				    resampleVoicePrecise_internal(voice, nullptr, chunkOutput, 0, aSamplerate, aResampler, tempScratch, actualTempScratchSize);

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
				// No delay offset: mix directly into output buffer
				mMixer->panAndExpand(voice, aBuffer, totalProduced, aBufferSize, voiceScratch, aChannels);
			}
			else
			{
				// Handle delay offset using temporary buffer
				memset(delayScratch, 0, aChannels * aBufferSize * sizeof(float));
				mMixer->panAndExpand(voice, delayScratch, totalProduced, aBufferSize, voiceScratch, aChannels);

				// Mix delayed content into output buffer at correct offset
				for (unsigned int ch = 0; ch < aChannels; ch++)
				{
					for (unsigned int sample = 0; sample < totalProduced; sample++)
					{
						if (sample + outputOffset < aSamplesToRead)
						{
							aBuffer[(sample + outputOffset) + ch * aBufferSize] += delayScratch[sample + ch * aBufferSize];
						}
					}
				}
			}
		}

		// Check if voice should be automatically stopped
		// Use small lookahead buffer check to account for chunked processing approach
		unsigned int remainingSamples = voice->mResampleBufferFill - voice->mResampleBufferPos;
		unsigned int lookaheadSamples = (aResampler == RESAMPLER_CATMULLROM) ? ResamplingConstants::CATMULLROM_LOOKAHEAD : ResamplingConstants::LINEAR_LOOKAHEAD;
		bool bufferNearEmpty = remainingSamples <= lookaheadSamples;

		// Stop voice if it has ended and auto-stop is enabled
		if (!(voice->mFlags & (AudioSourceInstance::LOOPING | AudioSourceInstance::DISABLE_AUTOSTOP)) && voice->hasEnded() && bufferNearEmpty)
		{
			stopVoice_internal(mActiveVoice[voiceIdx]);
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
	unsigned int *data = mActiveVoice.data() + mustlive;
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
	{
		static thread_local bool once = false;
		if (!once)
		{
			once = true;
			unsigned int u;
			u = _controlfp(0, 0);
			u = u & ~(_EM_INVALID | /*_EM_DENORMAL |*/ _EM_ZERODIVIDE | _EM_OVERFLOW /*| _EM_UNDERFLOW  | _EM_INEXACT*/);
			_controlfp(u, _MCW_EM);
		}
	}
#endif

	{
		// setFPUOptimizedRegs does some things that causes all math to consider really tiny values as zero, which
		// helps performance.
		static thread_local bool once = false;
		if (!once)
		{
			once = true;
			if (!(mFlags & NO_FPU_REGISTER_CHANGE))
				setFPUOptimizedRegs();
		}
	}

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
		const auto &currentVoice = mVoice[i];
		if (!currentVoice || (currentVoice->mFlags & AudioSourceInstance::PAUSED))
		{
			continue;
		}

		float volume[2];

		currentVoice->mActiveFader = 0;

		if (mGlobalVolumeFader.mActive > 0)
		{
			currentVoice->mActiveFader = 1;
		}

		currentVoice->mStreamTime += buffertime;
		currentVoice->mStreamPosition += (double)buffertime * (double)currentVoice->mOverallRelativePlaySpeed;

		// TODO: this is actually unstable, because mStreamTime depends on the relative
		// play speed.
		if (currentVoice->mRelativePlaySpeedFader.mActive > 0)
		{
			float speed = currentVoice->mRelativePlaySpeedFader.get(currentVoice->mStreamTime);
			setVoiceRelativePlaySpeed_internal(i, speed);
		}

		volume[0] = currentVoice->mOverallVolume;
		if (currentVoice->mVolumeFader.mActive > 0)
		{
			currentVoice->mSetVolume = currentVoice->mVolumeFader.get(currentVoice->mStreamTime);
			currentVoice->mActiveFader = 1;
			updateVoiceVolume_internal(i);
			mActiveVoiceDirty = true;
		}
		volume[1] = currentVoice->mOverallVolume;

		if (currentVoice->mPanFader.mActive > 0)
		{
			float pan = currentVoice->mPanFader.get(currentVoice->mStreamTime);
			setVoicePan_internal(i, pan);
			currentVoice->mActiveFader = 1;
		}

		if (currentVoice->mPauseScheduler.mActive)
		{
			currentVoice->mPauseScheduler.get(currentVoice->mStreamTime);
			if (currentVoice->mPauseScheduler.mActive == -1)
			{
				currentVoice->mPauseScheduler.mActive = 0;
				setVoicePause_internal(i, 1);
			}
		}

		if (currentVoice->mStopScheduler.mActive)
		{
			currentVoice->mStopScheduler.get(currentVoice->mStreamTime);
			if (currentVoice->mStopScheduler.mActive == -1)
			{
				currentVoice->mStopScheduler.mActive = 0;
				stopVoice_internal(i);
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
	mMixer->clip_samples(mOutputScratch.mData, mScratch.mData, aStride, mChannels, globalVolume[0], globalVolume[1], mPostClipScaler, mFlags & CLIP_ROUNDOFF);

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
					float absvol = (float)std::fabs(sample);
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
					float absvol = (float)std::fabs(sample);
					if (mVisualizationChannelVolume[j] < absvol)
						mVisualizationChannelVolume[j] = absvol;
					mVisualizationWaveData[i] += sample;
				}
			}
		}
	}
}

void Soloud::mix(void *aBuffer, unsigned int aSamples, SAMPLE_FORMAT aFormat)
{
	unsigned int stride = (aSamples + CPU_ALIGNMENT_MASK()) & ~CPU_ALIGNMENT_MASK();
	mix_internal(aSamples, stride);

	mMixer->interlace_samples(aBuffer, mScratch.mData, aSamples, stride, mChannels, aFormat);
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
