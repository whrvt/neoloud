/*
SoLoud audio engine - low-level mixing (internal) (SSE-optimized)
Copyright (c) 2013-2020 Jari Komppa
Copyright (c) 2025-2026 William Horvath

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

#include "soloud_cpu.h"

#if defined(SOLOUD_IS_X86) && !defined(DISABLE_SIMD)
#include "soloud.h"
#include "soloud_audiosource.h"
#include "soloud_ll_mixing.h"
#include "soloud_mixing_internal.h"

#ifdef _MSC_VER
#include <intrin.h>
#endif
#include <xmmintrin.h>

#include <climits> // _MAX
#include <cstring>

// NOLINTBEGIN(cppcoreguidelines-avoid-c-arrays, modernize-avoid-c-arrays, cppcoreguidelines-pro-type-member-init, hicpp-member-init)

namespace SoLoud::mixing
{

namespace
{
constexpr size_t SSE_OPTIMAL_CHUNK_SAMPLES = 4; // 4 floats per SSE register

// Used as a small temporary buffer for temporary SIMD register-sized data
template <typename T = float, size_t Size = SSE_OPTIMAL_CHUNK_SAMPLES>
struct alignas(SSE_ALIGNMENT_BYTES) TinyAlignedBuffer
{
	T data[Size];
};
} // namespace

void MixerSSE::interlace_samples(void *outputBuffer, const float *const rawBuffer, unsigned int aSamples, unsigned int stride, unsigned int aChannels,
                                 SAMPLE_FORMAT aFormat)
{
	using namespace mixing;
	unsigned int j = 0, c = 0;
	int i = 0;
	switch (aFormat)
	{
	case SAMPLE_FLOAT32: {
		float *buffer = static_cast<float *>(outputBuffer);
		for (j = 0; j < aChannels; j++)
		{
			c = j * stride;
			unsigned int outIdx = j;

			// Process 4 samples at a time with SSE
			unsigned int sampleQuads = aSamples / SSE_OPTIMAL_CHUNK_SAMPLES;
			unsigned int processedSamples = 0;

			for (unsigned int quad = 0; quad < sampleQuads; quad++)
			{
				__m128 samples = _mm_load_ps(&rawBuffer[c]);

				// Store samples to interleaved positions
				TinyAlignedBuffer temp;
				_mm_store_ps(temp.data, samples);

				for (unsigned int k = 0; k < SSE_OPTIMAL_CHUNK_SAMPLES; k++)
				{
					buffer[outIdx] = temp.data[k];
					outIdx += aChannels;
				}

				c += SSE_OPTIMAL_CHUNK_SAMPLES;
				processedSamples += SSE_OPTIMAL_CHUNK_SAMPLES;
			}

			// Handle remaining samples
			for (unsigned int remaining = processedSamples; remaining < aSamples; remaining++)
			{
				buffer[outIdx] = rawBuffer[c];
				outIdx += aChannels;
				c++;
			}
		}
	}
	break;

	case SAMPLE_UNSIGNED8: {
		unsigned char *buffer = static_cast<unsigned char *>(outputBuffer);
		__m128 scale = _mm_set_ps1(127.0f);
		__m128 offset = _mm_set_ps1(128.0f);
		__m128 zero = _mm_setzero_ps();
		__m128 max_val = _mm_set_ps1(255.0f);

		for (j = 0; j < aChannels; j++)
		{
			c = j * stride;
			unsigned int outIdx = j;

			unsigned int sampleQuads = aSamples / SSE_OPTIMAL_CHUNK_SAMPLES;
			unsigned int processedSamples = 0;

			for (unsigned int quad = 0; quad < sampleQuads; quad++)
			{
				__m128 samples = _mm_load_ps(&rawBuffer[c]);
				samples = _mm_mul_ps(samples, scale);
				samples = _mm_add_ps(samples, offset);

				// Clamp to [0, 255]
				samples = _mm_max_ps(samples, zero);
				samples = _mm_min_ps(samples, max_val);

				// Convert to integers and store
				__m128i int_samples = _mm_cvtps_epi32(samples);
				TinyAlignedBuffer<int> temp;

				_mm_store_si128((__m128i *)temp.data, int_samples);

				for (unsigned int k = 0; k < SSE_OPTIMAL_CHUNK_SAMPLES; k++)
				{
					buffer[outIdx] = (unsigned char)temp.data[k];
					outIdx += aChannels;
				}

				c += SSE_OPTIMAL_CHUNK_SAMPLES;
				processedSamples += SSE_OPTIMAL_CHUNK_SAMPLES;
			}

			// Handle remaining samples
			for (unsigned int remaining = processedSamples; remaining < aSamples; remaining++)
			{
				int sample = (int)(rawBuffer[c] * 127.0f + 128.0f);
				if (sample < 0)
					sample = 0;
				if (sample > 255)
					sample = 255;
				buffer[outIdx] = (unsigned char)sample;
				outIdx += aChannels;
				c++;
			}
		}
	}
	break;

	case SAMPLE_SIGNED16: {
		short *buffer = static_cast<short *>(outputBuffer);
		__m128 scale = _mm_set_ps1(0x7fff);

		for (j = 0; j < aChannels; j++)
		{
			c = j * stride;
			unsigned int outIdx = j;

			unsigned int sampleQuads = aSamples / SSE_OPTIMAL_CHUNK_SAMPLES;
			unsigned int processedSamples = 0;

			for (unsigned int quad = 0; quad < sampleQuads; quad++)
			{
				__m128 samples = _mm_load_ps(&rawBuffer[c]);
				samples = _mm_mul_ps(samples, scale);

				// Convert to int32 and store
				__m128i int_samples = _mm_cvtps_epi32(samples);
				TinyAlignedBuffer<int> temp;
				_mm_store_si128((__m128i *)temp.data, int_samples);

				for (unsigned int k = 0; k < SSE_OPTIMAL_CHUNK_SAMPLES; k++)
				{
					buffer[outIdx] = (short)temp.data[k];
					outIdx += aChannels;
				}

				c += SSE_OPTIMAL_CHUNK_SAMPLES;
				processedSamples += SSE_OPTIMAL_CHUNK_SAMPLES;
			}

			// Handle remaining samples
			for (unsigned int remaining = processedSamples; remaining < aSamples; remaining++)
			{
				buffer[outIdx] = (short)(rawBuffer[c] * 0x7fff);
				outIdx += aChannels;
				c++;
			}
		}
	}
	break;

	case SAMPLE_SIGNED24: {
		unsigned char *buffer = static_cast<unsigned char *>(outputBuffer);
		__m128 scale = _mm_set_ps1((float)(INT_MAX >> 8));
		__m128 min_val = _mm_set_ps1((float)(INT_MIN >> 8));
		__m128 max_val = _mm_set_ps1((float)(INT_MAX >> 8));

		for (j = 0; j < aChannels; j++)
		{
			c = j * stride;
			unsigned int outIdx = j * 3;

			unsigned int sampleQuads = aSamples / SSE_OPTIMAL_CHUNK_SAMPLES;
			unsigned int processedSamples = 0;

			for (unsigned int quad = 0; quad < sampleQuads; quad++)
			{
				__m128 samples = _mm_load_ps(&rawBuffer[c]);
				samples = _mm_mul_ps(samples, scale);

				// Clamp to 24-bit range
				samples = _mm_max_ps(samples, min_val);
				samples = _mm_min_ps(samples, max_val);

				// Convert to integers and store as 24-bit
				__m128i int_samples = _mm_cvtps_epi32(samples);
				TinyAlignedBuffer<int> temp;
				_mm_store_si128((__m128i *)temp.data, int_samples);

				for (unsigned int k = 0; k < SSE_OPTIMAL_CHUNK_SAMPLES; k++)
				{
					int sample = temp.data[k];
					buffer[outIdx] = (unsigned char)(sample & 0xff);
					buffer[outIdx + 1] = (unsigned char)((sample >> 8) & 0xff);
					buffer[outIdx + 2] = (unsigned char)((sample >> 16) & 0xff);
					outIdx += aChannels * 3;
				}

				c += SSE_OPTIMAL_CHUNK_SAMPLES;
				processedSamples += SSE_OPTIMAL_CHUNK_SAMPLES;
			}

			// Handle remaining samples
			for (unsigned int remaining = processedSamples; remaining < aSamples; remaining++)
			{
				int sample = (int)(rawBuffer[c] * (float)(INT_MAX >> 8));
				if (sample < (INT_MIN >> 8))
					sample = (INT_MIN >> 8);
				if (sample > (INT_MAX >> 8))
					sample = (INT_MAX >> 8);

				buffer[outIdx] = (unsigned char)(sample & 0xff);
				buffer[outIdx + 1] = (unsigned char)((sample >> 8) & 0xff);
				buffer[outIdx + 2] = (unsigned char)((sample >> 16) & 0xff);
				outIdx += aChannels * 3;
				c++;
			}
		}
	}
	break;

	case SAMPLE_SIGNED32: {
		int *buffer = static_cast<int *>(outputBuffer);
		for (j = 0; j < aChannels; j++)
		{
			c = j * stride;
			// unsigned int outIdx = j;

			// Use scalar for 32-bit to maintain precision
			for (i = j; i < aSamples * aChannels; i += aChannels)
			{
				double sample = (double)rawBuffer[c] * (double)INT_MAX;
				if (sample < (double)INT_MIN)
					sample = (double)INT_MIN;
				if (sample > (double)INT_MAX)
					sample = (double)INT_MAX;
				buffer[i] = (int)sample;
				c++;
			}
		}
	}
	break;
	}
}

void MixerSSE::clip_samples(const float *const aBuffer, float *aDestBuffer, unsigned int aSamples, unsigned int aChannels, float aVolume0, float aVolume1,
                            float aScaler, bool aRoundoff)
{
	using namespace ClippingConstants;

	float volumeDelta = (aVolume1 - aVolume0) / aSamples;
	float volume = aVolume0;
	unsigned int sampleQuads = (aSamples + 3) / 4; // Round up to process in SIMD groups

	if (aRoundoff)
	{
		// Roundoff clipping: smooth saturation curve instead of hard clipping
		__m128 negThreshold = _mm_load_ps1(&ROUNDOFF_NEG_THRESHOLD);
		__m128 posThreshold = _mm_load_ps1(&ROUNDOFF_POS_THRESHOLD);
		__m128 linearScale = _mm_load_ps1(&ROUNDOFF_LINEAR_SCALE);
		__m128 cubicScale = _mm_load_ps1(&ROUNDOFF_CUBIC_SCALE);
		__m128 negWall = _mm_load_ps1(&ROUNDOFF_NEG_WALL);
		__m128 posWall = _mm_load_ps1(&ROUNDOFF_POS_WALL);
		__m128 postScale = _mm_load_ps1(&aScaler);

		// Prepare volume ramp for SIMD processing
		TinyAlignedBuffer volumes;
		volumes.data[0] = volume;
		volumes.data[1] = volume + volumeDelta;
		volumes.data[2] = volume + volumeDelta * 2;
		volumes.data[3] = volume + volumeDelta * 3;
		volumeDelta *= SSE_OPTIMAL_CHUNK_SAMPLES;
		__m128 volDelta = _mm_load_ps1(&volumeDelta);

		unsigned int srcIdx = 0, dstIdx = 0;

		for (unsigned int channel = 0; channel < aChannels; channel++)
		{
			__m128 vol = _mm_load_ps(volumes.data);

			for (unsigned int quad = 0; quad < sampleQuads; quad++)
			{
				// Load 4 samples and apply volume ramp
				__m128 samples = _mm_load_ps(&aBuffer[srcIdx]);
				srcIdx += SSE_OPTIMAL_CHUNK_SAMPLES;
				samples = _mm_mul_ps(samples, vol);
				vol = _mm_add_ps(vol, volDelta);

				// Determine which samples are within the linear region
				__m128 aboveNegThreshold = _mm_cmpgt_ps(samples, negThreshold);
				__m128 belowPosThreshold = _mm_cmplt_ps(samples, posThreshold);

				// Apply roundoff curve: linear + cubic term
				__m128 linear = _mm_mul_ps(samples, linearScale);
				__m128 cubic = _mm_mul_ps(samples, samples);
				cubic = _mm_mul_ps(cubic, samples);
				cubic = _mm_mul_ps(cubic, cubicScale);
				samples = _mm_add_ps(linear, cubic);

				// Clamp values below negative threshold to negative wall
				__m128 lowMask = _mm_andnot_ps(aboveNegThreshold, negWall);
				__m128 lowKeep = _mm_and_ps(aboveNegThreshold, samples);
				samples = _mm_add_ps(lowMask, lowKeep);

				// Clamp values above positive threshold to positive wall
				__m128 highMask = _mm_andnot_ps(belowPosThreshold, posWall);
				__m128 highKeep = _mm_and_ps(belowPosThreshold, samples);
				samples = _mm_add_ps(highMask, highKeep);

				// Apply post-clip scaling and store
				samples = _mm_mul_ps(samples, postScale);
				_mm_store_ps(&aDestBuffer[dstIdx], samples);
				dstIdx += SSE_OPTIMAL_CHUNK_SAMPLES;
			}
		}
	}
	else
	{
		// Hard clipping: simple min/max bounds
		__m128 negBound = _mm_load_ps1(&HARD_CLIP_MIN);
		__m128 posBound = _mm_load_ps1(&HARD_CLIP_MAX);
		__m128 postScale = _mm_load_ps1(&aScaler);

		// Prepare volume ramp for SIMD processing
		TinyAlignedBuffer volumes;
		volumes.data[0] = volume;
		volumes.data[1] = volume + volumeDelta;
		volumes.data[2] = volume + volumeDelta * 2;
		volumes.data[3] = volume + volumeDelta * 3;
		volumeDelta *= SSE_OPTIMAL_CHUNK_SAMPLES;
		__m128 volDelta = _mm_load_ps1(&volumeDelta);

		unsigned int srcIdx = 0, dstIdx = 0;

		for (unsigned int channel = 0; channel < aChannels; channel++)
		{
			__m128 vol = _mm_load_ps(volumes.data);

			for (unsigned int quad = 0; quad < sampleQuads; quad++)
			{
				// Load 4 samples and apply volume ramp
				__m128 samples = _mm_load_ps(&aBuffer[srcIdx]);
				srcIdx += SSE_OPTIMAL_CHUNK_SAMPLES;
				samples = _mm_mul_ps(samples, vol);
				vol = _mm_add_ps(vol, volDelta);

				// Hard clipping to [-1, 1] range
				samples = _mm_max_ps(samples, negBound);
				samples = _mm_min_ps(samples, posBound);

				// Apply post-clip scaling and store
				samples = _mm_mul_ps(samples, postScale);
				_mm_store_ps(&aDestBuffer[dstIdx], samples);
				dstIdx += SSE_OPTIMAL_CHUNK_SAMPLES;
			}
		}
	}
}

void MixerSSE::panAndExpand(AudioSourceInstance *aVoice, float *aBuffer, unsigned int aSamplesToRead, unsigned int aBufferSize, float *aScratch,
                            unsigned int aChannels)
{
	// Only stereo output from 1 or 2 channels is currently accelerated. Delegate to the base implementation otherwise.
	if (aChannels != 2 || aVoice->mChannels > 2)
	{
		return Mixer::panAndExpand(aVoice, aBuffer, aSamplesToRead, aBufferSize, aScratch, aChannels);
	}

	using namespace ChannelMixingConstants;
	using namespace ClippingConstants;

	SOLOUD_ASSERT(((size_t)aBuffer & SSE_ALIGNMENT_MASK) == 0);
	SOLOUD_ASSERT(((size_t)aScratch & SSE_ALIGNMENT_MASK) == 0);
	SOLOUD_ASSERT(((size_t)aBufferSize & SSE_ALIGNMENT_MASK) == 0);

	// Calculate volume ramping for smooth transitions
	float currentPan[MAX_CHANNELS]; // Current speaker volumes
	float targetPan[MAX_CHANNELS];  // Target speaker volumes
	float panDelta[MAX_CHANNELS];   // Volume change per sample

	for (unsigned int ch = 0; ch < aChannels; ch++)
	{
		currentPan[ch] = aVoice->mCurrentChannelVolume[ch];
		targetPan[ch] = aVoice->mChannelVolume[ch] * aVoice->mOverallVolume;
		panDelta[ch] = (targetPan[ch] - currentPan[ch]) / aSamplesToRead;
	}

	// Channel mapping and mixing logic
	// Stereo output only

	switch (aVoice->mChannels)
	{
	case 2: // 2.0 -> 2.0: direct mapping with volume ramping
	{
		unsigned int processedSamples = 0;
		unsigned int sampleQuads = aSamplesToRead / SSE_OPTIMAL_CHUNK_SAMPLES;

		if (sampleQuads > 0)
		{
			// Prepare volume ramps for SIMD processing
			TinyAlignedBuffer pan0, pan1;
			for (int i = 0; i < SSE_OPTIMAL_CHUNK_SAMPLES; i++)
			{
				pan0.data[i] = currentPan[0] + panDelta[0] * (i + 1);
				pan1.data[i] = currentPan[1] + panDelta[1] * (i + 1);
			}

			float panDelta0_4x = panDelta[0] * SSE_OPTIMAL_CHUNK_SAMPLES;
			float panDelta1_4x = panDelta[1] * SSE_OPTIMAL_CHUNK_SAMPLES;
			__m128 pan0Delta = _mm_load_ps1(&panDelta0_4x);
			__m128 pan1Delta = _mm_load_ps1(&panDelta1_4x);
			__m128 p0 = _mm_load_ps(pan0.data);
			__m128 p1 = _mm_load_ps(pan1.data);

			for (unsigned int quad = 0; quad < sampleQuads; quad++)
			{
				unsigned int idx = processedSamples;

				// Load source samples
				__m128 src0 = _mm_load_ps(aScratch + idx);
				__m128 src1 = _mm_load_ps(aScratch + idx + aBufferSize);

				// Apply volume ramps
				__m128 mixed0 = _mm_mul_ps(src0, p0);
				__m128 mixed1 = _mm_mul_ps(src1, p1);

				// Add to existing output
				__m128 out0 = _mm_load_ps(aBuffer + idx);
				__m128 out1 = _mm_load_ps(aBuffer + idx + aBufferSize);
				out0 = _mm_add_ps(out0, mixed0);
				out1 = _mm_add_ps(out1, mixed1);

				// Store results
				_mm_store_ps(aBuffer + idx, out0);
				_mm_store_ps(aBuffer + idx + aBufferSize, out1);

				// Update volume ramps
				p0 = _mm_add_ps(p0, pan0Delta);
				p1 = _mm_add_ps(p1, pan1Delta);

				processedSamples += SSE_OPTIMAL_CHUNK_SAMPLES;
			}

			// Update current pan for remaining samples
			currentPan[0] += panDelta[0] * processedSamples;
			currentPan[1] += panDelta[1] * processedSamples;
		}

		// Process remaining samples
		for (unsigned int sample = processedSamples; sample < aSamplesToRead; sample++)
		{
			currentPan[0] += panDelta[0];
			currentPan[1] += panDelta[1];

			float leftSrc = aScratch[sample];
			float rightSrc = aScratch[aBufferSize + sample];

			aBuffer[sample] += leftSrc * currentPan[0];
			aBuffer[aBufferSize + sample] += rightSrc * currentPan[1];
		}
	}
	break;

	case 1: // 1.0 -> 2.0: distribute mono to stereo with panning
	{
		unsigned int processedSamples = 0;
		unsigned int sampleQuads = aSamplesToRead / SSE_OPTIMAL_CHUNK_SAMPLES;

		if (sampleQuads > 0)
		{
			// Prepare volume ramps for SIMD processing
			TinyAlignedBuffer pan0, pan1;
			for (int i = 0; i < SSE_OPTIMAL_CHUNK_SAMPLES; i++)
			{
				pan0.data[i] = currentPan[0] + panDelta[0] * (i + 1);
				pan1.data[i] = currentPan[1] + panDelta[1] * (i + 1);
			}

			float panDelta0_4x = panDelta[0] * SSE_OPTIMAL_CHUNK_SAMPLES;
			float panDelta1_4x = panDelta[1] * SSE_OPTIMAL_CHUNK_SAMPLES;
			__m128 pan0Delta = _mm_load_ps1(&panDelta0_4x);
			__m128 pan1Delta = _mm_load_ps1(&panDelta1_4x);
			__m128 p0 = _mm_load_ps(pan0.data);
			__m128 p1 = _mm_load_ps(pan1.data);

			for (unsigned int quad = 0; quad < sampleQuads; quad++)
			{
				unsigned int idx = processedSamples;

				// Load mono source
				__m128 src = _mm_load_ps(aScratch + idx);

				// Apply volume ramps to create stereo
				__m128 mixed0 = _mm_mul_ps(src, p0);
				__m128 mixed1 = _mm_mul_ps(src, p1);

				// Add to existing output
				__m128 out0 = _mm_load_ps(aBuffer + idx);
				__m128 out1 = _mm_load_ps(aBuffer + idx + aBufferSize);
				out0 = _mm_add_ps(out0, mixed0);
				out1 = _mm_add_ps(out1, mixed1);

				// Store results
				_mm_store_ps(aBuffer + idx, out0);
				_mm_store_ps(aBuffer + idx + aBufferSize, out1);

				// Update volume ramps
				p0 = _mm_add_ps(p0, pan0Delta);
				p1 = _mm_add_ps(p1, pan1Delta);

				processedSamples += SSE_OPTIMAL_CHUNK_SAMPLES;
			}

			// Update current pan for remaining samples
			currentPan[0] += panDelta[0] * processedSamples;
			currentPan[1] += panDelta[1] * processedSamples;
		}

		// Process remaining samples
		for (unsigned int sample = processedSamples; sample < aSamplesToRead; sample++)
		{
			currentPan[0] += panDelta[0];
			currentPan[1] += panDelta[1];

			float monoSrc = aScratch[sample];

			aBuffer[sample] += monoSrc * currentPan[0];
			aBuffer[aBufferSize + sample] += monoSrc * currentPan[1];
		}
	}
	break;
	}

	// Update voice state with final channel volumes
	for (unsigned int ch = 0; ch < aChannels; ch++)
	{
		aVoice->mCurrentChannelVolume[ch] = targetPan[ch];
	}
}

void MixerSSE::resample_channels(float **srcChannels, float *outputBuffer, unsigned int outputStride, unsigned int numChannels, unsigned int outputSamples,
                                 double srcPosition, double stepSize, unsigned int availableInput, unsigned int lookahead)
{
	using namespace ResamplingConstants;

	SOLOUD_ASSERT(srcChannels != nullptr);
	SOLOUD_ASSERT(outputBuffer != nullptr);
	SOLOUD_ASSERT(numChannels > 0);
	SOLOUD_ASSERT(availableInput >= lookahead);

	for (unsigned int ch = 0; ch < numChannels; ch++)
	{
		float *src = srcChannels[ch];
		float *dst = outputBuffer + ch * outputStride;

		unsigned int i = 0;

		// process 4 samples at a time with manual gathering
		double srcPos = srcPosition;
		unsigned int simdCount = outputSamples & ~3u;

		for (; i < simdCount; i += 4)
		{
			double pos0 = srcPos;
			double pos1 = srcPos + stepSize;
			double pos2 = srcPos + 2.0 * stepSize;
			double pos3 = srcPos + 3.0 * stepSize;

			unsigned int idx0 = (unsigned int)pos0;
			unsigned int idx1 = (unsigned int)pos1;
			unsigned int idx2 = (unsigned int)pos2;
			unsigned int idx3 = (unsigned int)pos3;

			float frac0 = (float)(pos0 - idx0);
			float frac1 = (float)(pos1 - idx1);
			float frac2 = (float)(pos2 - idx2);
			float frac3 = (float)(pos3 - idx3);

			__m128 frac = _mm_set_ps(frac3, frac2, frac1, frac0);
			__m128 result;

			if (lookahead == POINT_LOOKAHEAD)
			{
				result = _mm_set_ps(src[idx3], src[idx2], src[idx1], src[idx0]);
			}
			else if (lookahead == LINEAR_LOOKAHEAD)
			{
				unsigned int maxI = availableInput - 1;
				__m128 s0 = _mm_set_ps(src[idx3], src[idx2], src[idx1], src[idx0]);
				__m128 s1 = _mm_set_ps(src[(idx3 + 1 <= maxI) ? idx3 + 1 : idx3],
				                       src[(idx2 + 1 <= maxI) ? idx2 + 1 : idx2],
				                       src[(idx1 + 1 <= maxI) ? idx1 + 1 : idx1],
				                       src[(idx0 + 1 <= maxI) ? idx0 + 1 : idx0]);

				__m128 diff = _mm_sub_ps(s1, s0);
				result = _mm_add_ps(s0, _mm_mul_ps(diff, frac));
			}
			else // CATMULLROM_LOOKAHEAD
			{
				unsigned int maxI = availableInput - 1;

#define CLAMP_IDX_M1(idx) ((idx) >= 1 ? (idx) - 1 : (idx))
#define CLAMP_IDX_P1(idx) ((idx) + 1 <= maxI ? (idx) + 1 : maxI)
#define CLAMP_IDX_P2(idx) ((idx) + 2 <= maxI ? (idx) + 2 : maxI)

				__m128 p0 = _mm_set_ps(src[CLAMP_IDX_M1(idx3)], src[CLAMP_IDX_M1(idx2)], src[CLAMP_IDX_M1(idx1)], src[CLAMP_IDX_M1(idx0)]);
				__m128 p1 = _mm_set_ps(src[idx3], src[idx2], src[idx1], src[idx0]);
				__m128 p2 = _mm_set_ps(src[CLAMP_IDX_P1(idx3)], src[CLAMP_IDX_P1(idx2)], src[CLAMP_IDX_P1(idx1)], src[CLAMP_IDX_P1(idx0)]);
				__m128 p3 = _mm_set_ps(src[CLAMP_IDX_P2(idx3)], src[CLAMP_IDX_P2(idx2)], src[CLAMP_IDX_P2(idx1)], src[CLAMP_IDX_P2(idx0)]);

#undef CLAMP_IDX_M1
#undef CLAMP_IDX_P1
#undef CLAMP_IDX_P2

				__m128 t = frac;
				__m128 t2 = _mm_mul_ps(t, t);
				__m128 t3 = _mm_mul_ps(t2, t);

				__m128 cr_scale_sse = _mm_set1_ps(CATMULLROM_SCALE);
				__m128 cr_two_sse = _mm_set1_ps(CATMULLROM_LINEAR_COEFF);

				__m128 linear = _mm_mul_ps(cr_two_sse, p1);
				__m128 t_coeff = _mm_sub_ps(p2, p0);

				__m128 t2_coeff = _mm_mul_ps(_mm_set1_ps(CATMULLROM_QUAD_COEFFS[0]), p0);
				t2_coeff = _mm_add_ps(t2_coeff, _mm_mul_ps(_mm_set1_ps(CATMULLROM_QUAD_COEFFS[1]), p1));
				t2_coeff = _mm_add_ps(t2_coeff, _mm_mul_ps(_mm_set1_ps(CATMULLROM_QUAD_COEFFS[2]), p2));
				t2_coeff = _mm_add_ps(t2_coeff, _mm_mul_ps(_mm_set1_ps(CATMULLROM_QUAD_COEFFS[3]), p3));

				__m128 t3_coeff = _mm_mul_ps(_mm_set1_ps(CATMULLROM_CUBIC_COEFFS[0]), p0);
				t3_coeff = _mm_add_ps(t3_coeff, _mm_mul_ps(_mm_set1_ps(CATMULLROM_CUBIC_COEFFS[1]), p1));
				t3_coeff = _mm_add_ps(t3_coeff, _mm_mul_ps(_mm_set1_ps(CATMULLROM_CUBIC_COEFFS[2]), p2));
				t3_coeff = _mm_add_ps(t3_coeff, _mm_mul_ps(_mm_set1_ps(CATMULLROM_CUBIC_COEFFS[3]), p3));

				result = _mm_add_ps(linear, _mm_mul_ps(t_coeff, t));
				result = _mm_add_ps(result, _mm_mul_ps(t2_coeff, t2));
				result = _mm_add_ps(result, _mm_mul_ps(t3_coeff, t3));
				result = _mm_mul_ps(result, cr_scale_sse);
			}

			_mm_storeu_ps(dst + i, result);
			srcPos += 4.0 * stepSize;
		}

		// handle remainder after SIMD (delegate to base)
		Mixer::handleResampleRemainder(i, src, dst, outputSamples, srcPosition, stepSize, availableInput, lookahead);
	}

	return;
}

} // namespace SoLoud::mixing

// NOLINTEND(cppcoreguidelines-avoid-c-arrays, modernize-avoid-c-arrays, cppcoreguidelines-pro-type-member-init, hicpp-member-init)
#endif
