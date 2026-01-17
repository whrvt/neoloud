/*
SoLoud audio engine - low-level mixing (internal) (AVX-optimized)
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
#include <immintrin.h>
#include <xmmintrin.h>

#include <climits> // _MAX
#include <cstring>

// NOLINTBEGIN(cppcoreguidelines-avoid-c-arrays, modernize-avoid-c-arrays, cppcoreguidelines-pro-type-member-init, hicpp-member-init)

namespace SoLoud::mixing
{
namespace
{
constexpr size_t AVX_OPTIMAL_CHUNK_SAMPLES = 8; // 8 floats per AVX2 register

// Used as a small temporary buffer for temporary SIMD register-sized data
template <typename T = float, size_t Size = AVX_OPTIMAL_CHUNK_SAMPLES>
struct alignas(AVX_ALIGNMENT_BYTES) TinyAlignedBuffer
{
	T data[Size];
};
} // namespace

void MixerAVX::interlace_samples(void *outputBuffer, const float *const rawBuffer, unsigned int aSamples, unsigned int stride, unsigned int aChannels,
                                 SAMPLE_FORMAT aFormat)
{
	using namespace mixing;
	unsigned int j = 0, c = 0;
	switch (aFormat)
	{
	case SAMPLE_FLOAT32: {
		float *buffer = static_cast<float *>(outputBuffer);
		for (j = 0; j < aChannels; j++)
		{
			c = j * stride;
			unsigned int outIdx = j;

			// Process 8 samples at a time with AVX2
			unsigned int sampleOctuples = aSamples / AVX_OPTIMAL_CHUNK_SAMPLES;
			unsigned int processedSamples = 0;

			for (unsigned int octuple = 0; octuple < sampleOctuples; octuple++)
			{
				__m256 samples = _mm256_load_ps(&rawBuffer[c]);

				// Store samples to interleaved positions
				TinyAlignedBuffer temp;
				_mm256_store_ps(temp.data, samples);

				for (unsigned int k = 0; k < AVX_OPTIMAL_CHUNK_SAMPLES; k++)
				{
					buffer[outIdx] = temp.data[k];
					outIdx += aChannels;
				}

				c += AVX_OPTIMAL_CHUNK_SAMPLES;
				processedSamples += AVX_OPTIMAL_CHUNK_SAMPLES;
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
		__m256 scale = _mm256_set1_ps(127.0f);
		__m256 offset = _mm256_set1_ps(128.0f);
		__m256 zero = _mm256_setzero_ps();
		__m256 max_val = _mm256_set1_ps(255.0f);

		for (j = 0; j < aChannels; j++)
		{
			c = j * stride;
			unsigned int outIdx = j;

			unsigned int sampleOctuples = aSamples / AVX_OPTIMAL_CHUNK_SAMPLES;
			unsigned int processedSamples = 0;

			for (unsigned int octuple = 0; octuple < sampleOctuples; octuple++)
			{
				__m256 samples = _mm256_load_ps(&rawBuffer[c]);
				samples = _mm256_mul_ps(samples, scale);
				samples = _mm256_add_ps(samples, offset);

				// Clamp to [0, 255]
				samples = _mm256_max_ps(samples, zero);
				samples = _mm256_min_ps(samples, max_val);

				// Convert to integers and store
				__m256i int_samples = _mm256_cvtps_epi32(samples);
				TinyAlignedBuffer<int> temp;
				_mm256_store_si256((__m256i *)temp.data, int_samples);

				for (unsigned int k = 0; k < AVX_OPTIMAL_CHUNK_SAMPLES; k++)
				{
					buffer[outIdx] = (unsigned char)temp.data[k];
					outIdx += aChannels;
				}

				c += AVX_OPTIMAL_CHUNK_SAMPLES;
				processedSamples += AVX_OPTIMAL_CHUNK_SAMPLES;
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
		__m256 scale = _mm256_set1_ps(0x7fff);

		for (j = 0; j < aChannels; j++)
		{
			c = j * stride;
			unsigned int outIdx = j;

			unsigned int sampleOctuples = aSamples / AVX_OPTIMAL_CHUNK_SAMPLES;
			unsigned int processedSamples = 0;

			for (unsigned int octuple = 0; octuple < sampleOctuples; octuple++)
			{
				__m256 samples = _mm256_load_ps(&rawBuffer[c]);
				samples = _mm256_mul_ps(samples, scale);

				// Convert to int32 and store
				__m256i int_samples = _mm256_cvtps_epi32(samples);
				TinyAlignedBuffer<int> temp;
				_mm256_store_si256((__m256i *)temp.data, int_samples);

				for (unsigned int k = 0; k < AVX_OPTIMAL_CHUNK_SAMPLES; k++)
				{
					buffer[outIdx] = (short)temp.data[k];
					outIdx += aChannels;
				}

				c += AVX_OPTIMAL_CHUNK_SAMPLES;
				processedSamples += AVX_OPTIMAL_CHUNK_SAMPLES;
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
		__m256 scale = _mm256_set1_ps((float)(INT_MAX >> 8));
		__m256 min_val = _mm256_set1_ps((float)(INT_MIN >> 8));
		__m256 max_val = _mm256_set1_ps((float)(INT_MAX >> 8));

		for (j = 0; j < aChannels; j++)
		{
			c = j * stride;
			unsigned int outIdx = j * 3;

			unsigned int sampleOctuples = aSamples / AVX_OPTIMAL_CHUNK_SAMPLES;
			unsigned int processedSamples = 0;

			for (unsigned int octuple = 0; octuple < sampleOctuples; octuple++)
			{
				__m256 samples = _mm256_load_ps(&rawBuffer[c]);
				samples = _mm256_mul_ps(samples, scale);

				// Clamp to 24-bit range
				samples = _mm256_max_ps(samples, min_val);
				samples = _mm256_min_ps(samples, max_val);

				// Convert to integers and store as 24-bit
				__m256i int_samples = _mm256_cvtps_epi32(samples);
				TinyAlignedBuffer<int> temp;
				_mm256_store_si256((__m256i *)temp.data, int_samples);

				for (unsigned int k = 0; k < AVX_OPTIMAL_CHUNK_SAMPLES; k++)
				{
					int sample = temp.data[k];
					buffer[outIdx] = (unsigned char)(sample & 0xff);
					buffer[outIdx + 1] = (unsigned char)((sample >> 8) & 0xff);
					buffer[outIdx + 2] = (unsigned char)((sample >> 16) & 0xff);
					outIdx += aChannels * 3;
				}

				c += AVX_OPTIMAL_CHUNK_SAMPLES;
				processedSamples += AVX_OPTIMAL_CHUNK_SAMPLES;
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
		__m256d scale = _mm256_set1_pd((double)INT_MAX);
		__m256d min_val = _mm256_set1_pd((double)INT_MIN);
		__m256d max_val = _mm256_set1_pd((double)INT_MAX);

		for (j = 0; j < aChannels; j++)
		{
			c = j * stride;
			unsigned int outIdx = j;

			// Process 4 samples at a time for 32-bit (double precision needed)
			unsigned int sampleQuads = aSamples / 4;
			unsigned int processedSamples = 0;

			for (unsigned int quad = 0; quad < sampleQuads; quad++)
			{
				// Load 4 floats and convert to double for precision
				__m128 samples_f = _mm_load_ps(&rawBuffer[c]);
				__m256d samples = _mm256_cvtps_pd(samples_f);
				samples = _mm256_mul_pd(samples, scale);

				// Clamp to 32-bit range
				samples = _mm256_max_pd(samples, min_val);
				samples = _mm256_min_pd(samples, max_val);

				// Convert back to int
				TinyAlignedBuffer<double, 4> temp;
				_mm256_store_pd(temp.data, samples);

				for (unsigned int k = 0; k < 4; k++)
				{
					buffer[outIdx] = (int)temp.data[k];
					outIdx += aChannels;
				}

				c += 4;
				processedSamples += 4;
			}

			// Handle remaining samples
			for (unsigned int remaining = processedSamples; remaining < aSamples; remaining++)
			{
				double sample = (double)rawBuffer[c] * (double)INT_MAX;
				if (sample < (double)INT_MIN)
					sample = (double)INT_MIN;
				if (sample > (double)INT_MAX)
					sample = (double)INT_MAX;
				buffer[outIdx] = (int)sample;
				outIdx += aChannels;
				c++;
			}
		}
	}
	break;
	}
}

void MixerAVX::clip_samples(const float *const aBuffer, float *aDestBuffer, unsigned int aSamples, unsigned int aChannels, float aVolume0, float aVolume1,
                            float aScaler, bool aRoundoff)
{
	using namespace ClippingConstants;

	float volumeDelta = (aVolume1 - aVolume0) / aSamples;
	float volume = aVolume0;
	unsigned int sampleOctuples = (aSamples + 7) / 8; // Round up to process in AVX2 groups

	if (aRoundoff)
	{
		// Roundoff clipping: smooth saturation curve instead of hard clipping
		__m256 negThreshold = _mm256_broadcast_ss(&ROUNDOFF_NEG_THRESHOLD);
		__m256 posThreshold = _mm256_broadcast_ss(&ROUNDOFF_POS_THRESHOLD);
		__m256 linearScale = _mm256_broadcast_ss(&ROUNDOFF_LINEAR_SCALE);
		__m256 cubicScale = _mm256_broadcast_ss(&ROUNDOFF_CUBIC_SCALE);
		__m256 negWall = _mm256_broadcast_ss(&ROUNDOFF_NEG_WALL);
		__m256 posWall = _mm256_broadcast_ss(&ROUNDOFF_POS_WALL);
		__m256 postScale = _mm256_broadcast_ss(&aScaler);

		// Prepare volume ramp for AVX2 processing
		TinyAlignedBuffer volumes;
		for (int i = 0; i < AVX_OPTIMAL_CHUNK_SAMPLES; i++)
		{
			volumes.data[i] = volume + volumeDelta * i;
		}
		volumeDelta *= AVX_OPTIMAL_CHUNK_SAMPLES;
		__m256 volDelta = _mm256_broadcast_ss(&volumeDelta);

		unsigned int srcIdx = 0, dstIdx = 0;

		for (unsigned int channel = 0; channel < aChannels; channel++)
		{
			__m256 vol = _mm256_load_ps(volumes.data);

			for (unsigned int octuple = 0; octuple < sampleOctuples; octuple++)
			{
				// Load 8 samples and apply volume ramp
				__m256 samples = _mm256_load_ps(&aBuffer[srcIdx]);
				srcIdx += AVX_OPTIMAL_CHUNK_SAMPLES;
				samples = _mm256_mul_ps(samples, vol);
				vol = _mm256_add_ps(vol, volDelta);

				// Determine which samples are within the linear region
				__m256 aboveNegThreshold = _mm256_cmp_ps(samples, negThreshold, _CMP_GT_OQ);
				__m256 belowPosThreshold = _mm256_cmp_ps(samples, posThreshold, _CMP_LT_OQ);

				// Apply roundoff curve: linear + cubic term
				__m256 linear = _mm256_mul_ps(samples, linearScale);
				__m256 cubic = _mm256_mul_ps(samples, samples);
				cubic = _mm256_mul_ps(cubic, samples);
				cubic = _mm256_mul_ps(cubic, cubicScale);
				samples = _mm256_add_ps(linear, cubic);

				// Clamp values below negative threshold to negative wall
				samples = _mm256_blendv_ps(negWall, samples, aboveNegThreshold);

				// Clamp values above positive threshold to positive wall
				samples = _mm256_blendv_ps(posWall, samples, belowPosThreshold);

				// Apply post-clip scaling and store
				samples = _mm256_mul_ps(samples, postScale);
				_mm256_store_ps(&aDestBuffer[dstIdx], samples);
				dstIdx += AVX_OPTIMAL_CHUNK_SAMPLES;
			}
		}
	}
	else
	{
		// Hard clipping: simple min/max bounds
		__m256 negBound = _mm256_broadcast_ss(&HARD_CLIP_MIN);
		__m256 posBound = _mm256_broadcast_ss(&HARD_CLIP_MAX);
		__m256 postScale = _mm256_broadcast_ss(&aScaler);

		// Prepare volume ramp for AVX2 processing
		TinyAlignedBuffer volumes;
		for (int i = 0; i < AVX_OPTIMAL_CHUNK_SAMPLES; i++)
		{
			volumes.data[i] = volume + volumeDelta * i;
		}
		volumeDelta *= AVX_OPTIMAL_CHUNK_SAMPLES;
		__m256 volDelta = _mm256_broadcast_ss(&volumeDelta);

		unsigned int srcIdx = 0, dstIdx = 0;

		for (unsigned int channel = 0; channel < aChannels; channel++)
		{
			__m256 vol = _mm256_load_ps(volumes.data);

			for (unsigned int octuple = 0; octuple < sampleOctuples; octuple++)
			{
				// Load 8 samples and apply volume ramp
				__m256 samples = _mm256_load_ps(&aBuffer[srcIdx]);
				srcIdx += AVX_OPTIMAL_CHUNK_SAMPLES;
				samples = _mm256_mul_ps(samples, vol);
				vol = _mm256_add_ps(vol, volDelta);

				// Hard clipping to [-1, 1] range
				samples = _mm256_max_ps(samples, negBound);
				samples = _mm256_min_ps(samples, posBound);

				// Apply post-clip scaling and store
				samples = _mm256_mul_ps(samples, postScale);
				_mm256_store_ps(&aDestBuffer[dstIdx], samples);
				dstIdx += AVX_OPTIMAL_CHUNK_SAMPLES;
			}
		}
	}
}

void MixerAVX::panAndExpand(AudioSourceInstance *aVoice, float *aBuffer, unsigned int aSamplesToRead, unsigned int aBufferSize, float *aScratch,
                            unsigned int aChannels)
{
	// Only stereo output from 1 or 2 channels is currently accelerated. Delegate to the base implementation otherwise.
	if (aChannels != 2 || aVoice->mChannels > 2)
	{
		return Mixer::panAndExpand(aVoice, aBuffer, aSamplesToRead, aBufferSize, aScratch, aChannels);
	}

	using namespace ChannelMixingConstants;
	using namespace ClippingConstants;
	SOLOUD_ASSERT(((size_t)aBuffer & AVX_ALIGNMENT_MASK) == 0);
	SOLOUD_ASSERT(((size_t)aScratch & AVX_ALIGNMENT_MASK) == 0);
	SOLOUD_ASSERT(((size_t)aBufferSize & AVX_ALIGNMENT_MASK) == 0);

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
		unsigned int sampleOctuples = aSamplesToRead / AVX_OPTIMAL_CHUNK_SAMPLES;

		if (sampleOctuples > 0)
		{
			// Prepare volume ramps for AVX2 processing
			TinyAlignedBuffer pan0, pan1;
			for (int i = 0; i < AVX_OPTIMAL_CHUNK_SAMPLES; i++)
			{
				pan0.data[i] = currentPan[0] + panDelta[0] * (i + 1);
				pan1.data[i] = currentPan[1] + panDelta[1] * (i + 1);
			}

			float panDelta0_8x = panDelta[0] * AVX_OPTIMAL_CHUNK_SAMPLES;
			float panDelta1_8x = panDelta[1] * AVX_OPTIMAL_CHUNK_SAMPLES;
			__m256 pan0Delta = _mm256_broadcast_ss(&panDelta0_8x);
			__m256 pan1Delta = _mm256_broadcast_ss(&panDelta1_8x);
			__m256 p0 = _mm256_load_ps(pan0.data);
			__m256 p1 = _mm256_load_ps(pan1.data);

			for (unsigned int octuple = 0; octuple < sampleOctuples; octuple++)
			{
				unsigned int idx = processedSamples;

				// Load source samples
				__m256 src0 = _mm256_load_ps(aScratch + idx);
				__m256 src1 = _mm256_load_ps(aScratch + idx + aBufferSize);

				// Apply volume ramps
				__m256 mixed0 = _mm256_mul_ps(src0, p0);
				__m256 mixed1 = _mm256_mul_ps(src1, p1);

				// Add to existing output
				__m256 out0 = _mm256_load_ps(aBuffer + idx);
				__m256 out1 = _mm256_load_ps(aBuffer + idx + aBufferSize);
				out0 = _mm256_add_ps(out0, mixed0);
				out1 = _mm256_add_ps(out1, mixed1);

				// Store results
				_mm256_store_ps(aBuffer + idx, out0);
				_mm256_store_ps(aBuffer + idx + aBufferSize, out1);

				// Update volume ramps
				p0 = _mm256_add_ps(p0, pan0Delta);
				p1 = _mm256_add_ps(p1, pan1Delta);

				processedSamples += AVX_OPTIMAL_CHUNK_SAMPLES;
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
		unsigned int sampleOctuples = aSamplesToRead / AVX_OPTIMAL_CHUNK_SAMPLES;

		if (sampleOctuples > 0)
		{
			// Prepare volume ramps for AVX2 processing
			TinyAlignedBuffer pan0, pan1;
			for (int i = 0; i < AVX_OPTIMAL_CHUNK_SAMPLES; i++)
			{
				pan0.data[i] = currentPan[0] + panDelta[0] * (i + 1);
				pan1.data[i] = currentPan[1] + panDelta[1] * (i + 1);
			}

			float panDelta0_8x = panDelta[0] * AVX_OPTIMAL_CHUNK_SAMPLES;
			float panDelta1_8x = panDelta[1] * AVX_OPTIMAL_CHUNK_SAMPLES;
			__m256 pan0Delta = _mm256_broadcast_ss(&panDelta0_8x);
			__m256 pan1Delta = _mm256_broadcast_ss(&panDelta1_8x);
			__m256 p0 = _mm256_load_ps(pan0.data);
			__m256 p1 = _mm256_load_ps(pan1.data);

			for (unsigned int octuple = 0; octuple < sampleOctuples; octuple++)
			{
				unsigned int idx = processedSamples;

				// Load mono source
				__m256 src = _mm256_load_ps(aScratch + idx);

				// Apply volume ramps to create stereo
				__m256 mixed0 = _mm256_mul_ps(src, p0);
				__m256 mixed1 = _mm256_mul_ps(src, p1);

				// Add to existing output
				__m256 out0 = _mm256_load_ps(aBuffer + idx);
				__m256 out1 = _mm256_load_ps(aBuffer + idx + aBufferSize);
				out0 = _mm256_add_ps(out0, mixed0);
				out1 = _mm256_add_ps(out1, mixed1);

				// Store results
				_mm256_store_ps(aBuffer + idx, out0);
				_mm256_store_ps(aBuffer + idx + aBufferSize, out1);

				// Update volume ramps
				p0 = _mm256_add_ps(p0, pan0Delta);
				p1 = _mm256_add_ps(p1, pan1Delta);

				processedSamples += AVX_OPTIMAL_CHUNK_SAMPLES;
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

void MixerAVX::resample_channels(float **srcChannels, float *outputBuffer, unsigned int outputStride, unsigned int numChannels, unsigned int outputSamples,
                                 double srcPosition, double stepSize, unsigned int availableInput, unsigned int lookahead)
{
	using namespace ResamplingConstants;

	SOLOUD_ASSERT(srcChannels != nullptr);
	SOLOUD_ASSERT(outputBuffer != nullptr);
	SOLOUD_ASSERT(numChannels > 0);
	SOLOUD_ASSERT(availableInput >= lookahead);

	// Constants for index clamping
	__m256i zero_i = _mm256_setzero_si256();
	__m256i maxIdx = _mm256_set1_epi32((int)(availableInput - 1));

	// Step values for position advancement
	__m256d step4 = _mm256_set1_pd(stepSize * 4.0);
	__m256d step8 = _mm256_set1_pd(stepSize * 8.0);

	// Catmull-Rom coefficients
	__m256 cr_scale = _mm256_set1_ps(CATMULLROM_SCALE);
	__m256 cr_two = _mm256_set1_ps(CATMULLROM_LINEAR_COEFF);
	__m256 cr_q0 = _mm256_set1_ps(CATMULLROM_QUAD_COEFFS[0]);
	__m256 cr_q1 = _mm256_set1_ps(CATMULLROM_QUAD_COEFFS[1]);
	__m256 cr_q2 = _mm256_set1_ps(CATMULLROM_QUAD_COEFFS[2]);
	__m256 cr_q3 = _mm256_set1_ps(CATMULLROM_QUAD_COEFFS[3]);
	__m256 cr_c0 = _mm256_set1_ps(CATMULLROM_CUBIC_COEFFS[0]);
	__m256 cr_c1 = _mm256_set1_ps(CATMULLROM_CUBIC_COEFFS[1]);
	__m256 cr_c2 = _mm256_set1_ps(CATMULLROM_CUBIC_COEFFS[2]);
	__m256 cr_c3 = _mm256_set1_ps(CATMULLROM_CUBIC_COEFFS[3]);

	for (unsigned int ch = 0; ch < numChannels; ch++)
	{
		float *src = srcChannels[ch];
		float *dst = outputBuffer + ch * outputStride;

		unsigned int i = 0;

		__m256d pos_lo = _mm256_set_pd(srcPosition + 3.0 * stepSize, srcPosition + 2.0 * stepSize, srcPosition + stepSize, srcPosition);
		__m256d pos_hi = _mm256_add_pd(pos_lo, step4);

		unsigned int simdCount = outputSamples & ~7u;
		for (; i < simdCount; i += 8)
		{
			__m256d floor_lo = _mm256_floor_pd(pos_lo);
			__m256d floor_hi = _mm256_floor_pd(pos_hi);

			__m128i idx_lo_128 = _mm256_cvttpd_epi32(floor_lo);
			__m128i idx_hi_128 = _mm256_cvttpd_epi32(floor_hi);
			__m256i idx = _mm256_set_m128i(idx_hi_128, idx_lo_128);

			__m256d frac_lo_d = _mm256_sub_pd(pos_lo, floor_lo);
			__m256d frac_hi_d = _mm256_sub_pd(pos_hi, floor_hi);
			__m128 frac_lo_f = _mm256_cvtpd_ps(frac_lo_d);
			__m128 frac_hi_f = _mm256_cvtpd_ps(frac_hi_d);
			__m256 frac = _mm256_set_m128(frac_hi_f, frac_lo_f);

			__m256 result;

			if (lookahead == POINT_LOOKAHEAD)
			{
				result = _mm256_i32gather_ps(src, idx, 4);
			}
			else if (lookahead == LINEAR_LOOKAHEAD)
			{
				__m256i idx_p1 = _mm256_add_epi32(idx, _mm256_set1_epi32(1));
				idx_p1 = _mm256_min_epi32(idx_p1, maxIdx);

				__m256 s0 = _mm256_i32gather_ps(src, idx, 4);
				__m256 s1 = _mm256_i32gather_ps(src, idx_p1, 4);

				__m256 diff = _mm256_sub_ps(s1, s0);
				result = _mm256_add_ps(s0, _mm256_mul_ps(diff, frac));
			}
			else // CATMULLROM_LOOKAHEAD
			{
				__m256i idx_m1 = _mm256_sub_epi32(idx, _mm256_set1_epi32(1));
				idx_m1 = _mm256_max_epi32(idx_m1, zero_i);
				__m256i idx_p1 = _mm256_add_epi32(idx, _mm256_set1_epi32(1));
				idx_p1 = _mm256_min_epi32(idx_p1, maxIdx);
				__m256i idx_p2 = _mm256_add_epi32(idx, _mm256_set1_epi32(2));
				idx_p2 = _mm256_min_epi32(idx_p2, maxIdx);

				__m256 p0 = _mm256_i32gather_ps(src, idx_m1, 4);
				__m256 p1 = _mm256_i32gather_ps(src, idx, 4);
				__m256 p2 = _mm256_i32gather_ps(src, idx_p1, 4);
				__m256 p3 = _mm256_i32gather_ps(src, idx_p2, 4);

				__m256 t = frac;
				__m256 t2 = _mm256_mul_ps(t, t);
				__m256 t3 = _mm256_mul_ps(t2, t);

				__m256 linear = _mm256_mul_ps(cr_two, p1);
				__m256 t_coeff = _mm256_sub_ps(p2, p0);

				__m256 t2_coeff = _mm256_mul_ps(cr_q0, p0);
				t2_coeff = _mm256_add_ps(t2_coeff, _mm256_mul_ps(cr_q1, p1));
				t2_coeff = _mm256_add_ps(t2_coeff, _mm256_mul_ps(cr_q2, p2));
				t2_coeff = _mm256_add_ps(t2_coeff, _mm256_mul_ps(cr_q3, p3));

				__m256 t3_coeff = _mm256_mul_ps(cr_c0, p0);
				t3_coeff = _mm256_add_ps(t3_coeff, _mm256_mul_ps(cr_c1, p1));
				t3_coeff = _mm256_add_ps(t3_coeff, _mm256_mul_ps(cr_c2, p2));
				t3_coeff = _mm256_add_ps(t3_coeff, _mm256_mul_ps(cr_c3, p3));

				result = _mm256_add_ps(linear, _mm256_mul_ps(t3_coeff, t3));
				result = _mm256_add_ps(result, _mm256_mul_ps(t2_coeff, t2));
				result = _mm256_add_ps(result, _mm256_mul_ps(t_coeff, t));
				result = _mm256_mul_ps(result, cr_scale);
			}

			_mm256_storeu_ps(dst + i, result);

			pos_lo = _mm256_add_pd(pos_lo, step8);
			pos_hi = _mm256_add_pd(pos_hi, step8);
		}

		// handle remainder after SIMD (delegate to base)
		Mixer::handleResampleRemainder(i, src, dst, outputSamples, srcPosition, stepSize, availableInput, lookahead);
	}

	return;
}

} // namespace SoLoud::mixing

// NOLINTEND(cppcoreguidelines-avoid-c-arrays, modernize-avoid-c-arrays, cppcoreguidelines-pro-type-member-init, hicpp-member-init)
#endif
