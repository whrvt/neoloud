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

#include "soloud_intrin.h"
#include "soloud.h"
#include "soloud_audiosource.h"

#include <climits> // _MAX
#include <cstring>

namespace SoLoud
{

// Roundoff clipping algorithm constants
// These values provide smooth saturation instead of hard clipping
namespace ClippingConstants
{
// Input threshold bounds for roundoff clipping
constexpr float ROUNDOFF_NEG_THRESHOLD = -1.65f;
constexpr float ROUNDOFF_POS_THRESHOLD = 1.65f;

// Output saturation limits for roundoff clipping
constexpr float ROUNDOFF_NEG_WALL = -0.9862875f;
constexpr float ROUNDOFF_POS_WALL = 0.9862875f;

// Roundoff curve coefficients: output = LINEAR_SCALE * input + CUBIC_SCALE * input^3
// This creates a smooth S-curve that approaches the walls asymptotically
constexpr float ROUNDOFF_LINEAR_SCALE = 0.87f;
constexpr float ROUNDOFF_CUBIC_SCALE = -0.1f;

// Hard clipping bounds
constexpr float HARD_CLIP_MIN = -1.0f;
constexpr float HARD_CLIP_MAX = 1.0f;
} // namespace ClippingConstants

// Channel mixing coefficients for downmixing operations
namespace ChannelMixingConstants
{
// Coefficients for mixing multi-channel content to fewer channels
// These values are chosen to preserve perceived loudness while avoiding clipping

constexpr float MIX_8_TO_2_SCALE = 0.2f; // 8->2: Mix 5 channels per output
constexpr float MIX_6_TO_2_SCALE = 0.3f; // 6->2: Mix 4 channels per output
constexpr float MIX_4_TO_2_SCALE = 0.5f; // 4->2: Mix 2 channels per output

constexpr float CENTER_SUB_MIX_SCALE = 0.7f; // Center + sub mixing level
constexpr float SURROUND_MIX_SCALE = 0.5f;   // Surround channel mixing level
constexpr float QUAD_MIX_SCALE = 0.25f;      // 4-channel average mixing
} // namespace ChannelMixingConstants

AlignedFloatBuffer::AlignedFloatBuffer()
{
	mBasePtr = nullptr;
	mData = nullptr;
	mFloats = 0;
}

result AlignedFloatBuffer::init(unsigned int aFloats)
{
	delete[] mBasePtr;
	mBasePtr = nullptr;
	mData = nullptr;
	mFloats = aFloats;

	mBasePtr = new unsigned char[aFloats * sizeof(float) + SIMD_ALIGNMENT_BYTES];
	if (mBasePtr == nullptr)
		return OUT_OF_MEMORY;
	mData = (float *)(((size_t)mBasePtr + SIMD_ALIGNMENT_MASK) & ~SIMD_ALIGNMENT_MASK);

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
    : mActualData()
{
	unsigned char *basePtr = &mActualData[0];
	mData = (float *)(((size_t)basePtr + SIMD_ALIGNMENT_MASK) & ~SIMD_ALIGNMENT_MASK);
}

void interlace_samples(void *outputBuffer, const float *const rawBuffer, unsigned int aSamples, unsigned int stride, unsigned int aChannels,
                       detail::SAMPLE_FORMAT aFormat)
{
	using namespace detail;
	unsigned int j = 0, c = 0;

#if defined(SOLOUD_AVX_INTRINSICS)
	switch (aFormat)
	{
	case SAMPLE_FLOAT32: {
		float *buffer = static_cast<float *>(outputBuffer);
		for (j = 0; j < aChannels; j++)
		{
			c = j * stride;
			unsigned int outIdx = j;

			// Process 8 samples at a time with AVX2
			unsigned int sampleOctuples = aSamples / OPTIMAL_CHUNK_SAMPLES;
			unsigned int processedSamples = 0;

			for (unsigned int octuple = 0; octuple < sampleOctuples; octuple++)
			{
				__m256 samples = _mm256_load_ps(&rawBuffer[c]);

				// Store samples to interleaved positions
#ifdef _MSC_VER
				__declspec(align(32)) float temp[OPTIMAL_CHUNK_SAMPLES];
#else
				float temp[OPTIMAL_CHUNK_SAMPLES] __attribute__((aligned(32)));
#endif
				_mm256_store_ps(temp, samples);

				for (unsigned int k = 0; k < OPTIMAL_CHUNK_SAMPLES; k++)
				{
					buffer[outIdx] = temp[k];
					outIdx += aChannels;
				}

				c += OPTIMAL_CHUNK_SAMPLES;
				processedSamples += OPTIMAL_CHUNK_SAMPLES;
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

			unsigned int sampleOctuples = aSamples / OPTIMAL_CHUNK_SAMPLES;
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
#ifdef _MSC_VER
				__declspec(align(32)) int temp[OPTIMAL_CHUNK_SAMPLES];
#else
				int temp[OPTIMAL_CHUNK_SAMPLES] __attribute__((aligned(32)));
#endif
				_mm256_store_si256((__m256i *)temp, int_samples);

				for (unsigned int k = 0; k < OPTIMAL_CHUNK_SAMPLES; k++)
				{
					buffer[outIdx] = (unsigned char)temp[k];
					outIdx += aChannels;
				}

				c += OPTIMAL_CHUNK_SAMPLES;
				processedSamples += OPTIMAL_CHUNK_SAMPLES;
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

			unsigned int sampleOctuples = aSamples / OPTIMAL_CHUNK_SAMPLES;
			unsigned int processedSamples = 0;

			for (unsigned int octuple = 0; octuple < sampleOctuples; octuple++)
			{
				__m256 samples = _mm256_load_ps(&rawBuffer[c]);
				samples = _mm256_mul_ps(samples, scale);

				// Convert to int32 and store
				__m256i int_samples = _mm256_cvtps_epi32(samples);
#ifdef _MSC_VER
				__declspec(align(32)) int temp[OPTIMAL_CHUNK_SAMPLES];
#else
				int temp[OPTIMAL_CHUNK_SAMPLES] __attribute__((aligned(32)));
#endif
				_mm256_store_si256((__m256i *)temp, int_samples);

				for (unsigned int k = 0; k < OPTIMAL_CHUNK_SAMPLES; k++)
				{
					buffer[outIdx] = (short)temp[k];
					outIdx += aChannels;
				}

				c += OPTIMAL_CHUNK_SAMPLES;
				processedSamples += OPTIMAL_CHUNK_SAMPLES;
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

			unsigned int sampleOctuples = aSamples / OPTIMAL_CHUNK_SAMPLES;
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
#ifdef _MSC_VER
				__declspec(align(32)) int temp[OPTIMAL_CHUNK_SAMPLES];
#else
				int temp[OPTIMAL_CHUNK_SAMPLES] __attribute__((aligned(32)));
#endif
				_mm256_store_si256((__m256i *)temp, int_samples);

				for (unsigned int k = 0; k < OPTIMAL_CHUNK_SAMPLES; k++)
				{
					int sample = temp[k];
					buffer[outIdx] = (unsigned char)(sample & 0xff);
					buffer[outIdx + 1] = (unsigned char)((sample >> 8) & 0xff);
					buffer[outIdx + 2] = (unsigned char)((sample >> 16) & 0xff);
					outIdx += aChannels * 3;
				}

				c += OPTIMAL_CHUNK_SAMPLES;
				processedSamples += OPTIMAL_CHUNK_SAMPLES;
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
#ifdef _MSC_VER
				__declspec(align(32)) double temp[4];
#else
				double temp[4] __attribute__((aligned(32)));
#endif
				_mm256_store_pd(temp, samples);

				for (unsigned int k = 0; k < 4; k++)
				{
					buffer[outIdx] = (int)temp[k];
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

#elif defined(SOLOUD_SSE_INTRINSICS)
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
			unsigned int sampleQuads = aSamples / OPTIMAL_CHUNK_SAMPLES;
			unsigned int processedSamples = 0;

			for (unsigned int quad = 0; quad < sampleQuads; quad++)
			{
				__m128 samples = _mm_load_ps(&rawBuffer[c]);

				// Store samples to interleaved positions
#ifdef _MSC_VER
				__declspec(align(16)) float temp[OPTIMAL_CHUNK_SAMPLES];
#else
				float temp[OPTIMAL_CHUNK_SAMPLES] __attribute__((aligned(16)));
#endif
				_mm_store_ps(temp, samples);

				for (unsigned int k = 0; k < OPTIMAL_CHUNK_SAMPLES; k++)
				{
					buffer[outIdx] = temp[k];
					outIdx += aChannels;
				}

				c += OPTIMAL_CHUNK_SAMPLES;
				processedSamples += OPTIMAL_CHUNK_SAMPLES;
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

			unsigned int sampleQuads = aSamples / OPTIMAL_CHUNK_SAMPLES;
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
#ifdef _MSC_VER
				__declspec(align(16)) int temp[OPTIMAL_CHUNK_SAMPLES];
#else
				int temp[OPTIMAL_CHUNK_SAMPLES] __attribute__((aligned(16)));
#endif
				_mm_store_si128((__m128i *)temp, int_samples);

				for (unsigned int k = 0; k < OPTIMAL_CHUNK_SAMPLES; k++)
				{
					buffer[outIdx] = (unsigned char)temp[k];
					outIdx += aChannels;
				}

				c += OPTIMAL_CHUNK_SAMPLES;
				processedSamples += OPTIMAL_CHUNK_SAMPLES;
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

			unsigned int sampleQuads = aSamples / OPTIMAL_CHUNK_SAMPLES;
			unsigned int processedSamples = 0;

			for (unsigned int quad = 0; quad < sampleQuads; quad++)
			{
				__m128 samples = _mm_load_ps(&rawBuffer[c]);
				samples = _mm_mul_ps(samples, scale);

				// Convert to int32 and store
				__m128i int_samples = _mm_cvtps_epi32(samples);
#ifdef _MSC_VER
				__declspec(align(16)) int temp[OPTIMAL_CHUNK_SAMPLES];
#else
				int temp[OPTIMAL_CHUNK_SAMPLES] __attribute__((aligned(16)));
#endif
				_mm_store_si128((__m128i *)temp, int_samples);

				for (unsigned int k = 0; k < OPTIMAL_CHUNK_SAMPLES; k++)
				{
					buffer[outIdx] = (short)temp[k];
					outIdx += aChannels;
				}

				c += OPTIMAL_CHUNK_SAMPLES;
				processedSamples += OPTIMAL_CHUNK_SAMPLES;
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

			unsigned int sampleQuads = aSamples / OPTIMAL_CHUNK_SAMPLES;
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
#ifdef _MSC_VER
				__declspec(align(16)) int temp[OPTIMAL_CHUNK_SAMPLES];
#else
				int temp[OPTIMAL_CHUNK_SAMPLES] __attribute__((aligned(16)));
#endif
				_mm_store_si128((__m128i *)temp, int_samples);

				for (unsigned int k = 0; k < OPTIMAL_CHUNK_SAMPLES; k++)
				{
					int sample = temp[k];
					buffer[outIdx] = (unsigned char)(sample & 0xff);
					buffer[outIdx + 1] = (unsigned char)((sample >> 8) & 0xff);
					buffer[outIdx + 2] = (unsigned char)((sample >> 16) & 0xff);
					outIdx += aChannels * 3;
				}

				c += OPTIMAL_CHUNK_SAMPLES;
				processedSamples += OPTIMAL_CHUNK_SAMPLES;
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
			unsigned int outIdx = j;

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

#else // Fallback scalar implementation
	int i = 0;
	switch (aFormat)
	{
	case SAMPLE_FLOAT32: {
		float *buffer = static_cast<float *>(outputBuffer);
		for (j = 0; j < aChannels; j++)
		{
			c = j * stride;
			for (i = j; i < aSamples * aChannels; i += aChannels)
			{
				buffer[i] = rawBuffer[c];
				c++;
			}
		}
	}
	break;

	case SAMPLE_UNSIGNED8: {
		unsigned char *buffer = static_cast<unsigned char *>(outputBuffer);
		for (j = 0; j < aChannels; j++)
		{
			c = j * stride;
			for (i = j; i < aSamples * aChannels; i += aChannels)
			{
				int sample = (int)(rawBuffer[c] * 127.0f + 128.0f);
				if (sample < 0)
					sample = 0;
				if (sample > 255)
					sample = 255;
				buffer[i] = (unsigned char)sample;
				c++;
			}
		}
	}
	break;

	case SAMPLE_SIGNED16: {
		short *buffer = static_cast<short *>(outputBuffer);
		for (j = 0; j < aChannels; j++)
		{
			c = j * stride;
			for (i = j; i < aSamples * aChannels; i += aChannels)
			{
				buffer[i] = (short)(rawBuffer[c] * 0x7fff);
				c++;
			}
		}
	}
	break;

	case SAMPLE_SIGNED24: {
		unsigned char *buffer = static_cast<unsigned char *>(outputBuffer);
		for (j = 0; j < aChannels; j++)
		{
			c = j * stride;
			for (i = j; i < aSamples * aChannels; i += aChannels)
			{
				int sample = (int)(rawBuffer[c] * (float)(INT_MAX >> 8));
				if (sample < (INT_MIN >> 8))
					sample = (INT_MIN >> 8);
				if (sample > (INT_MAX >> 8))
					sample = (INT_MAX >> 8);

				unsigned int destIdx = i * 3;
				buffer[destIdx] = (unsigned char)(sample & 0xff);
				buffer[destIdx + 1] = (unsigned char)((sample >> 8) & 0xff);
				buffer[destIdx + 2] = (unsigned char)((sample >> 16) & 0xff);
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
#endif
}

#if defined(SOLOUD_AVX_INTRINSICS)
void clip_samples(const float *const aBuffer, float *aDestBuffer, unsigned int aSamples, unsigned int aChannels, float aVolume0, float aVolume1, float aScaler,
                  bool aRoundoff)
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
		TinyAlignedFloatBuffer volumes;
		for (int i = 0; i < OPTIMAL_CHUNK_SAMPLES; i++)
		{
			volumes.mData[i] = volume + volumeDelta * i;
		}
		volumeDelta *= OPTIMAL_CHUNK_SAMPLES;
		__m256 volDelta = _mm256_broadcast_ss(&volumeDelta);

		unsigned int srcIdx = 0, dstIdx = 0;

		for (unsigned int channel = 0; channel < aChannels; channel++)
		{
			__m256 vol = _mm256_load_ps(volumes.mData);

			for (unsigned int octuple = 0; octuple < sampleOctuples; octuple++)
			{
				// Load 8 samples and apply volume ramp
				__m256 samples = _mm256_load_ps(&aBuffer[srcIdx]);
				srcIdx += OPTIMAL_CHUNK_SAMPLES;
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
				dstIdx += OPTIMAL_CHUNK_SAMPLES;
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
		TinyAlignedFloatBuffer volumes;
		for (int i = 0; i < OPTIMAL_CHUNK_SAMPLES; i++)
		{
			volumes.mData[i] = volume + volumeDelta * i;
		}
		volumeDelta *= OPTIMAL_CHUNK_SAMPLES;
		__m256 volDelta = _mm256_broadcast_ss(&volumeDelta);

		unsigned int srcIdx = 0, dstIdx = 0;

		for (unsigned int channel = 0; channel < aChannels; channel++)
		{
			__m256 vol = _mm256_load_ps(volumes.mData);

			for (unsigned int octuple = 0; octuple < sampleOctuples; octuple++)
			{
				// Load 8 samples and apply volume ramp
				__m256 samples = _mm256_load_ps(&aBuffer[srcIdx]);
				srcIdx += OPTIMAL_CHUNK_SAMPLES;
				samples = _mm256_mul_ps(samples, vol);
				vol = _mm256_add_ps(vol, volDelta);

				// Hard clipping to [-1, 1] range
				samples = _mm256_max_ps(samples, negBound);
				samples = _mm256_min_ps(samples, posBound);

				// Apply post-clip scaling and store
				samples = _mm256_mul_ps(samples, postScale);
				_mm256_store_ps(&aDestBuffer[dstIdx], samples);
				dstIdx += OPTIMAL_CHUNK_SAMPLES;
			}
		}
	}
}

#elif defined(SOLOUD_SSE_INTRINSICS)
void clip_samples(const float *const aBuffer, float *aDestBuffer, unsigned int aSamples, unsigned int aChannels, float aVolume0, float aVolume1, float aScaler,
                  bool aRoundoff)
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
		TinyAlignedFloatBuffer volumes;
		volumes.mData[0] = volume;
		volumes.mData[1] = volume + volumeDelta;
		volumes.mData[2] = volume + volumeDelta * 2;
		volumes.mData[3] = volume + volumeDelta * 3;
		volumeDelta *= OPTIMAL_CHUNK_SAMPLES;
		__m128 volDelta = _mm_load_ps1(&volumeDelta);

		unsigned int srcIdx = 0, dstIdx = 0;

		for (unsigned int channel = 0; channel < aChannels; channel++)
		{
			__m128 vol = _mm_load_ps(volumes.mData);

			for (unsigned int quad = 0; quad < sampleQuads; quad++)
			{
				// Load 4 samples and apply volume ramp
				__m128 samples = _mm_load_ps(&aBuffer[srcIdx]);
				srcIdx += OPTIMAL_CHUNK_SAMPLES;
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
				dstIdx += OPTIMAL_CHUNK_SAMPLES;
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
		TinyAlignedFloatBuffer volumes;
		volumes.mData[0] = volume;
		volumes.mData[1] = volume + volumeDelta;
		volumes.mData[2] = volume + volumeDelta * 2;
		volumes.mData[3] = volume + volumeDelta * 3;
		volumeDelta *= OPTIMAL_CHUNK_SAMPLES;
		__m128 volDelta = _mm_load_ps1(&volumeDelta);

		unsigned int srcIdx = 0, dstIdx = 0;

		for (unsigned int channel = 0; channel < aChannels; channel++)
		{
			__m128 vol = _mm_load_ps(volumes.mData);

			for (unsigned int quad = 0; quad < sampleQuads; quad++)
			{
				// Load 4 samples and apply volume ramp
				__m128 samples = _mm_load_ps(&aBuffer[srcIdx]);
				srcIdx += OPTIMAL_CHUNK_SAMPLES;
				samples = _mm_mul_ps(samples, vol);
				vol = _mm_add_ps(vol, volDelta);

				// Hard clipping to [-1, 1] range
				samples = _mm_max_ps(samples, negBound);
				samples = _mm_min_ps(samples, posBound);

				// Apply post-clip scaling and store
				samples = _mm_mul_ps(samples, postScale);
				_mm_store_ps(&aDestBuffer[dstIdx], samples);
				dstIdx += OPTIMAL_CHUNK_SAMPLES;
			}
		}
	}
}

#else // Fallback implementation without SIMD

void clip_samples(const float *const aBuffer, float *aDestBuffer, unsigned int aSamples, unsigned int aChannels, float aVolume0, float aVolume1, float aScaler,
                  bool aRoundoff)
{
	using namespace ClippingConstants;

	float volumeDelta = (aVolume1 - aVolume0) / aSamples;
	float volume = aVolume0;
	unsigned int sampleQuads = (aSamples + 3) / 4; // Process in groups of 4 for consistency

	if (aRoundoff)
	{
		unsigned int srcIdx = 0, dstIdx = 0;

		for (unsigned int channel = 0; channel < aChannels; channel++)
		{
			volume = aVolume0;

			for (unsigned int quad = 0; quad < sampleQuads; quad++)
			{
				// Process 4 samples at a time
				for (int sample = 0; sample < OPTIMAL_CHUNK_SAMPLES; sample++)
				{
					float input = aBuffer[srcIdx++] * volume;
					volume += volumeDelta;

					// Apply roundoff clipping curve
					float output;
					if (input <= ROUNDOFF_NEG_THRESHOLD)
					{
						output = ROUNDOFF_NEG_WALL;
					}
					else if (input >= ROUNDOFF_POS_THRESHOLD)
					{
						output = ROUNDOFF_POS_WALL;
					}
					else
					{
						// Smooth curve: linear + cubic term
						output = ROUNDOFF_LINEAR_SCALE * input + ROUNDOFF_CUBIC_SCALE * input * input * input;
					}

					aDestBuffer[dstIdx++] = output * aScaler;
				}
			}
		}
	}
	else
	{
		unsigned int srcIdx = 0, dstIdx = 0;

		for (unsigned int channel = 0; channel < aChannels; channel++)
		{
			volume = aVolume0;

			for (unsigned int quad = 0; quad < sampleQuads; quad++)
			{
				// Process 4 samples at a time
				for (int sample = 0; sample < OPTIMAL_CHUNK_SAMPLES; sample++)
				{
					float input = aBuffer[srcIdx++] * volume;
					volume += volumeDelta;

					// Hard clipping to [-1, 1] range
					float output = (input <= HARD_CLIP_MIN) ? HARD_CLIP_MIN : (input >= HARD_CLIP_MAX) ? HARD_CLIP_MAX : input;

					aDestBuffer[dstIdx++] = output * aScaler;
				}
			}
		}
	}
}
#endif

void panAndExpand(AudioSourceInstance *aVoice, float *aBuffer, unsigned int aSamplesToRead, unsigned int aBufferSize, float *aScratch, unsigned int aChannels)
{
	using namespace ChannelMixingConstants;
	using namespace ClippingConstants;
#if defined(SOLOUD_AVX_INTRINSICS) || defined(SOLOUD_SSE_INTRINSICS)
	SOLOUD_ASSERT(((size_t)aBuffer & SIMD_ALIGNMENT_MASK) == 0);
	SOLOUD_ASSERT(((size_t)aScratch & SIMD_ALIGNMENT_MASK) == 0);
	SOLOUD_ASSERT(((size_t)aBufferSize & SIMD_ALIGNMENT_MASK) == 0);
#endif

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
	switch (aChannels)
	{
	case 1: // Mono output: sum all input channels
	{
		for (unsigned int srcCh = 0; srcCh < aVoice->mChannels; srcCh++)
		{
			unsigned int srcOffset = srcCh * aBufferSize;
			currentPan[0] = aVoice->mCurrentChannelVolume[0];

			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				currentPan[0] += panDelta[0];
				aBuffer[sample] += aScratch[srcOffset + sample] * currentPan[0];
			}
		}
	}
	break;

	case 2: // Stereo output
		switch (aVoice->mChannels)
		{
		case 8: // 8.0 -> 2.0: mix left/right groups with center and sub contribution
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				currentPan[0] += panDelta[0];
				currentPan[1] += panDelta[1];

				// Load all 8 channels
				float ch[8];
				for (int i = 0; i < 8; i++)
				{
					ch[i] = aScratch[aBufferSize * i + sample];
				}

				// Mix: FL + C + LFE + SL + BL for left, FR + C + LFE + SR + BR for right
				aBuffer[sample] += MIX_8_TO_2_SCALE * (ch[0] + ch[2] + ch[3] + ch[4] + ch[6]) * currentPan[0];
				aBuffer[aBufferSize + sample] += MIX_8_TO_2_SCALE * (ch[1] + ch[2] + ch[3] + ch[5] + ch[7]) * currentPan[1];
			}
			break;

		case 6: // 5.1 -> 2.0: mix left/right groups with center and sub contribution
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				currentPan[0] += panDelta[0];
				currentPan[1] += panDelta[1];

				// Load all 6 channels: FL, FR, C, LFE, SL, SR
				float ch[6];
				for (int i = 0; i < 6; i++)
				{
					ch[i] = aScratch[aBufferSize * i + sample];
				}

				// Mix: FL + C + LFE + SL for left, FR + C + LFE + SR for right
				aBuffer[sample] += MIX_6_TO_2_SCALE * (ch[0] + ch[2] + ch[3] + ch[4]) * currentPan[0];
				aBuffer[aBufferSize + sample] += MIX_6_TO_2_SCALE * (ch[1] + ch[2] + ch[3] + ch[5]) * currentPan[1];
			}
			break;

		case 4: // 4.0 -> 2.0: sum front/back pairs
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				currentPan[0] += panDelta[0];
				currentPan[1] += panDelta[1];

				float frontLeft = aScratch[sample];
				float frontRight = aScratch[aBufferSize + sample];
				float backLeft = aScratch[aBufferSize * 2 + sample];
				float backRight = aScratch[aBufferSize * 3 + sample];

				aBuffer[sample] += MIX_4_TO_2_SCALE * (frontLeft + backLeft) * currentPan[0];
				aBuffer[aBufferSize + sample] += MIX_4_TO_2_SCALE * (frontRight + backRight) * currentPan[1];
			}
			break;

		case 2: // 2.0 -> 2.0: direct mapping with volume ramping
#if defined(SOLOUD_AVX_INTRINSICS)
		{
			unsigned int processedSamples = 0;
			unsigned int sampleOctuples = aSamplesToRead / OPTIMAL_CHUNK_SAMPLES;

			if (sampleOctuples > 0)
			{
				// Prepare volume ramps for AVX2 processing
				TinyAlignedFloatBuffer pan0, pan1;
				for (int i = 0; i < OPTIMAL_CHUNK_SAMPLES; i++)
				{
					pan0.mData[i] = currentPan[0] + panDelta[0] * (i + 1);
					pan1.mData[i] = currentPan[1] + panDelta[1] * (i + 1);
				}

				float panDelta0_8x = panDelta[0] * OPTIMAL_CHUNK_SAMPLES;
				float panDelta1_8x = panDelta[1] * OPTIMAL_CHUNK_SAMPLES;
				__m256 pan0Delta = _mm256_broadcast_ss(&panDelta0_8x);
				__m256 pan1Delta = _mm256_broadcast_ss(&panDelta1_8x);
				__m256 p0 = _mm256_load_ps(pan0.mData);
				__m256 p1 = _mm256_load_ps(pan1.mData);

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

					processedSamples += OPTIMAL_CHUNK_SAMPLES;
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
#elif defined(SOLOUD_SSE_INTRINSICS)
		{
			unsigned int processedSamples = 0;
			unsigned int sampleQuads = aSamplesToRead / OPTIMAL_CHUNK_SAMPLES;

			if (sampleQuads > 0)
			{
				// Prepare volume ramps for SIMD processing
				TinyAlignedFloatBuffer pan0, pan1;
				for (int i = 0; i < OPTIMAL_CHUNK_SAMPLES; i++)
				{
					pan0.mData[i] = currentPan[0] + panDelta[0] * (i + 1);
					pan1.mData[i] = currentPan[1] + panDelta[1] * (i + 1);
				}

				float panDelta0_4x = panDelta[0] * OPTIMAL_CHUNK_SAMPLES;
				float panDelta1_4x = panDelta[1] * OPTIMAL_CHUNK_SAMPLES;
				__m128 pan0Delta = _mm_load_ps1(&panDelta0_4x);
				__m128 pan1Delta = _mm_load_ps1(&panDelta1_4x);
				__m128 p0 = _mm_load_ps(pan0.mData);
				__m128 p1 = _mm_load_ps(pan1.mData);

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

					processedSamples += OPTIMAL_CHUNK_SAMPLES;
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
#else // Fallback implementation
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				currentPan[0] += panDelta[0];
				currentPan[1] += panDelta[1];

				float leftSrc = aScratch[sample];
				float rightSrc = aScratch[aBufferSize + sample];

				aBuffer[sample] += leftSrc * currentPan[0];
				aBuffer[aBufferSize + sample] += rightSrc * currentPan[1];
			}
#endif
		break;

		case 1: // 1.0 -> 2.0: distribute mono to stereo with panning
#if defined(SOLOUD_AVX_INTRINSICS)
		{
			unsigned int processedSamples = 0;
			unsigned int sampleOctuples = aSamplesToRead / OPTIMAL_CHUNK_SAMPLES;

			if (sampleOctuples > 0)
			{
				// Prepare volume ramps for AVX2 processing
				TinyAlignedFloatBuffer pan0, pan1;
				for (int i = 0; i < OPTIMAL_CHUNK_SAMPLES; i++)
				{
					pan0.mData[i] = currentPan[0] + panDelta[0] * (i + 1);
					pan1.mData[i] = currentPan[1] + panDelta[1] * (i + 1);
				}

				float panDelta0_8x = panDelta[0] * OPTIMAL_CHUNK_SAMPLES;
				float panDelta1_8x = panDelta[1] * OPTIMAL_CHUNK_SAMPLES;
				__m256 pan0Delta = _mm256_broadcast_ss(&panDelta0_8x);
				__m256 pan1Delta = _mm256_broadcast_ss(&panDelta1_8x);
				__m256 p0 = _mm256_load_ps(pan0.mData);
				__m256 p1 = _mm256_load_ps(pan1.mData);

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

					processedSamples += OPTIMAL_CHUNK_SAMPLES;
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
#elif defined(SOLOUD_SSE_INTRINSICS)
		{
			unsigned int processedSamples = 0;
			unsigned int sampleQuads = aSamplesToRead / OPTIMAL_CHUNK_SAMPLES;

			if (sampleQuads > 0)
			{
				// Prepare volume ramps for SIMD processing
				TinyAlignedFloatBuffer pan0, pan1;
				for (int i = 0; i < OPTIMAL_CHUNK_SAMPLES; i++)
				{
					pan0.mData[i] = currentPan[0] + panDelta[0] * (i + 1);
					pan1.mData[i] = currentPan[1] + panDelta[1] * (i + 1);
				}

				float panDelta0_4x = panDelta[0] * OPTIMAL_CHUNK_SAMPLES;
				float panDelta1_4x = panDelta[1] * OPTIMAL_CHUNK_SAMPLES;
				__m128 pan0Delta = _mm_load_ps1(&panDelta0_4x);
				__m128 pan1Delta = _mm_load_ps1(&panDelta1_4x);
				__m128 p0 = _mm_load_ps(pan0.mData);
				__m128 p1 = _mm_load_ps(pan1.mData);

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

					processedSamples += OPTIMAL_CHUNK_SAMPLES;
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
#else // Fallback implementation
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				currentPan[0] += panDelta[0];
				currentPan[1] += panDelta[1];

				float monoSrc = aScratch[sample];

				aBuffer[sample] += monoSrc * currentPan[0];
				aBuffer[aBufferSize + sample] += monoSrc * currentPan[1];
			}
#endif
		break;
		}
		break;

	case 4: // Quadraphonic output
		switch (aVoice->mChannels)
		{
		case 8: // 8.0 -> 4.0: group channels and add center/sub
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				for (int ch = 0; ch < 4; ch++)
				{
					currentPan[ch] += panDelta[ch];
				}

				float ch[8];
				for (int i = 0; i < 8; i++)
				{
					ch[i] = aScratch[aBufferSize * i + sample];
				}

				float centerSub = (ch[2] + ch[3]) * CENTER_SUB_MIX_SCALE;

				aBuffer[sample] += ch[0] * currentPan[0] + centerSub;
				aBuffer[aBufferSize + sample] += ch[1] * currentPan[1] + centerSub;
				aBuffer[aBufferSize * 2 + sample] += SURROUND_MIX_SCALE * (ch[4] + ch[6]) * currentPan[2];
				aBuffer[aBufferSize * 3 + sample] += SURROUND_MIX_SCALE * (ch[5] + ch[7]) * currentPan[3];
			}
			break;

		case 6: // 5.1 -> 4.0: map channels and add center/sub
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				for (int ch = 0; ch < 4; ch++)
				{
					currentPan[ch] += panDelta[ch];
				}

				float ch[6];
				for (int i = 0; i < 6; i++)
				{
					ch[i] = aScratch[aBufferSize * i + sample];
				}

				float centerSub = (ch[2] + ch[3]) * CENTER_SUB_MIX_SCALE;

				aBuffer[sample] += ch[0] * currentPan[0] + centerSub;
				aBuffer[aBufferSize + sample] += ch[1] * currentPan[1] + centerSub;
				aBuffer[aBufferSize * 2 + sample] += ch[4] * currentPan[2];
				aBuffer[aBufferSize * 3 + sample] += ch[5] * currentPan[3];
			}
			break;

		case 4: // 4.0 -> 4.0: direct mapping
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				for (int ch = 0; ch < 4; ch++)
				{
					currentPan[ch] += panDelta[ch];
					float src = aScratch[aBufferSize * ch + sample];
					aBuffer[aBufferSize * ch + sample] += src * currentPan[ch];
				}
			}
			break;

		case 2: // 2.0 -> 4.0: duplicate to front and back
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				for (int ch = 0; ch < 4; ch++)
				{
					currentPan[ch] += panDelta[ch];
				}

				float left = aScratch[sample];
				float right = aScratch[aBufferSize + sample];

				aBuffer[sample] += left * currentPan[0];
				aBuffer[aBufferSize + sample] += right * currentPan[1];
				aBuffer[aBufferSize * 2 + sample] += left * currentPan[2];
				aBuffer[aBufferSize * 3 + sample] += right * currentPan[3];
			}
			break;

		case 1: // 1.0 -> 4.0: distribute mono to all channels
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				for (int ch = 0; ch < 4; ch++)
				{
					currentPan[ch] += panDelta[ch];
				}

				float monoSrc = aScratch[sample];

				for (int ch = 0; ch < 4; ch++)
				{
					aBuffer[aBufferSize * ch + sample] += monoSrc * currentPan[ch];
				}
			}
			break;
		}
		break;

	case 6: // 5.1 Surround output
		switch (aVoice->mChannels)
		{
		case 8: // 8.0 -> 5.1: combine side and back channels
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				for (int ch = 0; ch < 6; ch++)
				{
					currentPan[ch] += panDelta[ch];
				}

				float ch[8];
				for (int i = 0; i < 8; i++)
				{
					ch[i] = aScratch[aBufferSize * i + sample];
				}

				// FL, FR, C, LFE, SL+BL, SR+BR
				aBuffer[sample] += ch[0] * currentPan[0];
				aBuffer[aBufferSize + sample] += ch[1] * currentPan[1];
				aBuffer[aBufferSize * 2 + sample] += ch[2] * currentPan[2];
				aBuffer[aBufferSize * 3 + sample] += ch[3] * currentPan[3];
				aBuffer[aBufferSize * 4 + sample] += SURROUND_MIX_SCALE * (ch[4] + ch[6]) * currentPan[4];
				aBuffer[aBufferSize * 5 + sample] += SURROUND_MIX_SCALE * (ch[5] + ch[7]) * currentPan[5];
			}
			break;

		case 6: // 5.1 -> 5.1: direct mapping
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				for (int ch = 0; ch < 6; ch++)
				{
					currentPan[ch] += panDelta[ch];
					float src = aScratch[aBufferSize * ch + sample];
					aBuffer[aBufferSize * ch + sample] += src * currentPan[ch];
				}
			}
			break;

		case 4: // 4.0 -> 5.1: create center/sub from front channels
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				for (int ch = 0; ch < 6; ch++)
				{
					currentPan[ch] += panDelta[ch];
				}

				float ch[4];
				for (int i = 0; i < 4; i++)
				{
					ch[i] = aScratch[aBufferSize * i + sample];
				}

				// FL, FR, C (mix), LFE (mix), SL, SR
				aBuffer[sample] += ch[0] * currentPan[0];
				aBuffer[aBufferSize + sample] += ch[1] * currentPan[1];
				aBuffer[aBufferSize * 2 + sample] += SURROUND_MIX_SCALE * (ch[0] + ch[1]) * currentPan[2];
				aBuffer[aBufferSize * 3 + sample] += QUAD_MIX_SCALE * (ch[0] + ch[1] + ch[2] + ch[3]) * currentPan[3];
				aBuffer[aBufferSize * 4 + sample] += ch[2] * currentPan[4];
				aBuffer[aBufferSize * 5 + sample] += ch[3] * currentPan[5];
			}
			break;

		case 2: // 2.0 -> 5.1: create surround content
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				for (int ch = 0; ch < 6; ch++)
				{
					currentPan[ch] += panDelta[ch];
				}

				float left = aScratch[sample];
				float right = aScratch[aBufferSize + sample];

				// FL, FR, C (mix), LFE (mix), SL, SR
				aBuffer[sample] += left * currentPan[0];
				aBuffer[aBufferSize + sample] += right * currentPan[1];
				aBuffer[aBufferSize * 2 + sample] += SURROUND_MIX_SCALE * (left + right) * currentPan[2];
				aBuffer[aBufferSize * 3 + sample] += SURROUND_MIX_SCALE * (left + right) * currentPan[3];
				aBuffer[aBufferSize * 4 + sample] += left * currentPan[4];
				aBuffer[aBufferSize * 5 + sample] += right * currentPan[5];
			}
			break;

		case 1: // 1.0 -> 5.1: distribute mono to all channels
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				for (int ch = 0; ch < 6; ch++)
				{
					currentPan[ch] += panDelta[ch];
				}

				float monoSrc = aScratch[sample];

				for (int ch = 0; ch < 6; ch++)
				{
					aBuffer[aBufferSize * ch + sample] += monoSrc * currentPan[ch];
				}
			}
			break;
		}
		break;

	case 8: // 7.1 Surround output
		switch (aVoice->mChannels)
		{
		case 8: // 8.0 -> 7.1: direct mapping
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				for (int ch = 0; ch < 8; ch++)
				{
					currentPan[ch] += panDelta[ch];
					float src = aScratch[aBufferSize * ch + sample];
					aBuffer[aBufferSize * ch + sample] += src * currentPan[ch];
				}
			}
			break;

		case 6: // 5.1 -> 7.1: create side channels from surrounds
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				for (int ch = 0; ch < 8; ch++)
				{
					currentPan[ch] += panDelta[ch];
				}

				float ch[6];
				for (int i = 0; i < 6; i++)
				{
					ch[i] = aScratch[aBufferSize * i + sample];
				}

				// FL, FR, C, LFE, SL+blend, SR+blend, SL, SR
				aBuffer[sample] += ch[0] * currentPan[0];
				aBuffer[aBufferSize + sample] += ch[1] * currentPan[1];
				aBuffer[aBufferSize * 2 + sample] += ch[2] * currentPan[2];
				aBuffer[aBufferSize * 3 + sample] += ch[3] * currentPan[3];
				aBuffer[aBufferSize * 4 + sample] += SURROUND_MIX_SCALE * (ch[4] + ch[0]) * currentPan[4];
				aBuffer[aBufferSize * 5 + sample] += SURROUND_MIX_SCALE * (ch[5] + ch[1]) * currentPan[5];
				aBuffer[aBufferSize * 6 + sample] += ch[4] * currentPan[6];
				aBuffer[aBufferSize * 7 + sample] += ch[5] * currentPan[7];
			}
			break;

		case 4: // 4.0 -> 7.1: create center/sub and side channels
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				for (int ch = 0; ch < 8; ch++)
				{
					currentPan[ch] += panDelta[ch];
				}

				float ch[4];
				for (int i = 0; i < 4; i++)
				{
					ch[i] = aScratch[aBufferSize * i + sample];
				}

				// FL, FR, C (mix), LFE (mix), SL+FL, SR+FR, SL, SR
				aBuffer[sample] += ch[0] * currentPan[0];
				aBuffer[aBufferSize + sample] += ch[1] * currentPan[1];
				aBuffer[aBufferSize * 2 + sample] += SURROUND_MIX_SCALE * (ch[0] + ch[1]) * currentPan[2];
				aBuffer[aBufferSize * 3 + sample] += QUAD_MIX_SCALE * (ch[0] + ch[1] + ch[2] + ch[3]) * currentPan[3];
				aBuffer[aBufferSize * 4 + sample] += SURROUND_MIX_SCALE * (ch[0] + ch[2]) * currentPan[4];
				aBuffer[aBufferSize * 5 + sample] += SURROUND_MIX_SCALE * (ch[1] + ch[3]) * currentPan[5];
				aBuffer[aBufferSize * 6 + sample] += ch[2] * currentPan[4];
				aBuffer[aBufferSize * 7 + sample] += ch[3] * currentPan[5];
			}
			break;

		case 2: // 2.0 -> 7.1: create full surround field
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				for (int ch = 0; ch < 8; ch++)
				{
					currentPan[ch] += panDelta[ch];
				}

				float left = aScratch[sample];
				float right = aScratch[aBufferSize + sample];

				// FL, FR, C (mix), LFE (mix), SL, SR, BL, BR
				aBuffer[sample] += left * currentPan[0];
				aBuffer[aBufferSize + sample] += right * currentPan[1];
				aBuffer[aBufferSize * 2 + sample] += SURROUND_MIX_SCALE * (left + right) * currentPan[2];
				aBuffer[aBufferSize * 3 + sample] += SURROUND_MIX_SCALE * (left + right) * currentPan[3];
				aBuffer[aBufferSize * 4 + sample] += left * currentPan[4];
				aBuffer[aBufferSize * 5 + sample] += right * currentPan[5];
				aBuffer[aBufferSize * 6 + sample] += left * currentPan[6];
				aBuffer[aBufferSize * 7 + sample] += right * currentPan[7];
			}
			break;

		case 1: // 1.0 -> 7.1: distribute mono to all channels
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				for (int ch = 0; ch < 8; ch++)
				{
					currentPan[ch] += panDelta[ch];
				}

				float monoSrc = aScratch[sample];

				for (int ch = 0; ch < 8; ch++)
				{
					aBuffer[aBufferSize * ch + sample] += monoSrc * currentPan[ch];
				}
			}
			break;
		}
		break;
	}

	// Update voice state with final channel volumes
	for (unsigned int ch = 0; ch < aChannels; ch++)
	{
		aVoice->mCurrentChannelVolume[ch] = targetPan[ch];
	}
}

void resample_channels(float **srcChannels, float *outputBuffer, unsigned int outputStride, unsigned int numChannels, unsigned int outputSamples, double srcPosition,
                       double stepSize, unsigned int availableInput, unsigned int lookahead)
{
	using namespace ResamplingConstants;

	SOLOUD_ASSERT(srcChannels != nullptr);
	SOLOUD_ASSERT(outputBuffer != nullptr);
	SOLOUD_ASSERT(numChannels > 0);
	SOLOUD_ASSERT(availableInput >= lookahead);

#if defined(SOLOUD_AVX_INTRINSICS)
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
#endif

	for (unsigned int ch = 0; ch < numChannels; ch++)
	{
		float *src = srcChannels[ch];
		float *dst = outputBuffer + ch * outputStride;

		unsigned int i = 0;

#if defined(SOLOUD_AVX_INTRINSICS)
		// AVX2 path: process 8 samples at a time using gather instructions
		{
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
		}

#elif defined(SOLOUD_SSE_INTRINSICS)
		// SSE path: process 4 samples at a time with manual gathering
		{
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
		}
#endif

		// Scalar path: handles remainder after SIMD, or all samples if no SIMD available
		double srcPos = srcPosition + i * stepSize;
		for (; i < outputSamples; i++)
		{
			unsigned int intPos = (unsigned int)srcPos;
			float frac = (float)(srcPos - intPos);

			float sample;
			if (lookahead == POINT_LOOKAHEAD)
			{
				sample = src[intPos];
			}
			else if (lookahead == LINEAR_LOOKAHEAD)
			{
				float s0 = src[intPos];
				float s1 = (intPos + 1 < availableInput) ? src[intPos + 1] : s0;
				sample = s0 + (s1 - s0) * frac;
			}
			else // CATMULLROM_LOOKAHEAD
			{
				float p0 = (intPos >= 1) ? src[intPos - 1] : src[intPos];
				float p1 = src[intPos];
				float p2 = (intPos + 1 < availableInput) ? src[intPos + 1] : p1;
				float p3 = (intPos + 2 < availableInput) ? src[intPos + 2] : p1;

				float t = frac;
				float t2 = t * t;
				float t3 = t2 * t;

				sample =
				    CATMULLROM_SCALE *
				    ((CATMULLROM_LINEAR_COEFF * p1) + (-p0 + p2) * t +
				     (CATMULLROM_QUAD_COEFFS[0] * p0 + CATMULLROM_QUAD_COEFFS[1] * p1 + CATMULLROM_QUAD_COEFFS[2] * p2 + CATMULLROM_QUAD_COEFFS[3] * p3) * t2 +
				     (CATMULLROM_CUBIC_COEFFS[0] * p0 + CATMULLROM_CUBIC_COEFFS[1] * p1 + CATMULLROM_CUBIC_COEFFS[2] * p2 + CATMULLROM_CUBIC_COEFFS[3] * p3) * t3);
			}

			dst[i] = sample;
			srcPos += stepSize;
		}
	}

	return;
}

} // namespace SoLoud
