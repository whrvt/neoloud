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

// SIMD processing parameters
constexpr unsigned int SAMPLES_PER_SIMD_REGISTER = 4;
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
#ifndef SOLOUD_SSE_INTRINSICS
	mBasePtr = new unsigned char[aFloats * sizeof(float)];
	if (mBasePtr == NULL)
		return OUT_OF_MEMORY;
	mData = (float *)mBasePtr;
#else
	mBasePtr = new unsigned char[aFloats * sizeof(float) + 16];
	if (mBasePtr == nullptr)
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

TinyAlignedFloatBuffer::TinyAlignedFloatBuffer() : mActualData()
{
	unsigned char *basePtr = &mActualData[0];
	mData = (float *)(((size_t)basePtr + 15) & ~15);
}

#if defined(SOLOUD_SSE_INTRINSICS)
void clip_internal(const Soloud *aSoloud, AlignedFloatBuffer &aBuffer, AlignedFloatBuffer &aDestBuffer, unsigned int aSamples, float aVolume0, float aVolume1)
{
	using namespace ClippingConstants;

	float volumeDelta = (aVolume1 - aVolume0) / aSamples;
	float volume = aVolume0;
	unsigned int sampleQuads = (aSamples + 3) / 4; // Round up to process in SIMD groups

	if (aSoloud->mFlags & Soloud::CLIP_ROUNDOFF)
	{
		// Roundoff clipping: smooth saturation curve instead of hard clipping
		__m128 negThreshold = _mm_load_ps1(&ROUNDOFF_NEG_THRESHOLD);
		__m128 posThreshold = _mm_load_ps1(&ROUNDOFF_POS_THRESHOLD);
		__m128 linearScale = _mm_load_ps1(&ROUNDOFF_LINEAR_SCALE);
		__m128 cubicScale = _mm_load_ps1(&ROUNDOFF_CUBIC_SCALE);
		__m128 negWall = _mm_load_ps1(&ROUNDOFF_NEG_WALL);
		__m128 posWall = _mm_load_ps1(&ROUNDOFF_POS_WALL);
		__m128 postScale = _mm_load_ps1(&aSoloud->mPostClipScaler);

		// Prepare volume ramp for SIMD processing
		TinyAlignedFloatBuffer volumes;
		volumes.mData[0] = volume;
		volumes.mData[1] = volume + volumeDelta;
		volumes.mData[2] = volume + volumeDelta * 2;
		volumes.mData[3] = volume + volumeDelta * 3;
		volumeDelta *= SAMPLES_PER_SIMD_REGISTER;
		__m128 volDelta = _mm_load_ps1(&volumeDelta);

		unsigned int srcIdx = 0, dstIdx = 0;

		for (unsigned int channel = 0; channel < aSoloud->mChannels; channel++)
		{
			__m128 vol = _mm_load_ps(volumes.mData);

			for (unsigned int quad = 0; quad < sampleQuads; quad++)
			{
				// Load 4 samples and apply volume ramp
				__m128 samples = _mm_load_ps(&aBuffer.mData[srcIdx]);
				srcIdx += SAMPLES_PER_SIMD_REGISTER;
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
				_mm_store_ps(&aDestBuffer.mData[dstIdx], samples);
				dstIdx += SAMPLES_PER_SIMD_REGISTER;
			}
		}
	}
	else
	{
		// Hard clipping: simple min/max bounds
		__m128 negBound = _mm_load_ps1(&HARD_CLIP_MIN);
		__m128 posBound = _mm_load_ps1(&HARD_CLIP_MAX);
		__m128 postScale = _mm_load_ps1(&aSoloud->mPostClipScaler);

		// Prepare volume ramp for SIMD processing
		TinyAlignedFloatBuffer volumes;
		volumes.mData[0] = volume;
		volumes.mData[1] = volume + volumeDelta;
		volumes.mData[2] = volume + volumeDelta * 2;
		volumes.mData[3] = volume + volumeDelta * 3;
		volumeDelta *= SAMPLES_PER_SIMD_REGISTER;
		__m128 volDelta = _mm_load_ps1(&volumeDelta);

		unsigned int srcIdx = 0, dstIdx = 0;

		for (unsigned int channel = 0; channel < aSoloud->mChannels; channel++)
		{
			__m128 vol = _mm_load_ps(volumes.mData);

			for (unsigned int quad = 0; quad < sampleQuads; quad++)
			{
				// Load 4 samples and apply volume ramp
				__m128 samples = _mm_load_ps(&aBuffer.mData[srcIdx]);
				srcIdx += SAMPLES_PER_SIMD_REGISTER;
				samples = _mm_mul_ps(samples, vol);
				vol = _mm_add_ps(vol, volDelta);

				// Hard clipping to [-1, 1] range
				samples = _mm_max_ps(samples, negBound);
				samples = _mm_min_ps(samples, posBound);

				// Apply post-clip scaling and store
				samples = _mm_mul_ps(samples, postScale);
				_mm_store_ps(&aDestBuffer.mData[dstIdx], samples);
				dstIdx += SAMPLES_PER_SIMD_REGISTER;
			}
		}
	}
}

#else // Fallback implementation without SSE

void clip_internal(const Soloud *aSoloud, AlignedFloatBuffer &aBuffer, AlignedFloatBuffer &aDestBuffer, unsigned int aSamples, float aVolume0, float aVolume1)
{
	using namespace ClippingConstants;

	float volumeDelta = (aVolume1 - aVolume0) / aSamples;
	float volume = aVolume0;
	unsigned int sampleQuads = (aSamples + 3) / 4; // Process in groups of 4 for consistency

	if (aSoloud->mFlags & Soloud::CLIP_ROUNDOFF)
	{
		unsigned int srcIdx = 0, dstIdx = 0;

		for (unsigned int channel = 0; channel < aSoloud->mChannels; channel++)
		{
			volume = aVolume0;

			for (unsigned int quad = 0; quad < sampleQuads; quad++)
			{
				// Process 4 samples at a time
				for (int sample = 0; sample < SAMPLES_PER_SIMD_REGISTER; sample++)
				{
					float input = aBuffer.mData[srcIdx++] * volume;
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

					aDestBuffer.mData[dstIdx++] = output * aSoloud->mPostClipScaler;
				}
			}
		}
	}
	else
	{
		unsigned int srcIdx = 0, dstIdx = 0;

		for (unsigned int channel = 0; channel < aSoloud->mChannels; channel++)
		{
			volume = aVolume0;

			for (unsigned int quad = 0; quad < sampleQuads; quad++)
			{
				// Process 4 samples at a time
				for (int sample = 0; sample < SAMPLES_PER_SIMD_REGISTER; sample++)
				{
					float input = aBuffer.mData[srcIdx++] * volume;
					volume += volumeDelta;

					// Hard clipping to [-1, 1] range
					float output = (input <= HARD_CLIP_MIN) ? HARD_CLIP_MIN : (input >= HARD_CLIP_MAX) ? HARD_CLIP_MAX : input;

					aDestBuffer.mData[dstIdx++] = output * aSoloud->mPostClipScaler;
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
#ifdef SOLOUD_SSE_INTRINSICS
	SOLOUD_ASSERT(((size_t)aBuffer & 0xf) == 0);
	SOLOUD_ASSERT(((size_t)aScratch & 0xf) == 0);
	SOLOUD_ASSERT(((size_t)aBufferSize & 0xf) == 0);
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
#if defined(SOLOUD_SSE_INTRINSICS)
		{
			unsigned int processedSamples = 0;
			unsigned int sampleQuads = aSamplesToRead / SAMPLES_PER_SIMD_REGISTER;

			if (sampleQuads > 0)
			{
				// Prepare volume ramps for SIMD processing
				TinyAlignedFloatBuffer pan0, pan1;
				for (int i = 0; i < SAMPLES_PER_SIMD_REGISTER; i++)
				{
					pan0.mData[i] = currentPan[0] + panDelta[0] * (i + 1);
					pan1.mData[i] = currentPan[1] + panDelta[1] * (i + 1);
				}

				float panDelta0_4x = panDelta[0] * SAMPLES_PER_SIMD_REGISTER;
				float panDelta1_4x = panDelta[1] * SAMPLES_PER_SIMD_REGISTER;
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

					processedSamples += SAMPLES_PER_SIMD_REGISTER;
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
#if defined(SOLOUD_SSE_INTRINSICS)
		{
			unsigned int processedSamples = 0;
			unsigned int sampleQuads = aSamplesToRead / SAMPLES_PER_SIMD_REGISTER;

			if (sampleQuads > 0)
			{
				// Prepare volume ramps for SIMD processing
				TinyAlignedFloatBuffer pan0, pan1;
				for (int i = 0; i < SAMPLES_PER_SIMD_REGISTER; i++)
				{
					pan0.mData[i] = currentPan[0] + panDelta[0] * (i + 1);
					pan1.mData[i] = currentPan[1] + panDelta[1] * (i + 1);
				}

				float panDelta0_4x = panDelta[0] * SAMPLES_PER_SIMD_REGISTER;
				float panDelta1_4x = panDelta[1] * SAMPLES_PER_SIMD_REGISTER;
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

					processedSamples += SAMPLES_PER_SIMD_REGISTER;
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

} // namespace SoLoud
