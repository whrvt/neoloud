/*
SoLoud audio engine - LUFS K-weighting (SSE2-optimized)
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

#if defined(SOLOUD_SUPPORT_SSE2)
#include "soloud_config.h"
#include "soloud_loudness_internal.h"

#ifdef _MSC_VER
#include <intrin.h>
#endif
#include <emmintrin.h>
#include <xmmintrin.h>

#include <string.h>

// NOLINTBEGIN(cppcoreguidelines-avoid-c-arrays, modernize-avoid-c-arrays)

namespace SoLoud::Loudness
{
namespace
{
// Transpose strided channel data to interleaved layout: [ch * stride + i] -> [i * 4 + ch]
// Zeroes unused lanes (aChannels < 4).
inline void transposeToInterleaved(const float *aStrided, unsigned int aStride, unsigned int aChannels, unsigned int aSamples, float *aInterleaved)
{
	if (aChannels == 4)
	{
		for (unsigned int i = 0; i < aSamples; i++)
		{
			float *dst = aInterleaved + i * 4;
			dst[0] = aStrided[i];
			dst[1] = aStrided[aStride + i];
			dst[2] = aStrided[2 * aStride + i];
			dst[3] = aStrided[3 * aStride + i];
		}
	}
	else if (aChannels == 2)
	{
		for (unsigned int i = 0; i < aSamples; i++)
		{
			float *dst = aInterleaved + i * 4;
			dst[0] = aStrided[i];
			dst[1] = aStrided[aStride + i];
			dst[2] = 0.0f;
			dst[3] = 0.0f;
		}
	}
	else
	{
		for (unsigned int i = 0; i < aSamples; i++)
		{
			float *dst = aInterleaved + i * 4;
			for (unsigned int ch = 0; ch < aChannels; ch++)
				dst[ch] = aStrided[ch * aStride + i];
			for (unsigned int ch = aChannels; ch < 4; ch++)
				dst[ch] = 0.0f;
		}
	}
}

// Transpose interleaved layout back to strided: [i * 4 + ch] -> [ch * stride + i]
inline void transposeToStrided(const float *aInterleaved, unsigned int aStride, unsigned int aChannels, unsigned int aSamples, float *aStrided)
{
	if (aChannels == 2)
	{
		for (unsigned int i = 0; i < aSamples; i++)
		{
			const float *src = aInterleaved + i * 4;
			aStrided[i] = src[0];
			aStrided[aStride + i] = src[1];
		}
	}
	else if (aChannels == 4)
	{
		for (unsigned int i = 0; i < aSamples; i++)
		{
			const float *src = aInterleaved + i * 4;
			aStrided[i] = src[0];
			aStrided[aStride + i] = src[1];
			aStrided[2 * aStride + i] = src[2];
			aStrided[3 * aStride + i] = src[3];
		}
	}
	else
	{
		for (unsigned int i = 0; i < aSamples; i++)
		{
			const float *src = aInterleaved + i * 4;
			for (unsigned int ch = 0; ch < aChannels; ch++)
				aStrided[ch * aStride + i] = src[ch];
		}
	}
}

// Apply two cascaded biquad filters on interleaved [sample * 4 + ch] data using SSE.
// aData is modified in-place. aSamples is the number of samples to process.
inline void biquadInterleavedSSE(const BiquadCoeffs &aPreCoeffs, const BiquadCoeffs &aRlbCoeffs, __m128 &aPreX1, __m128 &aPreX2, __m128 &aPreY1, __m128 &aPreY2,
                                 __m128 &aRlbX1, __m128 &aRlbX2, __m128 &aRlbY1, __m128 &aRlbY2, float *aData, unsigned int aSamples)
{
	__m128 preB0 = _mm_set1_ps(aPreCoeffs.b0);
	__m128 preB1 = _mm_set1_ps(aPreCoeffs.b1);
	__m128 preB2 = _mm_set1_ps(aPreCoeffs.b2);
	__m128 preA1 = _mm_set1_ps(aPreCoeffs.a1);
	__m128 preA2 = _mm_set1_ps(aPreCoeffs.a2);

	__m128 rlbB0 = _mm_set1_ps(aRlbCoeffs.b0);
	__m128 rlbB1 = _mm_set1_ps(aRlbCoeffs.b1);
	__m128 rlbB2 = _mm_set1_ps(aRlbCoeffs.b2);
	__m128 rlbA1 = _mm_set1_ps(aRlbCoeffs.a1);
	__m128 rlbA2 = _mm_set1_ps(aRlbCoeffs.a2);

	for (unsigned int i = 0; i < aSamples; i++)
	{
		__m128 x = _mm_load_ps(aData + i * 4);

		// pre-filter biquad: y = b0*x + b1*x1 + b2*x2 - a1*y1 - a2*y2
		__m128 y = _mm_add_ps(_mm_add_ps(_mm_mul_ps(preB0, x), _mm_mul_ps(preB1, aPreX1)), _mm_mul_ps(preB2, aPreX2));
		y = _mm_sub_ps(y, _mm_add_ps(_mm_mul_ps(preA1, aPreY1), _mm_mul_ps(preA2, aPreY2)));
		aPreX2 = aPreX1;
		aPreX1 = x;
		aPreY2 = aPreY1;
		aPreY1 = y;

		// RLB biquad (input is pre-filter output)
		x = y;
		y = _mm_add_ps(_mm_add_ps(_mm_mul_ps(rlbB0, x), _mm_mul_ps(rlbB1, aRlbX1)), _mm_mul_ps(rlbB2, aRlbX2));
		y = _mm_sub_ps(y, _mm_add_ps(_mm_mul_ps(rlbA1, aRlbY1), _mm_mul_ps(rlbA2, aRlbY2)));
		aRlbX2 = aRlbX1;
		aRlbX1 = x;
		aRlbY2 = aRlbY1;
		aRlbY1 = y;

		_mm_store_ps(aData + i * 4, y);
	}
}

// Helper: pack biquad states [0..aChannels-1] into __m128 (zero-padding unused lanes)
inline __m128 packState(const BiquadState *aStates, unsigned int aChannels, float BiquadState::*aField)
{
	alignas(16) float tmp[4] = {0.0f, 0.0f, 0.0f, 0.0f};
	for (unsigned int ch = 0; ch < aChannels; ch++)
		tmp[ch] = aStates[ch].*aField;
	return _mm_load_ps(tmp);
}

// Helper: unpack __m128 back into biquad states [0..aChannels-1]
inline void unpackState(BiquadState *aStates, unsigned int aChannels, float BiquadState::*aField, __m128 aVal)
{
	alignas(16) float tmp[4];
	_mm_store_ps(tmp, aVal);
	for (unsigned int ch = 0; ch < aChannels; ch++)
		aStates[ch].*aField = tmp[ch];
}
} // namespace

void kWeightChunkSSE(const BiquadCoeffs &aPreCoeffs, const BiquadCoeffs &aRlbCoeffs, BiquadState *aPreState, BiquadState *aRlbState, const float *aBuffer,
                     unsigned int aSamplesRead, unsigned int aBufStride, unsigned int aChannels, float *aOut)
{
	// mono: no benefit from SSE (we should not have been called in the first place)
	SOLOUD_ASSERT(aChannels > 1);

	// interleaved temp buffer: [sample * 4 + ch], 16-byte aligned for SSE loads/stores
	// SAMPLE_GRANULARITY * 4 floats = 8KB, fits in L1 cache
	alignas(16) float interleaved[SAMPLE_GRANULARITY * 4];

	if (aChannels <= 4)
	{
		// transpose strided -> interleaved
		transposeToInterleaved(aBuffer, aBufStride, aChannels, aSamplesRead, interleaved);

		// pack biquad states into SSE registers
		__m128 preX1 = packState(aPreState, aChannels, &BiquadState::x1);
		__m128 preX2 = packState(aPreState, aChannels, &BiquadState::x2);
		__m128 preY1 = packState(aPreState, aChannels, &BiquadState::y1);
		__m128 preY2 = packState(aPreState, aChannels, &BiquadState::y2);
		__m128 rlbX1 = packState(aRlbState, aChannels, &BiquadState::x1);
		__m128 rlbX2 = packState(aRlbState, aChannels, &BiquadState::x2);
		__m128 rlbY1 = packState(aRlbState, aChannels, &BiquadState::y1);
		__m128 rlbY2 = packState(aRlbState, aChannels, &BiquadState::y2);

		// SSE biquad on contiguous interleaved data
		biquadInterleavedSSE(aPreCoeffs, aRlbCoeffs, preX1, preX2, preY1, preY2, rlbX1, rlbX2, rlbY1, rlbY2, interleaved, aSamplesRead);

		// unpack biquad states
		unpackState(aPreState, aChannels, &BiquadState::x1, preX1);
		unpackState(aPreState, aChannels, &BiquadState::x2, preX2);
		unpackState(aPreState, aChannels, &BiquadState::y1, preY1);
		unpackState(aPreState, aChannels, &BiquadState::y2, preY2);
		unpackState(aRlbState, aChannels, &BiquadState::x1, rlbX1);
		unpackState(aRlbState, aChannels, &BiquadState::x2, rlbX2);
		unpackState(aRlbState, aChannels, &BiquadState::y1, rlbY1);
		unpackState(aRlbState, aChannels, &BiquadState::y2, rlbY2);

		// transpose interleaved -> strided output
		transposeToStrided(interleaved, aBufStride, aChannels, aSamplesRead, aOut);
	}
	else
	{
		// 5-8 channels: process as two groups of 4

		// --- lo group: channels 0..3 ---
		transposeToInterleaved(aBuffer, aBufStride, 4, aSamplesRead, interleaved);

		__m128 preX1 = packState(aPreState, 4, &BiquadState::x1);
		__m128 preX2 = packState(aPreState, 4, &BiquadState::x2);
		__m128 preY1 = packState(aPreState, 4, &BiquadState::y1);
		__m128 preY2 = packState(aPreState, 4, &BiquadState::y2);
		__m128 rlbX1 = packState(aRlbState, 4, &BiquadState::x1);
		__m128 rlbX2 = packState(aRlbState, 4, &BiquadState::x2);
		__m128 rlbY1 = packState(aRlbState, 4, &BiquadState::y1);
		__m128 rlbY2 = packState(aRlbState, 4, &BiquadState::y2);

		biquadInterleavedSSE(aPreCoeffs, aRlbCoeffs, preX1, preX2, preY1, preY2, rlbX1, rlbX2, rlbY1, rlbY2, interleaved, aSamplesRead);

		unpackState(aPreState, 4, &BiquadState::x1, preX1);
		unpackState(aPreState, 4, &BiquadState::x2, preX2);
		unpackState(aPreState, 4, &BiquadState::y1, preY1);
		unpackState(aPreState, 4, &BiquadState::y2, preY2);
		unpackState(aRlbState, 4, &BiquadState::x1, rlbX1);
		unpackState(aRlbState, 4, &BiquadState::x2, rlbX2);
		unpackState(aRlbState, 4, &BiquadState::y1, rlbY1);
		unpackState(aRlbState, 4, &BiquadState::y2, rlbY2);

		transposeToStrided(interleaved, aBufStride, 4, aSamplesRead, aOut);

		// --- hi group: channels 4..N-1 ---
		unsigned int hiChannels = aChannels - 4;

		transposeToInterleaved(aBuffer + 4 * aBufStride, aBufStride, hiChannels, aSamplesRead, interleaved);

		preX1 = packState(aPreState + 4, hiChannels, &BiquadState::x1);
		preX2 = packState(aPreState + 4, hiChannels, &BiquadState::x2);
		preY1 = packState(aPreState + 4, hiChannels, &BiquadState::y1);
		preY2 = packState(aPreState + 4, hiChannels, &BiquadState::y2);
		rlbX1 = packState(aRlbState + 4, hiChannels, &BiquadState::x1);
		rlbX2 = packState(aRlbState + 4, hiChannels, &BiquadState::x2);
		rlbY1 = packState(aRlbState + 4, hiChannels, &BiquadState::y1);
		rlbY2 = packState(aRlbState + 4, hiChannels, &BiquadState::y2);

		biquadInterleavedSSE(aPreCoeffs, aRlbCoeffs, preX1, preX2, preY1, preY2, rlbX1, rlbX2, rlbY1, rlbY2, interleaved, aSamplesRead);

		unpackState(aPreState + 4, hiChannels, &BiquadState::x1, preX1);
		unpackState(aPreState + 4, hiChannels, &BiquadState::x2, preX2);
		unpackState(aPreState + 4, hiChannels, &BiquadState::y1, preY1);
		unpackState(aPreState + 4, hiChannels, &BiquadState::y2, preY2);
		unpackState(aRlbState + 4, hiChannels, &BiquadState::x1, rlbX1);
		unpackState(aRlbState + 4, hiChannels, &BiquadState::x2, rlbX2);
		unpackState(aRlbState + 4, hiChannels, &BiquadState::y1, rlbY1);
		unpackState(aRlbState + 4, hiChannels, &BiquadState::y2, rlbY2);

		transposeToStrided(interleaved, aBufStride, hiChannels, aSamplesRead, aOut + 4 * aBufStride);
	}
}
} // namespace SoLoud::Loudness

// NOLINTEND(cppcoreguidelines-avoid-c-arrays, modernize-avoid-c-arrays)
#endif
