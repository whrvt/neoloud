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

#ifndef SOLOUD_INTRIN_H
#define SOLOUD_INTRIN_H

#include "soloud_config.h"

#ifdef SOLOUD_AVX_INTRINSICS
#include <immintrin.h>
#elif defined(SOLOUD_SSE_INTRINSICS)
#include <xmmintrin.h>
#ifdef _M_IX86
#include <emmintrin.h>
#endif
#else
#include <cstddef>
#endif

namespace SoLoud
{
	// SIMD configuration constants - adjust based on available instruction sets
#ifdef SOLOUD_AVX_INTRINSICS
	// AVX2 configuration
	static constexpr size_t SIMD_ALIGNMENT_BYTES = 32;    // 32-byte alignment for AVX2
	static constexpr size_t SIMD_ALIGNMENT_MASK = 31;     // Mask for 32-byte alignment check
	static constexpr size_t OPTIMAL_CHUNK_SAMPLES = 8;    // 8 floats per AVX2 register
	static constexpr size_t MEMORY_ALIGNMENT_MASK = 7;    // Align to 8-float boundaries
	static constexpr size_t TINY_BUFFER_FLOAT_COUNT = 32; // Buffer space for temporary operations
#elif defined(SOLOUD_SSE_INTRINSICS)
	// SSE configuration
	static constexpr size_t SIMD_ALIGNMENT_BYTES = 16;    // 16-byte alignment for SSE
	static constexpr size_t SIMD_ALIGNMENT_MASK = 15;     // Mask for 16-byte alignment check (0xf)
	static constexpr size_t OPTIMAL_CHUNK_SAMPLES = 4;    // 4 floats per SSE register
	static constexpr size_t MEMORY_ALIGNMENT_MASK = 3;    // Align to 4-float boundaries
	static constexpr size_t TINY_BUFFER_FLOAT_COUNT = 16; // Buffer space for temporary operations
#else
	// Scalar fallback configuration
	static constexpr size_t SIMD_ALIGNMENT_BYTES = 4;     // Basic float alignment
	static constexpr size_t SIMD_ALIGNMENT_MASK = 3;      // Mask for 4-byte alignment check
	static constexpr size_t OPTIMAL_CHUNK_SAMPLES = 4;    // Still process 4 samples at a time even for scalar
	static constexpr size_t MEMORY_ALIGNMENT_MASK = 3;    // Align to 4-byte boundaries
	static constexpr size_t TINY_BUFFER_FLOAT_COUNT = 16; // Buffer space for temporary operations (keeping this the same as SSE, this is how it was before)
#endif
	// Class that handles aligned allocations to support vectorized operations
	class AlignedFloatBuffer
	{
	public:
		float *mData;            // SIMD-aligned pointer for vectorized operations
		unsigned char *mBasePtr; // Raw allocated pointer (for delete)
		unsigned int mFloats;             // Size of buffer in floats (without padding)

		// Constructor
		AlignedFloatBuffer();

		// Allocate and align buffer for specified number of floats
		result init(unsigned int aFloats);

		// Clear all data to zero
		void clear();

		// Destructor
		~AlignedFloatBuffer();
	};

	// Lightweight class that handles small aligned buffer to support vectorized operations
	// Used for temporary SIMD register-sized data
	class TinyAlignedFloatBuffer
	{
	public:
		float *mData;                                                                              // SIMD-aligned pointer
		unsigned char mActualData[sizeof(float) * TINY_BUFFER_FLOAT_COUNT + SIMD_ALIGNMENT_BYTES]; // Space for appropriate number of floats + alignment padding

		// Constructor - automatically aligns mData to proper boundary
		TinyAlignedFloatBuffer();
	};

	class AudioSourceInstance;
	class Soloud;

	/* Soloud::clip_internal defined in soloud.h for convenience access to class members */

	/**
	 * Convert channel-separated audio data to interleaved format with sample format conversion
	 *
	 * @param outputBuffer   Destination buffer for interleaved samples
	 * @param rawBuffer      Source buffer with channel-separated float data
	 * @param aSamples       Number of samples per channel to process
	 * @param stride         Size of each channel buffer (including padding for alignment)
	 * @param aChannels      Number of audio channels
	 * @param aFormat        Target sample format for conversion
	 *
	 * Converts from SoLoud's internal channel-separated float format to interleaved output
	 * in the specified sample format. Handles scaling, clamping, and type conversion for:
	 * - SAMPLE_FLOAT32: Direct copy (32-bit float)
	 * - SAMPLE_UNSIGNED8: Scale to [0,255] range (8-bit unsigned)
	 * - SAMPLE_SIGNED16: Scale to [-32767,32767] range (16-bit signed)
	 * - SAMPLE_SIGNED24: Scale to 24-bit range, store as 3-byte packed
	 * - SAMPLE_SIGNED32: Scale to full 32-bit signed integer range
	 *
	 */
	void interlace_samples(void *outputBuffer, const float *const &rawBuffer, const unsigned int &aSamples, const unsigned int &stride,
	                       const unsigned int &aChannels, const detail::SAMPLE_FORMAT &aFormat = detail::SAMPLE_FLOAT32);

	/**
	 * Pan and expand audio from source channel count to output channel count
	 *
	 * @param aVoice         Voice instance containing channel volume settings
	 * @param aBuffer        Output buffer (accumulative mixing)
	 * @param aSamplesToRead Number of samples to process
	 * @param aBufferSize    Size of each channel buffer (for stride calculation)
	 * @param aScratch       Input source data (separated by channel)
	 * @param aChannels      Number of output channels
	 *
	 * Handles all common channel configurations:
	 * - 1.0 (mono) -> 1.0, 2.0, 4.0, 5.1, 7.1
	 * - 2.0 (stereo) -> 1.0, 2.0, 4.0, 5.1, 7.1
	 * - 4.0 (quad) -> 1.0, 2.0, 4.0, 5.1, 7.1
	 * - 5.1 (surround) -> 1.0, 2.0, 4.0, 5.1, 7.1
	 * - 7.1 (surround) -> 1.0, 2.0, 4.0, 5.1, 7.1
	 *
	 * Provides smooth volume ramping to avoid audio artifacts.
	 * Uses intelligent downmixing coefficients to preserve perceived loudness.
	 * Optimized with AVX2 intrinsics for stereo output paths when available,
	 * falls back to SSE intrinsics or scalar implementation.
	 */
	void panAndExpand(AudioSourceInstance * aVoice, float *aBuffer, unsigned int aSamplesToRead, unsigned int aBufferSize, float *aScratch, unsigned int aChannels);
}

#endif
