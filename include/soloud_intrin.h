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

#ifdef SOLOUD_SSE_INTRINSICS
#include <xmmintrin.h>

#include <cstddef>
#ifdef _M_IX86
#include <emmintrin.h>
#endif
#endif

namespace SoLoud
{
	// Class that handles aligned allocations to support vectorized operations
	class AlignedFloatBuffer
	{
	public:
		float *mData;            // 16-byte aligned pointer for SIMD operations
		unsigned char *mBasePtr; // Raw allocated pointer (for delete)
		int mFloats;             // Size of buffer in floats (without padding)

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
	// Used for temporary SIMD register-sized data (4 floats)
	class TinyAlignedFloatBuffer
	{
	public:
		float *mData;                                       // 16-byte aligned pointer
		unsigned char mActualData[sizeof(float) * 16 + 16]; // Space for 16 floats + alignment padding

		// Constructor - automatically aligns mData to 16-byte boundary
		TinyAlignedFloatBuffer();
	};

	class AudioSourceInstance;
	class Soloud;

	/**
	 * Apply volume scaling and clipping to audio buffer
	 *
	 * @param aSoloud       SoLoud instance (for configuration flags)
	 * @param aBuffer       Input buffer with source samples
	 * @param aDestBuffer   Output buffer for processed samples
	 * @param aSamples      Number of samples to process
	 * @param aVolume0      Starting volume level
	 * @param aVolume1      Ending volume level (for smooth ramping)
	 *
	 * Supports two clipping modes:
	 * - Hard clipping: Simple [-1,1] bounds
	 * - Roundoff clipping: Smooth saturation curve that approaches limits asymptotically
	 *
	 * Uses SSE intrinsics when available for 4x performance improvement
	 */
	void clip_internal(const Soloud *aSoloud, AlignedFloatBuffer &aBuffer, AlignedFloatBuffer &aDestBuffer, unsigned int aSamples, float aVolume0, float aVolume1);

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
	 * Optimized with SSE intrinsics for stereo output paths.
	 */
	void panAndExpand(AudioSourceInstance * aVoice, float *aBuffer, unsigned int aSamplesToRead, unsigned int aBufferSize, float *aScratch, unsigned int aChannels);
}

#endif
