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

namespace SoLoud
{
	// Class that handles aligned allocations to support vectorized operations
	class AlignedFloatBuffer
	{
	public:
		float *mData;            // aligned pointer
		unsigned char *mBasePtr; // raw allocated pointer (for delete)
		int mFloats;             // size of buffer (w/out padding)

		// ctor
		AlignedFloatBuffer();
		// Allocate and align buffer
		result init(unsigned int aFloats);
		// Clear data to zero.
		void clear();
		// dtor
		~AlignedFloatBuffer();
	};

	// Lightweight class that handles small aligned buffer to support vectorized operations
	class TinyAlignedFloatBuffer
	{
	public:
		float *mData; // aligned pointer
		unsigned char mActualData[sizeof(float) * 16 + 16];

		// ctor
		TinyAlignedFloatBuffer();
	};

	class AudioSourceInstance;
	void clip_internal(const Soloud *aSoloud, AlignedFloatBuffer &aBuffer, AlignedFloatBuffer &aDestBuffer, unsigned int aSamples, float aVolume0, float aVolume1);
	void panAndExpand(AudioSourceInstance * aVoice, float *aBuffer, unsigned int aSamplesToRead, unsigned int aBufferSize, float *aScratch, unsigned int aChannels);
}

#endif
