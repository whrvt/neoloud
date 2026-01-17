
/*
SoLoud audio engine - miscellaneous SIMD-related things
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

#include "soloud.h"
#include "soloud_cpu.h"
#include <cstring>

namespace SoLoud
{

// This is publicly-accessible, declared in include/soloud.h, but uses the correct alignment depending on the CPU features.
// Does not rely on the global constants being defined, since this is a standalone class.
result AlignedFloatBuffer::init(unsigned int aFloats)
{
	delete[] mBasePtr;
	mBasePtr = nullptr;
	mData = nullptr;
	mFloats = aFloats;

	size_t alignmentBytes, alignmentMask;
#if !defined(DISABLE_SIMD)
	const unsigned int CPUType = detectCPUextensions();
	if (CPUType & CPUFEATURE_AVX2)
	{
		alignmentBytes = AVX_ALIGNMENT_BYTES;
		alignmentMask = AVX_ALIGNMENT_MASK;
	}
	else if (CPUType & CPUFEATURE_SSE2)
	{
		alignmentBytes = SSE_ALIGNMENT_BYTES;
		alignmentMask = SSE_ALIGNMENT_MASK;
	}
	else
#endif
	{
		alignmentBytes = SCALAR_ALIGNMENT_BYTES;
		alignmentMask = SCALAR_ALIGNMENT_MASK;
	}

	mBasePtr = new unsigned char[aFloats * sizeof(float) + alignmentBytes];
	if (mBasePtr == nullptr)
		return OUT_OF_MEMORY;
	mData = (float *)(((size_t)mBasePtr + alignmentMask) & ~alignmentMask);

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

} // namespace SoLoud