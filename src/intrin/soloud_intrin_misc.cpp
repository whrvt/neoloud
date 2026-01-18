
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
#include <new>

namespace SoLoud
{

namespace
{

std::align_val_t alignment{0};

template <typename R>
R *allocateAligned(size_t sizeBytes)
{
	if (alignment == static_cast<std::align_val_t>(0))
		initCPUFeatures();

	alignment = static_cast<std::align_val_t>(CPU_ALIGNMENT_BYTES());
	return static_cast<R *>(::operator new(sizeBytes * sizeof(R), alignment));
}

void freeAligned(void *data)
{
	if (!data)
		return;

	SOLOUD_ASSERT(alignment != static_cast<std::align_val_t>(0));

	::operator delete(data, alignment);
}

} // namespace

// This is publicly-accessible, declared in include/soloud.h, but uses the correct alignment depending on the CPU features.
result AlignedFloatBuffer::init(unsigned int aFloats)
{
	if (mData)
		freeAligned(mData);
	mData = nullptr;
	mFloats = aFloats;

	if (!(mData = allocateAligned<float>(mFloats)))
		return OUT_OF_MEMORY;

	return SO_NO_ERROR;
}

void AlignedFloatBuffer::clear()
{
	if (mData) // memset(NULL,...) is UB
		memset(mData, 0, sizeof(float) * mFloats);
}

AlignedFloatBuffer::~AlignedFloatBuffer()
{
	freeAligned(mData);
}

} // namespace SoLoud