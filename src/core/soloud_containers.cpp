/*
SoLoud audio engine - helper containers
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

#include "soloud_config.h"
#include "soloud_cpu.h"
#include "soloud_containers.h"
#include "soloud_error.h"

#include <string.h>
#include <new>

namespace SoLoud
{

namespace
{

std::align_val_t getAlignment()
{
	// magic static, so concurrent first-time allocations from multiple threads are safe
	static const std::align_val_t alignment{[] {
		initCPUFeatures();
		return static_cast<std::align_val_t>(CPU_ALIGNMENT_BYTES());
	}()};
	return alignment;
}

template <typename R>
R *allocateAligned(size_t aCount)
{
	// with -fno-exceptions, allocation failure terminates instead of throwing, so this never returns null
	return static_cast<R *>(::operator new(aCount * sizeof(R), getAlignment()));
}

void freeAligned(void *aData)
{
	if (!aData)
		return;

	::operator delete(aData, getAlignment());
}

} // namespace

// This is publicly-accessible, declared in include/soloud_containers.h, but uses the correct alignment depending on the CPU features.
result AlignedFloatBuffer::init(unsigned int aFloats)
{
	freeAligned(mData);
	mFloats = aFloats;
	mData = allocateAligned<float>(mFloats);
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

InlineFloatArray::~InlineFloatArray()
{
	freeAligned(mHeap);
}

void InlineFloatArray::ensureCapacity(size_t aFloats)
{
	if (aFloats <= mCapacity)
		return;

	float *bigger = allocateAligned<float>(aFloats);
	freeAligned(mHeap);
	mHeap = bigger;
	mCapacity = aFloats;
}

} // namespace SoLoud
