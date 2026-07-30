/*
SoLoud audio engine - custom container types
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

#ifndef SOLOUD_CONTAINERS_H
#define SOLOUD_CONTAINERS_H

#include <stddef.h>
#include "soloud_config.h"

namespace SoLoud
{

// Class that handles aligned allocations to support vectorized operations
// Implementation in src/core/soloud_containers.cpp
class AlignedFloatBuffer
{
public:
	float *mData{nullptr};   // SIMD-aligned pointer for vectorized operations
	unsigned int mFloats{0}; // Size of buffer in floats (without padding)

	AlignedFloatBuffer() = default;
	AlignedFloatBuffer(unsigned int aFloats) { init(aFloats); }
	~AlignedFloatBuffer();

	// Not copy/moveable.
	AlignedFloatBuffer(const AlignedFloatBuffer &) = delete;
	AlignedFloatBuffer &operator=(const AlignedFloatBuffer &) = delete;
	AlignedFloatBuffer(AlignedFloatBuffer &&) = delete;
	AlignedFloatBuffer &operator=(AlignedFloatBuffer &&) = delete;

	// Allocate and align buffer for specified number of floats
	result init(unsigned int aFloats);

	// Clear all data to zero
	void clear();
};

// Array that doesn't allocate for small amounts of floats
// Mainly intended for use by AudioSourceInstance's source resampling buffers
// NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init) - mInline starts uninitialized, users zero out what they use
struct InlineFloatArray
{
	// inline 2 channels
	static constexpr const size_t INLINE_CHANNELS{2};
	static constexpr const size_t INLINE_FLOATS{((size_t)SAMPLE_GRANULARITY * 3) * INLINE_CHANNELS};

	InlineFloatArray() = default;
	~InlineFloatArray();

	// Not copy/moveable.
	InlineFloatArray(const InlineFloatArray &) = delete;
	InlineFloatArray &operator=(const InlineFloatArray &) = delete;
	InlineFloatArray(InlineFloatArray &&) = delete;
	InlineFloatArray &operator=(InlineFloatArray &&) = delete;

	[[nodiscard]] float *get() { return mHeap ? mHeap : mInline; }
	[[nodiscard]] const float *get() const { return mHeap ? mHeap : mInline; }
	[[nodiscard]] size_t capacity() const { return mCapacity; }

	// grow-only; newly acquired storage is uninitialized, existing contents are not preserved
	void ensureCapacity(size_t aFloats);

private:
	// use max (AVX) alignment for inline floats, only check heap allocation alignment dynamically
	alignas(32) float mInline[INLINE_FLOATS];
	float *mHeap{nullptr};
	size_t mCapacity{INLINE_FLOATS};
};

} // namespace SoLoud

#endif
