/*
SoLoud audio engine - CPU detection/SIMD intrinsics utilities
Copyright (c) 2013-2020 Jari Komppa
Copyright (c) 2026 William Horvath

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

#ifndef SOLOUD_CPUDETECT_H
#define SOLOUD_CPUDETECT_H

#include <cstddef>
#include <cstdint>

#if defined(__i386__) || defined(__x86_64__) || defined(_M_IX86) || defined(_M_X64)
#define SOLOUD_IS_X86
#if defined(__x86_64__) || defined(_M_X64)
#define SOLOUD_IS_X86_64
#endif
#endif

namespace SoLoud
{
enum CPUFeature : unsigned int
{
	CPUFEATURE_SSE = (1 << 0),
	CPUFEATURE_SSE2 = (1 << 1),
	CPUFEATURE_AVX = (1 << 2),
	CPUFEATURE_AVX2 = (1 << 3),
	CPUFEATURE_FMA = (1 << 4),
};

// SIMD configuration constants

// AVX2 configuration
inline constexpr size_t AVX_ALIGNMENT_BYTES = 32;      // 32-byte alignment for AVX2
inline constexpr size_t AVX_ALIGNMENT_MASK = 31;       // Mask for 32-byte alignment check
inline constexpr size_t AVX_MEMORY_ALIGNMENT_MASK = 7; // Align to 8-float boundaries

// SSE configuration
inline constexpr size_t SSE_ALIGNMENT_BYTES = 16;      // 16-byte alignment for SSE
inline constexpr size_t SSE_ALIGNMENT_MASK = 15;       // Mask for 16-byte alignment check (0xf)
inline constexpr size_t SSE_MEMORY_ALIGNMENT_MASK = 3; // Align to 4-float boundaries

// Scalar fallback configuration
inline constexpr size_t SCALAR_ALIGNMENT_BYTES = sizeof(void *);
inline constexpr size_t SCALAR_ALIGNMENT_MASK = sizeof(void *) - 1;
inline constexpr size_t SCALAR_MEMORY_ALIGNMENT_MASK = 3;

// Alignment requirements for buffers.
// Initialized after calling initCPUFeatures(), which is also done at the start of the Soloud constructor.
#if defined(SOLOUD_DISABLE_SIMD)
inline constexpr size_t CPU_ALIGNMENT_BYTES()
{
	return sizeof(void *);
}
inline constexpr size_t CPU_ALIGNMENT_MASK()
{
	return sizeof(void *) - 1;
}
inline constexpr size_t CPU_MEMORY_ALIGNMENT_MASK()
{
	return 3;
} // 4-float boundaries
#else

#ifdef __GNUC__
#define SL_PURECALL __attribute__((pure))
#else
#define SL_PURECALL
#endif

size_t CPU_ALIGNMENT_BYTES() SL_PURECALL;
size_t CPU_ALIGNMENT_MASK() SL_PURECALL;
size_t CPU_MEMORY_ALIGNMENT_MASK() SL_PURECALL;
#endif

/**
 * Initializes alignment constants, depending on the detected runtime architecture.
 *
 */
void initCPUFeatures();

/**
 * Checks which instruction set extensions are supported by the CPU.
 *
 * @return                  A bitmask of supported CPU extensions (CPUFeature enum).
 */
unsigned int detectCPUextensions();

/**
 * General-purpose function that sets floating point state for the currently detected CPU.
 *
 * For x86 (all): sets flush denormals (_DN_FLUSH)
 * For SSE: sets denorm clear to zero (CTZ) and denorms are zero (DAZ)
 * For aarch64 and arm32: sets flush-to-zero
 *
 * Note: this state is thread-local, and only needs to be set once.
 */
void setFPUOptimizedRegs();

} // namespace SoLoud

#endif
