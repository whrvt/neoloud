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

#include "soloud_cpu.h"

#include <cfloat>
#include <cstdint>
#include <mutex>

#if defined(__GNUC__) && (defined(__i386__) || defined(__x86_64__))
#include <cpuid.h>
#elif defined(_M_IX86) || defined(_M_X64)
#include <intrin.h>
#endif

// copied from the cpuid header
#define CPUBIT_SSE (1 << 25)
#define CPUBIT_SSE2 (1 << 26)
#define CPUBIT_OSXSAVE (1 << 27)
#define CPUBIT_AVX (1 << 28)
#define CPUBIT_AVX2 (1 << 5)
#define CPUBIT_FMA (1 << 12)

namespace SoLoud
{

namespace
{
std::once_flag once_detect;
#if !defined(SOLOUD_DISABLE_SIMD)
size_t REAL_CPU_ALIGNMENT_BYTES{SCALAR_ALIGNMENT_BYTES};
size_t REAL_CPU_ALIGNMENT_MASK{SCALAR_ALIGNMENT_MASK};
size_t REAL_CPU_MEMORY_ALIGNMENT_MASK{SCALAR_MEMORY_ALIGNMENT_MASK};
#endif
} // namespace

#if !defined(SOLOUD_DISABLE_SIMD)
// Don't do any actual work here and rely on initCPUFeatures being called beforehand, for speed.
size_t CPU_ALIGNMENT_BYTES()
{
	return REAL_CPU_ALIGNMENT_BYTES;
}
size_t CPU_ALIGNMENT_MASK()
{
	return REAL_CPU_ALIGNMENT_MASK;
}
size_t CPU_MEMORY_ALIGNMENT_MASK()
{
	return REAL_CPU_MEMORY_ALIGNMENT_MASK;
}
#endif

void initCPUFeatures()
{
	static bool once = false;
	if (!once)
	{
#if !defined(SOLOUD_DISABLE_SIMD)
		const unsigned int CPUType = detectCPUextensions();
		if (CPUType & CPUFEATURE_AVX2)
		{
			REAL_CPU_ALIGNMENT_BYTES = AVX_ALIGNMENT_BYTES;
			REAL_CPU_ALIGNMENT_MASK = AVX_ALIGNMENT_MASK;
			REAL_CPU_MEMORY_ALIGNMENT_MASK = AVX_MEMORY_ALIGNMENT_MASK;
		}
		else if (CPUType & CPUFEATURE_SSE2)
		{
			REAL_CPU_ALIGNMENT_BYTES = SSE_ALIGNMENT_BYTES;
			REAL_CPU_ALIGNMENT_MASK = SSE_ALIGNMENT_MASK;
			REAL_CPU_MEMORY_ALIGNMENT_MASK = SSE_MEMORY_ALIGNMENT_MASK;
		}
#endif
		once = true;
	}

	// They are already initialized to scalar values by default.
}

unsigned int detectCPUextensions()
{
	static unsigned int detectedExtensions{0};
	// Nothing we can do here if we're not running on x86, at the moment.
#if defined(SOLOUD_IS_X86)
	std::call_once(once_detect, [](void) -> void {
		unsigned int res = 0;

#if defined(__GNUC__)
		unsigned int eax, ebx, ecx, edx;

		if (!__get_cpuid(1, &eax, &ebx, &ecx, &edx))
			return; // no CPUID support

		if (edx & CPUBIT_SSE)
			res |= CPUFEATURE_SSE;
		if (edx & CPUBIT_SSE2)
			res |= CPUFEATURE_SSE2;

		if ((ecx & CPUBIT_OSXSAVE) && (ecx & CPUBIT_AVX))
		{
			unsigned int xcr0_lo, xcr0_hi;
			__asm__ __volatile__("xgetbv" : "=a"(xcr0_lo), "=d"(xcr0_hi) : "c"(0));
			if ((xcr0_lo & 0x6) == 0x6)
			{
				res |= CPUFEATURE_AVX;
				if (ecx & CPUBIT_FMA)
					res |= CPUFEATURE_FMA;

				// check for AVX2 support (CPUID function 7, subfunction 0)
				unsigned int eax7, ebx7, ecx7, edx7;
				if (__get_cpuid_count(7, 0, &eax7, &ebx7, &ecx7, &edx7))
				{
					if (ebx7 & CPUBIT_AVX2)
						res |= CPUFEATURE_AVX2;
				}
			}
		}
#else // must be _MSC_VER
	int reg[4]{-1, 0, 0, 0};

	__cpuid(reg, 0);
	if ((unsigned int)reg[0] == 0)
		return; // no CPUID support

	__cpuid(reg, 1);
	if ((unsigned int)reg[3] & CPUBIT_SSE)
		res |= CPUFEATURE_SSE;
	if ((unsigned int)reg[3] & CPUBIT_SSE2)
		res |= CPUFEATURE_SSE2;

	if (((unsigned int)reg[2] & CPUBIT_OSXSAVE) && ((unsigned int)reg[2] & CPUBIT_AVX))
	{
		unsigned long long xcr0 = _xgetbv(0);
		if ((xcr0 & 0x6) == 0x6)
		{
			res |= CPUFEATURE_AVX;
			if ((unsigned int)reg[2] & CPUBIT_FMA)
				res |= CPUFEATURE_FMA;

			// check for AVX2 support (CPUID function 7, subfunction 0)
			__cpuidex(reg, 7, 0);
			if ((unsigned int)reg[1] & CPUBIT_AVX2)
				res |= CPUFEATURE_AVX2;
		}
	}

#endif // defined(__GNUC__)
		detectedExtensions = res;
	});

#endif // defined(SOLOUD_IS_X86)
	return detectedExtensions;
}

#if defined(SOLOUD_SUPPORT_SSE2) || defined(SOLOUD_IS_X86_64)
// src/intrin/sse/soloud_misc_sse.cpp
extern void setCTZDAZ();
#endif

void setFPUOptimizedRegs()
{
// Currently, MSVC builds will not do anything for the inline assembly parts (because it's not supported on MSVC)
#if !(defined(_MSC_VER) && !defined(__clang__))

// flush-to-zero (FTZ) for arm32/arm64
#if defined(__arm__) || defined(_ARM_)
	{
		__asm__ __volatile__("vmsr fpscr,%0" ::"r"(1 << 24));
	}
#endif

#if defined(_ARM64_) || defined(__aarch64__) || defined(__arm64__)
	{
		uint64_t fpcr;
		__asm__ __volatile__("mrs %0, fpcr" : "=r"(fpcr));
		fpcr |= (1ULL << 24); // FZ
		__asm__ __volatile__("msr fpcr, %0" : : "r"(fpcr));
	}
#endif

#endif // !(defined(_MSC_VER) && !defined(__clang__))

#ifdef _MCW_DN
	{
		_controlfp(_DN_FLUSH, _MCW_DN);
	}
#endif

#if defined(SOLOUD_SUPPORT_SSE2) || defined(SOLOUD_IS_X86_64)
	// For SSE (external compilation unit)
	if (detectCPUextensions() & CPUFEATURE_SSE)
	{
		setCTZDAZ();
	}
#endif
}

} // namespace SoLoud
