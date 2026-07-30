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

#pragma once

// NOLINTBEGIN(cppcoreguidelines-avoid-c-arrays, modernize-avoid-c-arrays)

#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "soloud.h"
#include "soloud_wav.h"

// This option is useful while developing tests:
// #define NO_LASTKNOWN_CHECK

extern int gErrorCount;
extern int gTests;
extern int gVerbose;

extern float gLastKnownScratch[2048];
extern FILE *gLastKnownFile;
extern int gLastKnownWrite;

static const unsigned int BACKEND = SoLoud::Soloud::NULLDRIVER;

#if !(defined(_MSC_VER) && !defined(__clang__))
#if defined(__i386__) || defined(__x86_64__)
#define COMPILER_BARRIER \
	do \
	{ \
		__asm__ __volatile__("" ::: "memory"); \
	} while (0)
#else
#define COMPILER_BARRIER __atomic_thread_fence(__ATOMIC_ACQUIRE)
#endif
#else
#define COMPILER_BARRIER
#endif

#ifdef __GNUC__
#define PFUNC __PRETTY_FUNCTION__
#else
#define PFUNC __FUNCTION__
#endif

#define PRINTINFO(...) gVerbose ? SoLoud::logStdout(__VA_ARGS__) : (void)0;

#define QUOTE_MACRO(x) #x
#define STRINGIZE_MACRO(x) QUOTE_MACRO(x)

#define CHECK_RES(x) \
	gTests++; \
	COMPILER_BARRIER; \
	if ((x)) \
	{ \
		gErrorCount++; \
		SoLoud::logStdout(__FILE__ ":" STRINGIZE_MACRO(__LINE__) ":%s: %s\n", PFUNC, soloud.getErrorString((x))); \
	}
#define CHECK(x) \
	gTests++; \
	COMPILER_BARRIER; \
	if (!(x)) \
	{ \
		gErrorCount++; \
		SoLoud::logStdout(__FILE__ ":" STRINGIZE_MACRO(__LINE__) ":%s: Check \"" #x "\" fail\n", PFUNC); \
	}
#define CHECK_BUF_NONZERO(x, n) \
	gTests++; \
	COMPILER_BARRIER; \
	{ \
		int i, zero = 1; \
		for (i = 0; i < (n); i++) \
			if ((x)[i] != 0) \
				zero = 0; \
		if (zero) \
		{ \
			gErrorCount++; \
			SoLoud::logStdout(__FILE__ ":" STRINGIZE_MACRO(__LINE__) ":%s: buffer not nonzero\n", PFUNC, STRINGIZE_MACRO(__LINE__), PFUNC); \
		} \
	}
#define CHECK_BUF_ZERO(x, n) \
	gTests++; \
	COMPILER_BARRIER; \
	{ \
		int i, diff = 0; \
		for (i = 0; i < (n); i++) \
			if ((x)[i] != 0) \
				diff++; \
		if (diff) \
		{ \
			gErrorCount++; \
			SoLoud::logStdout(__FILE__ ":" STRINGIZE_MACRO(__LINE__) ":%s: buffer not zero (%d / %d)\n", PFUNC, diff, (n)); \
		} \
	}
#define CHECK_BUF_DIFF(x, y, n) \
	gTests++; \
	COMPILER_BARRIER; \
	{ \
		int i, same = 1; \
		for (i = 0; i < (n); i++) \
			if (std::abs((x)[i] - (y)[i]) > 0.00001) \
				same = 0; \
		if (same) \
		{ \
			gErrorCount++; \
			SoLoud::logStdout(__FILE__ ":" STRINGIZE_MACRO(__LINE__) ":%s: buffers are equal\n", PFUNC, STRINGIZE_MACRO(__LINE__), PFUNC); \
		} \
	}
#define CHECK_BUF_SAME(x, y, n) \
	gTests++; \
	COMPILER_BARRIER; \
	{ \
		int i, diff = 0; \
		for (i = 0; i < (n); i++) \
			if (std::abs((x)[i] - (y)[i]) > 0.00001) \
				diff++; \
		if (diff) \
		{ \
			gErrorCount++; \
			SoLoud::logStdout(__FILE__ ":" STRINGIZE_MACRO(__LINE__) ":%s: buffers differ (%d / %d)\n", PFUNC, diff, (n)); \
		} \
	}
#define CHECK_BUF_SAME_LASTKNOWN(x, n) \
	gTests++; \
	COMPILER_BARRIER; \
	{ \
		int i, diff = 0, ofs = 0; \
		float maxdiff = 0.0f; \
		for (i = 0; i < (n); i++) \
		{ \
			if (std::abs((x)[i] - gLastKnownScratch[i]) > 0.00001) \
				diff++; \
			if (std::abs((x)[i] - gLastKnownScratch[i]) > maxdiff) \
			{ \
				ofs = i; \
				maxdiff = (float)std::abs((x)[i] - gLastKnownScratch[i]); \
			} \
		} \
		if (diff) \
		{ \
			gErrorCount++; \
			SoLoud::logStdout(__FILE__ ":" STRINGIZE_MACRO(__LINE__) ":%s: output differs from last known (%d / %d) maxdiff %1.5f at ofs %d\n", PFUNC, diff, (n), \
			                  maxdiff, ofs); \
		} \
	}
#define CHECK_BUF_GTE(x, y, n) \
	gTests++; \
	COMPILER_BARRIER; \
	{ \
		int i, lt = 0; \
		float absdiff = 0; \
		for (i = 0; i < (n); i++) \
		{ \
			absdiff = std::abs((x)[i]) - std::abs((y)[i]); \
			if (absdiff < 0) \
			{ \
				lt = 1; \
				if (gVerbose > 1) \
				{ \
					SoLoud::logStdout(__FILE__ ":" STRINGIZE_MACRO(__LINE__) ":%s: " #x "[%d] (%f) not bigger than buffer " #y "[%d] (%f) (abs: %f)\n", PFUNC, i, \
					                  (x)[i], i, (y)[i], absdiff); \
				} \
			} \
		} \
		if (gVerbose <= 1 && lt) \
		{ \
			gErrorCount++; \
			SoLoud::logStdout(__FILE__ ":" STRINGIZE_MACRO(__LINE__) ":%s: buffer " #x " magnitude not bigger than buffer " #y "\n", PFUNC); \
		} \
	}
#define CHECKLASTKNOWN(x, n) \
	if (gLastKnownWrite && gLastKnownFile) \
	{ \
		fwrite((x), 1, (n) * sizeof(float), gLastKnownFile); \
	} \
	else if (gLastKnownFile) \
	{ \
		fread(gLastKnownScratch, 1, (n) * sizeof(float), gLastKnownFile); \
		CHECK_BUF_SAME_LASTKNOWN((x), n); \
	}

// helper functions
void generateTestWave(SoLoud::Wav &aWav);
void generateTimingTestWave(SoLoud::Wav &aWav, unsigned int sampleCount, unsigned int sampleRate);
long getmsec();

// test functions
void testMisc();
void testGetters();
void testVis();
void testPlay();
void test3d();
void testFilters();
void testCore();
void testSpeech();
void testLoudness();
void testRelativePlaySpeedTiming();
void testMixer();
void testSpeedThings();
void benchLoudness();

// NOLINTEND(cppcoreguidelines-avoid-c-arrays, modernize-avoid-c-arrays)
