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

#ifndef SOLOUD_CONFIG_H
#define SOLOUD_CONFIG_H

#ifdef SOLOUD_NO_ASSERTS
#define SOLOUD_ASSERT(x)
#else
#ifdef _MSC_VER
#include <stdio.h> // for sprintf in asserts
#ifndef VC_EXTRALEAN
#define VC_EXTRALEAN
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h> // only needed for OutputDebugStringA, should be solved somehow.
#define SOLOUD_ASSERT(x) \
	if (!(x)) \
	{ \
		char temp[200]; \
		sprintf(temp, "%s(%d): assert(%s) failed.\n", __FILE__, __LINE__, #x); \
		OutputDebugStringA(temp); \
		__debugbreak(); \
	}
#else
#include <assert.h> // assert
#define SOLOUD_ASSERT(x) assert(x)
#endif
#endif

#ifndef M_PI
#define M_PI 3.14159265359
#endif

#if !defined(DISABLE_SIMD)
#if defined(__SSE2__)
#define SOLOUD_SSE_INTRINSICS
#endif
#endif

#define SOLOUD_VERSION 202506

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
// Configuration defines

// Maximum number of filters per stream
#ifndef FILTERS_PER_STREAM
#define FILTERS_PER_STREAM 8
#endif

// Number of samples to process on one go
#ifndef SAMPLE_GRANULARITY
#define SAMPLE_GRANULARITY 512
#endif

// Maximum number of concurrent voices (hard limit is 4095)
#ifndef VOICE_COUNT
#define VOICE_COUNT 1024
#endif

// 1)mono, 2)stereo 4)quad 6)5.1 8)7.1
#ifndef MAX_CHANNELS
#define MAX_CHANNELS 8
#endif

// Default resampler for both main and bus mixers
#ifndef SOLOUD_DEFAULT_RESAMPLER
#define SOLOUD_DEFAULT_RESAMPLER SoLoud::Soloud::RESAMPLER_LINEAR
#endif

//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

// Typedefs have to be made before the includes, as the
// includes depend on them.
namespace SoLoud
{
	class Soloud;
	typedef void (*mutexCallFunction)(void *aMutexPtr);
	typedef void (*soloudCallFunction)(Soloud *aSoloud);
	typedef unsigned int result;
	typedef result (*soloudResultFunction)(Soloud *aSoloud);
	typedef unsigned int handle;
	typedef double time;
}; // namespace SoLoud

#endif
