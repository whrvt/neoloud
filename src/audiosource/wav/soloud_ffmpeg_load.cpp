/*
SoLoud audio engine - ffmpeg library loader/unloader
Copyright (c) 2013-2020 Jari Komppa
Copyright (c) 2025 William Horvath (ffmpeg interface)

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

#if defined(WITH_FFMPEG) && __has_include(<libavcodec/avcodec.h>) && ((defined(_WIN32) || defined(_WIN64)) || defined(__linux__))
#pragma message("building with ffmpeg support")

#include <mutex>
#include <string>

#if defined(_WIN32) || defined(_WIN64)
#define LNAME(x, ver) #x "-" #ver ".dll"
#else
#define LNAME(x, ver) "lib" #x ".so." #ver
#endif

#ifdef WITH_SDL3
#include <SDL3/SDL_loadso.h>
typedef SDL_SharedObject(*OBJHANDLE);
#define LIBLOAD(x) SDL_LoadObject(x)
#define LIBFUNC(x, y) SDL_LoadFunction(x, y)
#define LIBFREE(x) SDL_UnloadObject(x)
#define LIBGETERROR() std::string(SDL_GetError())
#elif defined(_WIN32) || defined(_WIN64)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef NOWINRES
#define NOWINRES
#endif
#ifndef NOSERVICE
#define NOSERVICE
#endif
#ifndef NOMCX
#define NOMCX
#endif
#ifndef NOCRYPT
#define NOCRYPT
#endif
#ifndef NOMETAFILE
#define NOMETAFILE
#endif
#ifndef MMNOSOUND
#define MMNOSOUND
#endif
#ifndef VC_EXTRALEAN
#define VC_EXTRALEAN
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
typedef HMODULE OBJHANDLE;
#define LIBLOAD(x) LoadLibraryA(x)
#define LIBFUNC(x, y) (void *)GetProcAddress(x, y)
#define LIBFREE(x) FreeLibrary(x)
#define LIBGETERROR() \
	([]() -> std::string { \
		DWORD errorCode = GetLastError(); \
		if (errorCode == 0) \
			return std::string("No error"); \
		LPSTR messageBuffer = nullptr; \
		DWORD size = FormatMessageA(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS, NULL, errorCode, \
		                            MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), (LPSTR) & messageBuffer, 0, NULL); \
		if (size == 0) \
			return std::string("Unknown error (code: ") + std::to_string(errorCode) + ")"; \
		std::string result(messageBuffer, size); \
		LocalFree(messageBuffer); \
		while (!result.empty() && (result.back() == '\n' || result.back() == '\r')) \
			result.pop_back(); \
		return result; \
	}())
#else
#include <dlfcn.h>
typedef void(*OBJHANDLE);
#define LIBLOAD(x) dlopen(x, RTLD_NOW)
#define LIBFUNC(x, y) dlsym(x, y)
#define LIBFREE(x) dlclose(x)
#define LIBGETERROR() std::string(dlerror())
#endif

#include "soloud_ffmpeg_load.h"

using namespace SoLoud::FFmpeg::FFmpegLoader::FFmpegFuncs;

namespace SoLoud::FFmpeg::FFmpegLoader
{

// library handles
static OBJHANDLE s_libAvutil = nullptr;
static OBJHANDLE s_libSwresample = nullptr;
static OBJHANDLE s_libAvcodec = nullptr;
static OBJHANDLE s_libAvformat = nullptr;

static std::once_flag s_init_flag;
static bool s_init_success{false};
static bool s_available{false};

// initialization state
static std::string s_errorDetails = "";

namespace FFmpegFuncs
{
// generate function pointer definitions
#define DEFINE_FFMPEG_FUNCTION(name) name##_t name = nullptr;
ALL_FFMPEG_FUNCTIONS(DEFINE_FFMPEG_FUNCTION)
} // namespace FFmpegFuncs

template <typename T>
static T loadFunction(OBJHANDLE lib, const char *funcName)
{
	if (!lib)
		return nullptr;
	return reinterpret_cast<T>(LIBFUNC(lib, funcName));
}

static bool loadAvutilFunctions()
{
	if (!s_libAvutil)
	{
		s_errorDetails += "avutil library not loaded\n";
		return false;
	}

	bool success = true;
	int failed_count = 0;

#define LOAD_AVUTIL_FUNCTION(name) \
	name = loadFunction<name##_t>(s_libAvutil, #name); \
	if (!(name)) \
	{ \
		s_errorDetails += std::string("Missing avutil function: ") + #name + "\n"; \
		failed_count++; \
		success = false; \
	}

	AVUTIL_FUNCTIONS(LOAD_AVUTIL_FUNCTION)

	if (!success)
		s_errorDetails += "Failed to load " + std::to_string(failed_count) + " avutil functions\n";

	return success;
}

static bool loadAvcodecFunctions()
{
	if (!s_libAvcodec)
	{
		s_errorDetails += "avcodec library not loaded\n";
		return false;
	}

	bool success = true;
	int failed_count = 0;

#define LOAD_AVCODEC_FUNCTION(name) \
	name = loadFunction<name##_t>(s_libAvcodec, #name); \
	if (!(name)) \
	{ \
		s_errorDetails += std::string("Missing avcodec function: ") + #name + "\n"; \
		failed_count++; \
		success = false; \
	}

	AVCODEC_FUNCTIONS(LOAD_AVCODEC_FUNCTION)

	if (!success)
		s_errorDetails += "Failed to load " + std::to_string(failed_count) + " avcodec functions\n";

	return success;
}

static bool loadAvformatFunctions()
{
	if (!s_libAvformat)
	{
		s_errorDetails += "avformat library not loaded\n";
		return false;
	}

	bool success = true;
	int failed_count = 0;

#define LOAD_AVFORMAT_FUNCTION(name) \
	name = loadFunction<name##_t>(s_libAvformat, #name); \
	if (!(name)) \
	{ \
		s_errorDetails += std::string("Missing avformat function: ") + #name + "\n"; \
		failed_count++; \
		success = false; \
	}

	AVFORMAT_FUNCTIONS(LOAD_AVFORMAT_FUNCTION)

	if (!success)
		s_errorDetails += "Failed to load " + std::to_string(failed_count) + " avformat functions\n";

	return success;
}

static bool loadSwresampleFunctions()
{
	if (!s_libSwresample)
	{
		s_errorDetails += "swresample library not loaded\n";
		return false;
	}

	bool success = true;
	int failed_count = 0;

#define LOAD_SWRESAMPLE_FUNCTION(name) \
	name = loadFunction<name##_t>(s_libSwresample, #name); \
	if (!(name)) \
	{ \
		s_errorDetails += std::string("Missing swresample function: ") + #name + "\n"; \
		failed_count++; \
		success = false; \
	}

	SWRESAMPLE_FUNCTIONS(LOAD_SWRESAMPLE_FUNCTION)

	if (!success)
		s_errorDetails += "Failed to load " + std::to_string(failed_count) + " swresample functions\n";

	return success;
}

// silence, you fool
static void empty_log_callback(void * /**/, int /**/, const char * /**/, va_list /**/)
{
	;
}

static bool init_locked()
{
	s_available = false;
	s_errorDetails = "";

	cleanup();

	s_libAvutil = LIBLOAD("./" LNAME(avutil, 59));
	if (!(s_libAvutil = LIBLOAD(LNAME(avutil, 59))))
	{
		s_errorDetails = "Failed to load libavutil-59 (error: " + LIBGETERROR() + ")";
		return false;
	}

	s_libSwresample = LIBLOAD("./" LNAME(swresample, 5));
	if (!(s_libSwresample = LIBLOAD(LNAME(swresample, 5))))
	{
		s_errorDetails = "Failed to load libswresample-5 (error: " + LIBGETERROR() + ")";
		cleanup();
		return false;
	}

	s_libAvcodec = LIBLOAD("./" LNAME(avcodec, 61));
	if (!(s_libAvcodec = LIBLOAD(LNAME(avcodec, 61))))
	{
		s_errorDetails = "Failed to load libavcodec-61 (error: " + LIBGETERROR() + ")";
		cleanup();
		return false;
	}

	s_libAvformat = LIBLOAD("./" LNAME(avformat, 61));
	if (!(s_libAvformat = LIBLOAD(LNAME(avformat, 61))))
	{
		s_errorDetails = "Failed to load libavformat-61 (error: " + LIBGETERROR() + ")";
		cleanup();
		return false;
	}

	bool functionsLoaded = true;
	functionsLoaded &= loadAvutilFunctions();
	functionsLoaded &= loadSwresampleFunctions();
	functionsLoaded &= loadAvcodecFunctions();
	functionsLoaded &= loadAvformatFunctions();

	if (!functionsLoaded)
	{
		cleanup();
		return false;
	}

#ifndef _DEBUG
	if (av_log_set_callback)
		av_log_set_callback(empty_log_callback);
#endif

	// verify basic functionality
	if (av_malloc && av_free)
	{
		void *test_ptr = av_malloc(64);
		if (test_ptr)
		{
			av_free(test_ptr);
		}
	}
	else
	{
		s_errorDetails += "Critical functions av_malloc/av_free not loaded\n";
		cleanup();
		return false;
	}

	s_available = true;
	s_errorDetails = "";
	return true;
}

bool init()
{
	std::call_once(s_init_flag, []() { s_init_success = init_locked(); });

	return s_init_success;
}

void cleanup()
{
	s_available = false;

#define RESET_FUNCTION(name) name = nullptr;
	ALL_FFMPEG_FUNCTIONS(RESET_FUNCTION)

	if (s_libAvformat)
	{
		LIBFREE(s_libAvformat);
		s_libAvformat = nullptr;
	}
	if (s_libAvcodec)
	{
		LIBFREE(s_libAvcodec);
		s_libAvcodec = nullptr;
	}
	if (s_libSwresample)
	{
		LIBFREE(s_libSwresample);
		s_libSwresample = nullptr;
	}
	if (s_libAvutil)
	{
		LIBFREE(s_libAvutil);
		s_libAvutil = nullptr;
	}
}

bool isAvailable()
{
	return s_available;
}

std::string getErrorDetails()
{
	return s_errorDetails;
}

} // namespace SoLoud::FFmpeg::FFmpegLoader

#else
#pragma message("building without ffmpeg support")
#endif
