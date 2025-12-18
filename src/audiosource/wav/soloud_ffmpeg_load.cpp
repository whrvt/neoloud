/*
SoLoud audio engine - ffmpeg library loader/unloader
Copyright (c) 2025 William Horvath

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

#include <array>
#include <atomic>
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

#define VC_EXTRALEAN
#define WIN32_LEAN_AND_MEAN
#include <windows.h>

typedef HMODULE OBJHANDLE;
#define LIBLOAD(x) LoadLibraryA(x)
#define LIBFUNC(x, y) (void *)GetProcAddress(x, y)
#define LIBFREE(x) FreeLibrary(x)
// clang-format off
#define LIBGETERROR() \
	([]() -> std::string { \
		DWORD errorCode = GetLastError(); \
		if (errorCode == 0) return std::string("No error"); \
		LPSTR messageBuffer = nullptr; \
		DWORD size = FormatMessageA(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS, \
		                            NULL, errorCode, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), (LPSTR)&messageBuffer, 0, NULL); \
		if (size == 0) return std::string("Unknown error (code: ") + std::to_string(errorCode) + ")"; \
		std::string result(messageBuffer, size); \
		LocalFree(messageBuffer); \
		while (!result.empty() && (result.back() == '\n' || result.back() == '\r')) result.pop_back(); \
		return result; \
	}())
// clang-format on
#else
#include <dlfcn.h>
typedef void(*OBJHANDLE);
#define LIBLOAD(x) dlopen(x, RTLD_NOW)
#define LIBFUNC(x, y) dlsym(x, y)
#define LIBFREE(x) dlclose(x)
#define LIBGETERROR() std::string(dlerror())
#endif

#include "soloud_config.h" // for SoLoud::logStderr user-customizable output macro
#include "soloud_ffmpeg_load.h"

#ifdef _MSC_VER
#include <string.h>
#ifndef strncasecmp
#define strncasecmp _strnicmp
#endif
#endif

using namespace SoLoud::FFmpeg::FFmpegLoader::FFmpegFuncs;

namespace SoLoud::FFmpeg::FFmpegLoader
{
namespace FFmpegFuncs
{
// generate function pointer definitions
#define DEFINE_FFMPEG_FUNCTION(name) name##_t name{nullptr};
ALL_FFMPEG_FUNCTIONS(DEFINE_FFMPEG_FUNCTION)
} // namespace FFmpegFuncs

namespace
{ // anon

// library handles
OBJHANDLE s_libavutil{nullptr};
OBJHANDLE s_libswresample{nullptr};
OBJHANDLE s_libavcodec{nullptr};
OBJHANDLE s_libavformat{nullptr};

std::mutex s_init_mutex;
std::atomic s_initialized{false};
bool s_available{false};

// initialization state
std::string s_errorDetails{""};

// TODO: This should probably be consolidated and documented instead of duplicating it everywhere, but probably good enough for now.
int parse_log_level_from_env()
{
	static const bool is_debug =
#ifdef _DEBUG
	    true
#else
	    false
#endif
	    ;

	// default for debug builds, 0 (effectively disabled besides crashes) otherwise
	static const int default_level = is_debug ? av_log_get_level() : 0;

	const char *env = getenv("SOLOUD_DEBUG");
	if (!env || !*env || *env == '0')
		return default_level;

	// TODO: check if this down-shifting makes any sense in practice
	if (strncasecmp(env, "debug", sizeof("debug") - 1) == 0)
		return is_debug ? AV_LOG_VERBOSE : AV_LOG_INFO;
	if (strncasecmp(env, "info", sizeof("info") - 1) == 0)
		return is_debug ? AV_LOG_INFO : AV_LOG_WARNING;
	if (strncasecmp(env, "warn", sizeof("warn") - 1) == 0)
		return is_debug ? AV_LOG_WARNING : AV_LOG_ERROR;
	if (strncasecmp(env, "error", sizeof("error") - 1) == 0)
		return is_debug ? AV_LOG_ERROR : AV_LOG_FATAL;
	if (strncasecmp(env, "none", sizeof("none") - 1) == 0)
		return 0;

	// unknown
	return default_level;
}

// The log callback must be thread-safe, according to the ffmpeg docs
std::mutex ff_log_mutex;
int ff_log_print_prefix; // this makes no sense but is required

void ff_log_callback(void *avcl, int level, const char *fmt, va_list vl)
{
	if (level >= 0)
		level &= 0xff; // mask off "tint"

	if (level > av_log_get_level())
		return;

	std::array<char, 1024> log_line;
	std::lock_guard lk(ff_log_mutex);

	int written = av_log_format_line2(avcl, level, fmt, vl, log_line.data(), static_cast<int>(log_line.size()), &ff_log_print_prefix);
	if (written <= 0)
		return;

	// sanitize control characters that could mess with terminal
	for (char *p = log_line.data(); *p; ++p)
	{
		auto c = static_cast<unsigned char>(*p);
		if (c < 0x08 || (c > 0x0D && c < 0x20))
			*p = '?';
	}

	log_line.back() = '\0';
	SoLoud::logStderr("%s", log_line.data());
}

// assumes that the init mutex is held
void cleanup_internal()
{
#define RESET_FUNCTION(name) name = nullptr;
	ALL_FFMPEG_FUNCTIONS(RESET_FUNCTION)
#undef RESET_FUNCTION

#define RESET_LIB(libname) \
	if (!!(s_lib##libname)) \
	{ \
		LIBFREE(s_lib##libname); \
		s_lib##libname = nullptr; \
	}
	RESET_LIB(avformat);
	RESET_LIB(avcodec);
	RESET_LIB(swresample);
	RESET_LIB(avutil);
#undef RESET_LIB
}

bool init_internal()
{
// a monstrous usage of macros... but it does reduce duplication
#define LOAD_LIB_FUNCS_BODY(fname) \
	failed_count += !((fname) = reinterpret_cast<fname##_t>(LIBFUNC(current_library_outer_macro, #fname))); \
	if (!(fname)) \
		s_errorDetails += missing_prefix_outer_macro + #fname + "\n";

#define LOAD_FULL_FF_LIB(ff_libname_lower, ff_lib_version, ff_lib_funcs_xmacro) \
	[](void) -> bool { /* first load the library */ \
		               const char *trypath1 = "./" LNAME(ff_libname_lower, ff_lib_version); \
		               const char *trypath2 = LNAME(ff_libname_lower, ff_lib_version); \
		               if (!(s_lib##ff_libname_lower = LIBLOAD(trypath1)) && !(s_lib##ff_libname_lower = LIBLOAD(trypath2))) \
		               { \
			               s_errorDetails += "Failed to load " #ff_libname_lower "-" #ff_lib_version " (error: " + LIBGETERROR() + ")\n"; \
			               return false; \
		               } /* then load the functions using the given x-macro list */ \
		               const std::string missing_prefix_outer_macro = "Missing " #ff_libname_lower " function: "; /* for passing into the x-macro expansion */ \
		               auto current_library_outer_macro = s_lib##ff_libname_lower; \
		               int failed_count = 0; \
		               ff_lib_funcs_xmacro(LOAD_LIB_FUNCS_BODY); \
		               if (failed_count > 0) \
			               s_errorDetails += "Failed to load " + std::to_string(failed_count) + " " #ff_libname_lower " functions\n"; \
		               return failed_count == 0; \
	}()

	// load all libraries and functions here
	if (!LOAD_FULL_FF_LIB(avutil, 59, AVUTIL_FUNCTIONS) ||        //
	    !LOAD_FULL_FF_LIB(swresample, 5, SWRESAMPLE_FUNCTIONS) || //
	    !LOAD_FULL_FF_LIB(avcodec, 61, AVCODEC_FUNCTIONS) ||      //
	    !LOAD_FULL_FF_LIB(avformat, 61, AVFORMAT_FUNCTIONS))
	{
		cleanup_internal();
		return false;
	}

#undef LOAD_LIB_FUNCS_BODY
#undef LOAD_FULL_FF_LIB

	// we would have bailed at this point if any functions/lib failed to load

	// sanity check
	void *test_ptr = av_malloc(64);
	if (test_ptr)
	{
		av_free(test_ptr);
	}
	else
	{
		s_errorDetails += "Dysfunctional av_malloc\n";
		cleanup_internal();
		return false;
	}

	// set up log level + callback
	av_log_set_level(parse_log_level_from_env());
	av_log_set_callback(ff_log_callback);

	s_errorDetails = "";
	return true;
}
} // namespace

bool isAvailable()
{
	if (s_initialized.load(std::memory_order_relaxed))
		return s_available;

	std::lock_guard lk(s_init_mutex);
	if (s_initialized.load(std::memory_order_acquire))
		return s_available;

	s_available = init_internal();
	if (!s_available)
		SoLoud::logStderr("Failed to load FFmpeg %s\n", s_errorDetails.c_str());

	s_initialized.store(true, std::memory_order_release);

	return s_available;
}

// NOTE: currently unused
// soloud will load ffmpeg as needed but doesn't unload it
void cleanup()
{
	if (!s_initialized.load(std::memory_order_acquire))
		return;

	std::lock_guard lk(s_init_mutex);

	// reset these first
	s_available = false;
	s_initialized.store(false, std::memory_order_release);

	cleanup_internal();
}

std::string getErrorDetails()
{
	return s_errorDetails;
}

} // namespace SoLoud::FFmpeg::FFmpegLoader

#else
#pragma message("building without ffmpeg support")
#endif
