/*
Openmpt module for SoLoud audio engine
Copyright (c) 2016 Jari Komppa

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
#include <stdlib.h>
#include <math.h>

typedef void * (*dll_openmpt_module_create_from_memory)(const void * filedata, size_t filesize, void *logfunc, void * user, void * ctls);
typedef void (*dll_openmpt_module_destroy)(void * mod);
typedef int (*dll_openmpt_module_read_float_stereo)(void * mod, int samplerate, size_t count, float * left, float * right);
typedef void (*dll_openmpt_module_set_repeat_count)(void* mod, int repeat_count);

static dll_openmpt_module_create_from_memory d_openmpt_module_create_from_memory = NULL;
static dll_openmpt_module_destroy d_openmpt_module_destroy = NULL;
static dll_openmpt_module_read_float_stereo d_openmpt_module_read_float_stereo = NULL;
static dll_openmpt_module_set_repeat_count d_openmpt_module_set_repeat_count = NULL;

#ifdef WINDOWS_VERSION
#include <windows.h>

static HMODULE openDll()
{
	HMODULE res = LoadLibraryA("libopenmpt.dll");
	return res;
}

static void* getDllProc(HMODULE aDllHandle, const char *aProcName)
{
	return (void*)GetProcAddress(aDllHandle, (LPCSTR)aProcName);
}

#elif defined(__vita__)

static void * openDll()
{
	return NULL;
}

static void* getDllProc(void * aLibrary, const char *aProcName)
{
	return NULL;
}

#else
#include <dlfcn.h> // dll functions

static void * openDll()
{
	void * res;
	res = dlopen("libopenmpt.so", RTLD_LAZY);
	return res;
}

static void* getDllProc(void * aLibrary, const char *aProcName)
{
	return dlsym(aLibrary, aProcName);
}

#endif

static int load_dll()
{
#ifdef WINDOWS_VERSION
	HMODULE dll = NULL;
#else
	void * dll = NULL;
#endif

	if (d_openmpt_module_create_from_memory != NULL)
	{
		return 1;
	}

	dll = openDll();

	if (dll)
	{
		d_openmpt_module_create_from_memory = (dll_openmpt_module_create_from_memory)getDllProc(dll, "openmpt_module_create_from_memory");
		d_openmpt_module_destroy = (dll_openmpt_module_destroy)getDllProc(dll, "openmpt_module_destroy");
		d_openmpt_module_read_float_stereo = (dll_openmpt_module_read_float_stereo)getDllProc(dll, "openmpt_module_read_float_stereo");
		d_openmpt_module_set_repeat_count = (dll_openmpt_module_set_repeat_count)getDllProc(dll, "openmpt_module_set_repeat_count");



		if (d_openmpt_module_create_from_memory &&
			d_openmpt_module_destroy &&
			d_openmpt_module_read_float_stereo &&
			d_openmpt_module_set_repeat_count)
		{
			return 1;
		}
	}
	d_openmpt_module_create_from_memory = NULL;
	return 0;
}

void *openmpt_module_create_from_memory(const void * filedata, size_t filesize, void *logfunc, void * user, void * ctls)
{
	if (load_dll())
		return d_openmpt_module_create_from_memory(filedata, filesize, logfunc, user, ctls);
	return 0;
}

void openmpt_module_destroy(void * mod)
{
	if (load_dll())
		d_openmpt_module_destroy(mod);
}

int openmpt_module_read_float_stereo(void * mod, int samplerate, size_t count, float * left, float * right)
{
	if (load_dll())
		return d_openmpt_module_read_float_stereo(mod, samplerate, count, left, right);
	return 0;
}

void openmpt_module_set_repeat_count(void* mod, int repeat_count)
{
	if (load_dll())
		d_openmpt_module_set_repeat_count(mod, repeat_count);
}
