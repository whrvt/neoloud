/*
Extension to Miniaudio adding an ASIO backend.
Choice of UNLICENSE or MIT-0. See license statements at the end of this file.

This library expects a C++ compiler, as it relies on the ASIO SDK.

The Steinberg ASIO SDK is licensed under GNU GPLv3, which is a viral license unlike this library's.
Download it here: https://www.steinberg.net/developers/asiosdk-open/
*/
#if defined(MINIAUDIO_ASIO_IMPLEMENTATION) || defined(MINIAUDIO_IMPLEMENTATION) || defined(MA_IMPLEMENTATION)
#ifndef miniaudio_asio_c
#define miniaudio_asio_c

#pragma comment(lib, "Advapi32.lib")
#pragma comment(lib, "ole32.lib")
#pragma comment(lib, "User32.lib")

#include "iasiodrv.h"
#include "asiodrivers.h"

struct ma_asio_state_t
{
	ASIOBufferInfo *buffers;
	long numInputChannels;
	long numOutputChannels;
	void *interleavedPcmBuffer;

	double sampleRate;
	ma_format format;
	ma_uint32 bytesPerFrame;
	ma_uint32 periodSizeInFrames;

	ma_device *pDevice;
	ASIOCallbacks callbacks;

	ma_context *context;
	char **deviceNames;
	int nbDevices;
};

/*
    Due to limitations of the ASIO SDK, we can only load one driver at a time.
    Because of this, it's simpler to make g_asio than extend ma_device.
*/
static ma_asio_state_t g_asio;
extern IASIO *theAsioDriver;
extern AsioDrivers *asioDrivers;

static inline double asiosr_to_double(ASIOSampleRate sr)
{
	static_assert(sizeof(ASIOSampleRate) == sizeof(double), "ASIOSampleRate must be 8 bytes");
	double v;
	MA_COPY_MEMORY(&v, &sr, sizeof(v));
	return v;
}

static inline ASIOSampleRate double_to_asiosr(double v)
{
	static_assert(sizeof(ASIOSampleRate) == sizeof(double), "ASIOSampleRate must be 8 bytes");
	ASIOSampleRate sr;
	MA_COPY_MEMORY(&sr, &v, sizeof(sr));
	return sr;
}

static long ma_asio_clamp(long minSize, long maxSize, long defaultSize, long granularity, long wantedSize)
{
	if (wantedSize == -1)
		return defaultSize;
	if (wantedSize < minSize)
		return minSize;
	if (wantedSize > maxSize)
		return maxSize;
	if (granularity == 0)
		return wantedSize;

	if (granularity == -1)
	{
		/* Buffer lengths are only allowed in powers of 2 */
		for (int oksize = minSize; oksize <= maxSize; oksize *= 2)
		{
			if (oksize == wantedSize)
			{
				return wantedSize;
			}
			else if (oksize > wantedSize)
			{
				oksize /= 2;
				return oksize;
			}
		}

		/* Unreachable */
		return defaultSize;
	}
	else
	{
		/* Buffer lengths are only allowed in multiples of granularity */
		wantedSize -= minSize;
		wantedSize = (wantedSize / granularity) * granularity; /* hopefully not optimized out */
		wantedSize += minSize;
		return wantedSize;
	}
}

static void ma_asio_buffer_switch(long doubleBufferIndex, ASIOBool directProcess)
{
	(void)directProcess;

	if (!g_asio.pDevice)
		return;

	ma_device_handle_backend_data_callback(g_asio.pDevice, g_asio.interleavedPcmBuffer, NULL, g_asio.periodSizeInFrames);

	void *outputChannels[MA_MAX_CHANNELS];
	for (int i = 0; i < g_asio.numOutputChannels; i++)
	{
		outputChannels[i] = g_asio.buffers[g_asio.numInputChannels + i].buffers[doubleBufferIndex];
	}
	ma_deinterleave_pcm_frames(g_asio.format, g_asio.numOutputChannels, g_asio.periodSizeInFrames, g_asio.interleavedPcmBuffer, outputChannels);

	ASIOOutputReady();
}

static void ma_asio_sample_rate_did_change(ASIOSampleRate sRate)
{
	g_asio.sampleRate = asiosr_to_double(sRate);
	/* TODO: reset device? since miniaudio has no way to handle this... */
}

static ASIOTime *ma_asio_buffer_switch_time_info(ASIOTime *timeInfo, long doubleBufferIndex, ASIOBool directProcess)
{
	(void)timeInfo;
	ma_asio_buffer_switch(doubleBufferIndex, directProcess);
	return 0;
}

static long ma_asio_message(long selector, long value, void *message, double *opt)
{
	(void)message;
	(void)opt;

	/* TODO: handle kAsioLatenciesChanged */

	switch (selector)
	{
	case kAsioSelectorSupported:
		if (value == kAsioEngineVersion || value == kAsioResetRequest || value == kAsioSupportsTimeInfo)
			return 1;
		return 0;

	case kAsioEngineVersion:
		return 2;

	case kAsioResetRequest:
		/* TODO: reset device */
		return 0;

	case kAsioSupportsTimeInfo:
		return 1;

	default:
		return 0;
	}
}

static ma_result ma_device_start__asio(ma_device *pDevice)
{
	MA_ASSERT(pDevice != NULL);
	(void)pDevice;

	/* Fill output buffers with null bytes before ASIOStart() */
	for (int i = 0; i < g_asio.numInputChannels + g_asio.numOutputChannels; i++)
	{
		if (g_asio.buffers[i].isInput)
			continue;

		void *bufA = g_asio.buffers[i].buffers[0];
		if (bufA)
			memset(bufA, 0, g_asio.periodSizeInFrames * g_asio.bytesPerFrame);

		void *bufB = g_asio.buffers[i].buffers[1];
		if (bufB)
			memset(bufB, 0, g_asio.periodSizeInFrames * g_asio.bytesPerFrame);
	}

	ASIOError err = ASIOStart();
	if (err == ASE_OK)
	{
		return MA_SUCCESS;
	}

	ma_log_postf(ma_device_get_log(pDevice), MA_LOG_LEVEL_ERROR, "[ASIO] ASIOStart() failed with error %ld\n", err);
	return MA_ERROR;
}

static ma_result ma_device_uninit__asio(ma_device *pDevice)
{
	MA_ASSERT(pDevice != NULL);
	(void)pDevice;

	if (theAsioDriver)
	{
		ASIOStop();
		ASIODisposeBuffers();
		ASIOExit();
	}

	free(g_asio.buffers);
	g_asio.buffers = NULL;

	free(g_asio.interleavedPcmBuffer);
	g_asio.interleavedPcmBuffer = NULL;

	g_asio.pDevice = NULL;

	return MA_SUCCESS;
}

static ma_result ma_device_stop__asio(ma_device *pDevice)
{
	MA_ASSERT(pDevice != NULL);
	(void)pDevice;

	return (ASIOStop() == ASE_OK) ? MA_SUCCESS : MA_ERROR;
}

static ma_result ma_device_init__asio(ma_device *pDevice,
                                      const ma_device_config *pConfig,
                                      ma_device_descriptor *pDescriptorPlayback,
                                      ma_device_descriptor *pDescriptorCapture)
{
	MA_ASSERT(pDevice != NULL);
	MA_ASSERT(pDescriptorPlayback != NULL);
	(void)pDescriptorCapture;

	long minSize = 0, maxSize = 0, preferredSize = 0, granularity = 0;
	ASIOError err;
	ASIOSampleRate asioSR;
	size_t numChannels;

	ASIODriverInfo drvInfo = {};
	drvInfo.sysRef = g_asio.context->dsound.hWnd;

	ASIOChannelInfo channelInfo = {};
	channelInfo.channel = 0;
	channelInfo.isInput = ASIOFalse;

	const char *driverName = g_asio.deviceNames[0];
	if (pDescriptorPlayback->pDeviceID != NULL)
	{
		driverName = (const char *)pDescriptorPlayback->pDeviceID->custom.p;
		MA_ASSERT(driverName != NULL);
	}

	/* loadDriver takes in a char*, but doesn't modify it */
	if (!asioDrivers->loadDriver((char *)driverName))
	{
		goto error;
	}

	/* Hacky to pass window handle in dsound, but it's optional anyway */
	if (ASIOInit(&drvInfo) != ASE_OK)
	{
		goto error;
	}

	if (ASIOGetChannels(&g_asio.numInputChannels, &g_asio.numOutputChannels) != ASE_OK)
	{
		goto error;
	}

	if (ASIOGetBufferSize(&minSize, &maxSize, &preferredSize, &granularity) != ASE_OK)
	{
		goto error;
	}
	g_asio.periodSizeInFrames = preferredSize;
	if (pConfig->periodSizeInFrames > 0)
	{
		g_asio.periodSizeInFrames = ma_asio_clamp(minSize, maxSize, preferredSize, granularity, pConfig->periodSizeInFrames);
	}

	/* Some devices will accept 32khz but die until they are set back to 44.1khz... */
	asioSR = double_to_asiosr(44100.0);
	if (pConfig->sampleRate > 44100)
	{
		asioSR = double_to_asiosr((double)pConfig->sampleRate);
	}
	if (ASIOSetSampleRate(asioSR) != ASE_OK)
	{
		asioSR = double_to_asiosr(44100.0);
		if (ASIOSetSampleRate(asioSR) != ASE_OK)
		{
			ma_log_postf(ma_device_get_log(pDevice), MA_LOG_LEVEL_ERROR, "[ASIO] Failed to set sample rate!\n");
			goto error;
		}
	}
	g_asio.sampleRate = asiosr_to_double(asioSR);
	if (ASIOGetSampleRate(&asioSR) == ASE_OK && asiosr_to_double(asioSR) > 0.0)
	{
		/* Some devices will return a sample rate of 0, don't trust ASIOGetSampleRate. */
		g_asio.sampleRate = asiosr_to_double(asioSR);
	}
	ma_log_postf(ma_device_get_log(pDevice), MA_LOG_LEVEL_INFO, "[ASIO] Current sample rate: %u Hz\n", (ma_uint32)g_asio.sampleRate);

	/* TODO: figure out which of these we actually support, and fix byte order for the ones we don't yet */
	ASIOGetChannelInfo(&channelInfo);
	switch (channelInfo.type)
	{
	case ASIOSTInt16LSB:
	case ASIOSTInt16MSB:
		g_asio.bytesPerFrame = 2;
		g_asio.format = ma_format_s16;
		break;

	case ASIOSTInt24LSB:
	case ASIOSTInt24MSB:
		g_asio.bytesPerFrame = 3;
		g_asio.format = ma_format_s24;
		break;

	case ASIOSTInt32LSB:
	case ASIOSTInt32MSB:
	case ASIOSTInt32LSB16:
	case ASIOSTInt32LSB18:
	case ASIOSTInt32LSB20:
	case ASIOSTInt32LSB24:
	case ASIOSTInt32MSB16:
	case ASIOSTInt32MSB18:
	case ASIOSTInt32MSB20:
	case ASIOSTInt32MSB24:
		g_asio.bytesPerFrame = 4;
		g_asio.format = ma_format_s32;
		break;

	case ASIOSTFloat32LSB:
	case ASIOSTFloat32MSB:
		g_asio.bytesPerFrame = 4;
		g_asio.format = ma_format_f32;
		break;

	default:
		/* f64/DSD not supported */
		ma_log_post(ma_device_get_log(pDevice), MA_LOG_LEVEL_ERROR, "[ASIO] Unsupported output buffer format\n");
		goto error;
	}

	numChannels = g_asio.numInputChannels + g_asio.numOutputChannels;
	g_asio.buffers = (ASIOBufferInfo *)(malloc(sizeof(ASIOBufferInfo) * numChannels));
	if (!g_asio.buffers)
	{
		goto error;
	}
	for (int i = 0; i < g_asio.numInputChannels; i++)
	{
		g_asio.buffers[i].isInput = ASIOTrue;
		g_asio.buffers[i].channelNum = i;
		g_asio.buffers[i].buffers[0] = NULL;
		g_asio.buffers[i].buffers[1] = NULL;
	}
	for (int i = 0; i < g_asio.numOutputChannels; i++)
	{
		g_asio.buffers[g_asio.numInputChannels + i].isInput = ASIOFalse;
		g_asio.buffers[g_asio.numInputChannels + i].channelNum = i;
		g_asio.buffers[g_asio.numInputChannels + i].buffers[0] = NULL;
		g_asio.buffers[g_asio.numInputChannels + i].buffers[1] = NULL;
	}

	g_asio.callbacks.bufferSwitch = ma_asio_buffer_switch;
	g_asio.callbacks.sampleRateDidChange = ma_asio_sample_rate_did_change;
	g_asio.callbacks.asioMessage = ma_asio_message;
	g_asio.callbacks.bufferSwitchTimeInfo = ma_asio_buffer_switch_time_info;

	err = ASIOCreateBuffers(g_asio.buffers, numChannels, g_asio.periodSizeInFrames, &g_asio.callbacks);
	if (err != ASE_OK)
	{
		ma_log_postf(ma_device_get_log(pDevice), MA_LOG_LEVEL_ERROR, "[ASIO] ASIOCreateBuffers failed (code %ld)\n", err);
		free(g_asio.buffers);
		goto error;
	}

	g_asio.interleavedPcmBuffer = malloc(g_asio.periodSizeInFrames * g_asio.bytesPerFrame * g_asio.numOutputChannels);
	if (!g_asio.interleavedPcmBuffer)
	{
		ASIODisposeBuffers();
		free(g_asio.buffers);
		goto error;
	}

	pDescriptorPlayback->format = g_asio.format;
	pDescriptorPlayback->channels = g_asio.numOutputChannels;
	pDescriptorPlayback->sampleRate = g_asio.sampleRate;
	pDescriptorPlayback->periodSizeInFrames = g_asio.periodSizeInFrames;
	ma_channel_map_init_standard(
	    ma_standard_channel_map_default, pDescriptorPlayback->channelMap, ma_countof(pDescriptorPlayback->channelMap), g_asio.numOutputChannels);

	g_asio.pDevice = pDevice;
	return MA_SUCCESS;

error:
	ASIOExit();
	return MA_ERROR;
}

static ma_result ma_context_get_device_info__asio(ma_context *pContext, ma_device_type type, const ma_device_id *pDeviceID, ma_device_info *pDeviceInfo)
{
	MA_ASSERT(pContext != NULL);
	MA_ASSERT(pDeviceInfo != NULL);
	(void)pDeviceID;
	(void)type; /* always playback */

	/*
	    To query device info, we need to load the driver.
	    Which can be an issue with ASIO only allowing one driver at a time...
	    So we restrict the API to only allow get_device_info on the currently loaded device.
	*/
	if (g_asio.pDevice == NULL)
	{
		return MA_NO_DEVICE;
	}

	MA_ZERO_OBJECT(pDeviceInfo);
	if (!asioDrivers->getCurrentDriverName(pDeviceInfo->name))
	{
		return MA_NO_DEVICE;
	}
	pDeviceInfo->nativeDataFormatCount = 1;
	pDeviceInfo->nativeDataFormats[0].format = g_asio.format;
	pDeviceInfo->nativeDataFormats[0].channels = g_asio.numOutputChannels;
	pDeviceInfo->nativeDataFormats[0].sampleRate = g_asio.sampleRate;
	pDeviceInfo->nativeDataFormats[0].flags = 0;

	return MA_SUCCESS;
}

static ma_result ma_context_enumerate_devices__asio(ma_context *pContext, ma_enum_devices_callback_proc callback, void *pUserData)
{
	MA_ASSERT(pContext != NULL);

	for (int i = 0; i < g_asio.nbDevices; i++)
	{
		ma_device_info info;
		MA_ZERO_OBJECT(&info);

		info.id.custom.p = g_asio.deviceNames[i];
		ma_strcpy_s(info.name, sizeof(info.name), g_asio.deviceNames[i]);
		info.isDefault = MA_FALSE;

		if (!callback(pContext, ma_device_type_playback, &info, pUserData))
		{
			return MA_ERROR;
		}
	}

	return MA_SUCCESS;
}

static ma_result ma_context_uninit__asio(ma_context *pContext)
{
	MA_ASSERT(pContext != NULL);

	for (int i = 0; i < g_asio.nbDevices; i++)
	{
		free(g_asio.deviceNames[i]);
	}
	free(g_asio.deviceNames);
	g_asio.deviceNames = NULL;
	g_asio.nbDevices = 0;

	delete asioDrivers;
	asioDrivers = NULL;

	return MA_SUCCESS;
}

static ma_result ma_context_init__asio(ma_context *pContext, const ma_context_config *pConfig, ma_backend_callbacks *pCallbacks)
{
	if (asioDrivers != NULL)
	{
		return MA_ALREADY_EXISTS;
	}
	asioDrivers = new AsioDrivers();
	g_asio.context = pContext;

	g_asio.nbDevices = asioDrivers->asioGetNumDev();
	if (g_asio.nbDevices == 0)
	{
		if (pConfig && pConfig->pLog)
			ma_log_post(pConfig->pLog, MA_LOG_LEVEL_WARNING, "No ASIO drivers found.\n");
		goto end;
	}

	g_asio.deviceNames = (char **)malloc(g_asio.nbDevices * sizeof(char *));
	for (int i = 0; i < g_asio.nbDevices; i++)
	{
		g_asio.deviceNames[i] = (char *)malloc(MA_MAX_DEVICE_NAME_LENGTH);
		asioDrivers->asioGetDriverName(i, g_asio.deviceNames[i], MA_MAX_DEVICE_NAME_LENGTH);
	}

end:
	pCallbacks->onContextInit = ma_context_init__asio;
	pCallbacks->onContextUninit = ma_context_uninit__asio;
	pCallbacks->onContextEnumerateDevices = ma_context_enumerate_devices__asio;
	pCallbacks->onContextGetDeviceInfo = ma_context_get_device_info__asio;
	pCallbacks->onDeviceInit = ma_device_init__asio;
	pCallbacks->onDeviceUninit = ma_device_uninit__asio;
	pCallbacks->onDeviceStart = ma_device_start__asio;
	pCallbacks->onDeviceStop = ma_device_stop__asio;
	return MA_SUCCESS;
}

#endif /* miniaudio_asio_c */
#endif

/*
This software is available as a choice of the following licenses. Choose
whichever you prefer.

===============================================================================
ALTERNATIVE 1 - Public Domain (www.unlicense.org)
===============================================================================
This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or distribute this
software, either in source code form or as a compiled binary, for any purpose,
commercial or non-commercial, and by any means.

In jurisdictions that recognize copyright laws, the author or authors of this
software dedicate any and all copyright interest in the software to the public
domain. We make this dedication for the benefit of the public at large and to
the detriment of our heirs and successors. We intend this dedication to be an
overt act of relinquishment in perpetuity of all present and future rights to
this software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <http://unlicense.org/>

===============================================================================
ALTERNATIVE 2 - MIT No Attribution
===============================================================================
Copyright 2026 kiwec

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
