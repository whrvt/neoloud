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
#include "soloud.h"

#if !defined(WITH_MINIAUDIO)

namespace SoLoud
{
result miniaudio_init(SoLoud::Soloud *aSoloud, unsigned int aFlags, unsigned int aSamplerate, unsigned int aBuffer)
{
	return NOT_IMPLEMENTED;
}
} // namespace SoLoud

#else

#define MINIAUDIO_IMPLEMENTATION

#ifdef __linux__
#define MA_DEFAULT_PERIOD_SIZE_IN_MILLISECONDS_LOW_LATENCY 5 // smaller period size results in less crackling when auto-negotiated? idk why
#endif

#if defined(_WIN32) || defined(_WIN64) || defined(_MSC_VER)
#define MA_COINIT_VALUE \
	0x2 // apartment-threaded (there's some weird workaround that miniaudio does with multithreaded init which i'm not sure is required for our use case)
#define MA_DEFAULT_PERIOD_SIZE_IN_MILLISECONDS_LOW_LATENCY 1 // this seems to safely auto-renegotiate to the lowest period size anyways
#endif

// disable unneeded miniaudio features, we just want playback functionality
#define MA_NO_NULL // no null audio backend
#define MA_NO_DECODING
#define MA_NO_ENCODING
#define MA_NO_WAV
#define MA_NO_FLAC
#define MA_NO_MP3
#define MA_NO_RESOURCE_MANAGER
#define MA_NO_NODE_GRAPH
#define MA_NO_ENGINE

#include "miniaudio.h"

#include "soloud_internal.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace SoLoud
{
struct MiniaudioData
{
	ma_context context{};
	ma_device device{};
	ma_log log{};
	bool contextInitialized{false};
	bool deviceInitialized{false};
	bool logInitialized{false};
	ma_uint32 maxLogLevel{0};

	MiniaudioData() = default;
};

namespace // static
{
void soloud_miniaudio_log_callback(void *pUserData, ma_uint32 level, const char *pMessage)
{
	ma_uint32 maxLevel = *static_cast<ma_uint32 *>(pUserData);

	if (level > maxLevel)
		return; // don't log if level is more verbose than our max

	const char *levelStr = nullptr;
	switch (level)
	{
	case MA_LOG_LEVEL_DEBUG:
		levelStr = "DEBUG";
		break;
	case MA_LOG_LEVEL_INFO:
		levelStr = "INFO";
		break;
	case MA_LOG_LEVEL_WARNING:
		levelStr = "WARNING";
		break;
	case MA_LOG_LEVEL_ERROR:
		levelStr = "ERROR";
		break;
	default:
		levelStr = "UNKNOWN";
		break;
	}

	fprintf(stderr, "[MiniAudio %s] %s", levelStr, pMessage); // it seems that the log messages are pre-newline-terminated
}

ma_uint32 parse_log_level_from_env()
{
	static const char *env = getenv("MINIAUDIO_DEBUG");
	if (!env)
	{
#ifdef _DEBUG
		return MA_LOG_LEVEL_WARNING;
#else
		return 0; // none
#endif
	}

	if (strcmp(env, "debug") == 0)
		return MA_LOG_LEVEL_DEBUG;
	if (strcmp(env, "info") == 0)
		return MA_LOG_LEVEL_INFO;
	if (strcmp(env, "warn") == 0)
		return MA_LOG_LEVEL_WARNING;
	if (strcmp(env, "error") == 0)
		return MA_LOG_LEVEL_ERROR;
	if (strcmp(env, "none") == 0)
		return 0;

	// default fallback for unrecognized values
#ifdef _DEBUG
	return MA_LOG_LEVEL_WARNING;
#else
	return 0;
#endif
}

void soloud_miniaudio_audiomixer(ma_device *pDevice, void *pOutput, const void * /*pInput*/, ma_uint32 frameCount)
{
	auto *soloud = static_cast<SoLoud::Soloud *>(pDevice->pUserData);

	switch (pDevice->playback.internalFormat)
	{
	case ma_format_f32:
		soloud->mix(static_cast<float *>(pOutput), frameCount);
		break;
	case ma_format_s16:
		soloud->mixSigned16(static_cast<short *>(pOutput), frameCount);
		break;
	case ma_format_u8:
		soloud->mixUnsigned8(static_cast<unsigned char *>(pOutput), frameCount);
		break;
	case ma_format_s24:
		soloud->mixSigned24(static_cast<unsigned char *>(pOutput), frameCount);
		break;
	case ma_format_s32:
		soloud->mixSigned32(static_cast<int *>(pOutput), frameCount);
		break;
	default:
		// fallback to float if format is unknown or unsupported
		soloud->mix(static_cast<float *>(pOutput), frameCount);
		break;
	}
}

void soloud_miniaudio_deinit(SoLoud::Soloud *aSoloud)
{
	auto *data = static_cast<MiniaudioData *>(aSoloud->mBackendData);
	if (data)
	{
		if (data->deviceInitialized)
		{
			ma_device_uninit(&data->device);
			data->deviceInitialized = false;
		}
		if (data->contextInitialized)
		{
			ma_context_uninit(&data->context);
			data->contextInitialized = false;
		}
		if (data->logInitialized)
		{
			ma_log_uninit(&data->log);
			data->logInitialized = false;
		}
		delete data;
		aSoloud->mBackendData = nullptr;
	}
}

result soloud_miniaudio_pause(SoLoud::Soloud *aSoloud)
{
	auto *data = static_cast<MiniaudioData *>(aSoloud->mBackendData);
	if (data && data->deviceInitialized)
	{
		if (ma_device_stop(&data->device) != MA_SUCCESS)
			return UNKNOWN_ERROR;
	}
	return SO_NO_ERROR;
}

result soloud_miniaudio_resume(SoLoud::Soloud *aSoloud)
{
	auto *data = static_cast<MiniaudioData *>(aSoloud->mBackendData);
	if (data && data->deviceInitialized)
	{
		if (ma_device_start(&data->device) != MA_SUCCESS)
			return UNKNOWN_ERROR;
	}
	return SO_NO_ERROR;
}

} // namespace

result miniaudio_init(SoLoud::Soloud *aSoloud, unsigned int aFlags, unsigned int aSamplerate, unsigned int aBuffer, unsigned int aChannels)
{
	auto *data = new MiniaudioData();
	aSoloud->mBackendData = data;

	// setup logging
	data->maxLogLevel = parse_log_level_from_env();
	if (data->maxLogLevel > 0)
	{
		ma_result logResult = ma_log_init(nullptr, &data->log);
		if (logResult == MA_SUCCESS)
		{
			data->logInitialized = true;
			ma_log_callback callback = ma_log_callback_init(soloud_miniaudio_log_callback, &data->maxLogLevel);
			ma_log_register_callback(&data->log, callback);
		}
	}

	ma_context_config contextConfig = ma_context_config_init();
	contextConfig.threadPriority = ma_thread_priority_highest;
	if (data->logInitialized)
		contextConfig.pLog = &data->log;

	// this can be set+reordered if we want a different priority
	// const ma_backend backends[] = {ma_backend_wasapi, ma_backend_dsound, ma_backend_winmm,    ma_backend_coreaudio,  ma_backend_sndio,
	//                                ma_backend_audio4, ma_backend_oss,    ma_backend_alsa,     ma_backend_pulseaudio, ma_backend_jack,
	//                                ma_backend_aaudio, ma_backend_opensl, ma_backend_webaudio, ma_backend_custom,     ma_backend_null};

	// ma_result result = ma_context_init(&backends[0], sizeof(backends) / sizeof((backends)[0]), &contextConfig, &data->context);
	ma_result result = ma_context_init(nullptr, 0, &contextConfig, &data->context);
	if (result != MA_SUCCESS)
	{
		if (data->logInitialized)
			ma_log_uninit(&data->log);
		delete data;
		aSoloud->mBackendData = nullptr;
		switch (result)
		{
		case MA_INVALID_ARGS:
			return INVALID_PARAMETER;
		case MA_OUT_OF_MEMORY:
			return OUT_OF_MEMORY;
		case MA_NO_BACKEND:
			return NOT_IMPLEMENTED;
		default:
			return UNKNOWN_ERROR;
		}
	}
	data->contextInitialized = true;

	// configure device
	ma_device_config config = ma_device_config_init(ma_device_type_playback);
	config.playback.format = ma_format_unknown; // default device format, we'll do the conversion ourselves
	config.playback.channels = aChannels;
	config.dataCallback = soloud_miniaudio_audiomixer;
	config.pUserData = (void *)aSoloud;
	config.noPreSilencedOutputBuffer = true;
	config.noClip = true;
	config.performanceProfile = ma_performance_profile_low_latency;

	if (aSamplerate > 0) // respect miniaudio default (avoids extra resampling if we use the device native sample rate)
		config.sampleRate = aSamplerate;
	else
		config.sampleRate = 0;

	if (aBuffer > 0) // respect miniaudio default
		config.periodSizeInFrames = aBuffer;
	else
	{
		config.periodSizeInFrames = 0;
#if defined(_WIN32) || defined(_WIN64) || defined(_MSC_VER)
		config.periodSizeInMilliseconds = 1; // probably redundant with the define above
#endif
	}

	// backend-specific settings
	config.wasapi.noAutoConvertSRC = true; // soloud handles resampling
	config.wasapi.noDefaultQualitySRC = true;
	config.alsa.noAutoFormat = true;
	config.alsa.noAutoChannels = true;
	config.alsa.noAutoResample = true;

	result = ma_device_init(&data->context, &config, &data->device);
	if (result != MA_SUCCESS)
	{
		ma_context_uninit(&data->context);
		if (data->logInitialized)
			ma_log_uninit(&data->log);
		delete data;
		aSoloud->mBackendData = nullptr;

		switch (result)
		{
		case MA_INVALID_ARGS:
			return INVALID_PARAMETER;
		case MA_OUT_OF_MEMORY:
			return OUT_OF_MEMORY;
		case MA_FORMAT_NOT_SUPPORTED:
			return INVALID_PARAMETER;
		case MA_DEVICE_NOT_INITIALIZED:
			// return UNKNOWN_ERROR;
		case MA_DEVICE_ALREADY_INITIALIZED:
			// return UNKNOWN_ERROR;
		case MA_DEVICE_NOT_STARTED:
			// return UNKNOWN_ERROR;
		case MA_DEVICE_NOT_STOPPED:
			// return UNKNOWN_ERROR;
		default:
			return UNKNOWN_ERROR;
		}
	}
	data->deviceInitialized = true;

	// use the actual device configuration that was negotiated
	unsigned int actualSampleRate = data->device.sampleRate;
	unsigned int actualBufferSize = data->device.playback.internalPeriodSizeInFrames;
	unsigned int actualChannels = data->device.playback.channels;
	// actual format available as data->device.playback.internalFormat

	aSoloud->postinit_internal(actualSampleRate, actualBufferSize, aFlags, actualChannels);

	aSoloud->mBackendCleanupFunc = soloud_miniaudio_deinit;
	aSoloud->mBackendPauseFunc = soloud_miniaudio_pause;
	aSoloud->mBackendResumeFunc = soloud_miniaudio_resume;

	result = ma_device_start(&data->device);
	if (result != MA_SUCCESS)
	{
		soloud_miniaudio_deinit(aSoloud);
		return UNKNOWN_ERROR;
	}

	aSoloud->mBackendString = "MiniAudio";
	return SO_NO_ERROR;
}
}; // namespace SoLoud
#endif
