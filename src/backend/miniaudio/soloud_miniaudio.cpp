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

#include "soloud_internal.h"

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

#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <future>
#include <mutex>
#include <array>

namespace SoLoud
{
using namespace detail; // SAMPLE_FORMAT

struct MiniaudioData
{
	ma_context context{};
	ma_device device{};
	ma_log log{};
	bool contextInitialized{false};
	bool deviceInitialized{false};
	bool logInitialized{false};

	std::atomic<bool> deviceValid{true};
	std::mutex deviceMutex; // protects device operations during shutdown

	struct logCbData {
		ma_uint32 maxLogLevel{0};
		std::mutex logMutex; // logging messages can come from multiple threads
	} logData{};

	Soloud *soloudInstance{nullptr};

	// flag to track if initialization has timed out to prevent background thread interference
	std::atomic<bool> initAbandoned{false};
};

namespace // static
{

void soloud_miniaudio_log_callback(void *pUserData, ma_uint32 level, const char *pMessage)
{
	auto &data = *static_cast<MiniaudioData::logCbData*>(pUserData);
	ma_uint32 maxLevel = data.maxLogLevel;

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

	std::lock_guard<std::mutex> lock(data.logMutex);

	fprintf(stderr, "[MiniAudio %s] %s", levelStr, pMessage); // it seems that the log messages are pre-newline-terminated
}

ma_uint32 parse_log_level_from_env()
{
	const char *env = getenv("SOLOUD_DEBUG");
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

void soloud_miniaudio_notification_callback(const ma_device_notification *pNotification)
{
	auto *data = static_cast<MiniaudioData *>(pNotification->pDevice->pUserData);
	if (!data)
		return;

	switch (pNotification->type)
	{
	case ma_device_notification_type_stopped:
		if (data->logData.maxLogLevel >= MA_LOG_LEVEL_INFO)
			fprintf(stderr, "[MiniAudio INFO] Device stopped\n");
		data->deviceValid.store(false);
		break;

	case ma_device_notification_type_rerouted:
		if (data->logData.maxLogLevel >= MA_LOG_LEVEL_INFO)
			fprintf(stderr, "[MiniAudio INFO] Device rerouted\n");
		// device should continue working after reroute
		break;

	case ma_device_notification_type_interruption_began:
		if (data->logData.maxLogLevel >= MA_LOG_LEVEL_INFO)
			fprintf(stderr, "[MiniAudio INFO] Device interruption began\n");
		data->deviceValid.store(false);
		break;

	case ma_device_notification_type_interruption_ended:
		if (data->logData.maxLogLevel >= MA_LOG_LEVEL_INFO)
			fprintf(stderr, "[MiniAudio INFO] Device interruption ended\n");
		data->deviceValid.store(true);
		break;

	case ma_device_notification_type_started:
		if (data->logData.maxLogLevel >= MA_LOG_LEVEL_INFO)
			fprintf(stderr, "[MiniAudio INFO] Device started\n");
		data->deviceValid.store(true);
		break;

	case ma_device_notification_type_unlocked:
		if (data->logData.maxLogLevel >= MA_LOG_LEVEL_INFO)
			fprintf(stderr, "[MiniAudio INFO] Device unlocked\n");
		break;
	}
}

void soloud_miniaudio_audiomixer(ma_device *pDevice, void *pOutput, const void * /*pInput*/, ma_uint32 frameCount)
{
	auto *data = static_cast<MiniaudioData *>(pDevice->pUserData);
	auto *soloud = data->soloudInstance;

	// check if device is valid
	if (!data->deviceValid.load(std::memory_order_relaxed))
	{
		// output silence when device is invalid
		if (pOutput)
		{
			size_t bytes = static_cast<unsigned long>(frameCount) * pDevice->playback.channels;
			switch (pDevice->playback.internalFormat)
			{
			case ma_format_s16:
				bytes *= sizeof(short);
				memset(pOutput, 0, bytes);
				break;
			case ma_format_u8:
				bytes *= sizeof(unsigned char);
				memset(pOutput, 128, bytes);
				break;
			case ma_format_s24:
				bytes *= 3;
				memset(pOutput, 0, bytes);
				break;
			case ma_format_s32:
				bytes *= sizeof(int);
				memset(pOutput, 0, bytes);
				break;
			case ma_format_f32:
			default:
				bytes *= sizeof(float);
				memset(pOutput, 0, bytes);
				break;
			}
		}
		return;
	}

	// device is valid, proceed with normal mixing
	const ma_format &maFormat = pDevice->playback.internalFormat;
	const SAMPLE_FORMAT outputFormat = maFormat == ma_format_f32   ? SAMPLE_FLOAT32
	                                   : maFormat == ma_format_s16 ? SAMPLE_SIGNED16
	                                   : maFormat == ma_format_u8  ? SAMPLE_UNSIGNED8
	                                   : maFormat == ma_format_s24 ? SAMPLE_SIGNED24
	                                   : maFormat == ma_format_s32 ? SAMPLE_SIGNED32
	                                                               : SAMPLE_FLOAT32;

	soloud->mix(pOutput, frameCount, outputFormat);
}

void soloud_miniaudio_deinit(SoLoud::Soloud *aSoloud)
{
	auto *data = static_cast<MiniaudioData *>(aSoloud->mBackendData);
	if (data)
	{
		// mark as abandoned to prevent any ongoing background operations from interfering
		data->initAbandoned.store(true);
		data->deviceValid.store(false);

		// stop and uninitialize device
		{
			std::lock_guard<std::mutex> lock(data->deviceMutex);
			if (data->deviceInitialized)
			{
				ma_device_stop(&data->device);
				ma_device_uninit(&data->device);
				data->deviceInitialized = false;
			}

			if (data->contextInitialized)
			{
				ma_context_uninit(&data->context);
				data->contextInitialized = false;
			}
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
		std::lock_guard<std::mutex> lock(data->deviceMutex);
		if (data->deviceInitialized)
		{
			if (ma_device_stop(&data->device) != MA_SUCCESS)
				return UNKNOWN_ERROR;
		}
	}
	return SO_NO_ERROR;
}

result soloud_miniaudio_resume(SoLoud::Soloud *aSoloud)
{
	auto *data = static_cast<MiniaudioData *>(aSoloud->mBackendData);
	if (data && data->deviceInitialized)
	{
		std::lock_guard<std::mutex> lock(data->deviceMutex);
		if (data->deviceInitialized)
		{
			if (ma_device_start(&data->device) != MA_SUCCESS)
				return UNKNOWN_ERROR;
		}
	}
	return SO_NO_ERROR;
}

// helper function to attempt device initialization with a specific backend and timeout
ma_result try_backend_with_timeout(MiniaudioData *data, ma_backend backend, const ma_device_config &config)
{
	// create context with specific backend
	ma_context_config contextConfig = ma_context_config_init();
	contextConfig.threadPriority = ma_thread_priority_highest;
	if (data->logInitialized)
		contextConfig.pLog = &data->log;

	ma_result contextResult = ma_context_init(std::array <ma_backend, 1>{backend}.data(), 1, &contextConfig, &data->context);
	if (contextResult != MA_SUCCESS)
		return contextResult;
	data->contextInitialized = true;

	// use timeout for device initialization
	auto deviceInitFuture = std::async(std::launch::deferred, [&data, &config]() -> ma_result {
		ma_result initResult = ma_device_init(&data->context, &config, &data->device);

		// check if we've been abandoned due to timeout
		if (data->initAbandoned.load())
		{
			// we timed out - clean up the device if initialization succeeded
			if (initResult == MA_SUCCESS)
				ma_device_uninit(&data->device);
			return MA_ERROR;
		}

		return initResult;
	});

	// wait for up to 2 seconds for device initialization
	if (deviceInitFuture.wait_for(std::chrono::seconds(6)) == std::future_status::timeout)
	{
		// timeout occurred - mark as abandoned and cleanup
		data->initAbandoned.store(true);

		if (data->logData.maxLogLevel >= MA_LOG_LEVEL_WARNING)
			fprintf(stderr, "[MiniAudio WARNING] Backend '%s' timed out after 2 seconds\n", ma_get_backend_name(backend));

		// cleanup context
		ma_context_uninit(&data->context);
		data->contextInitialized = false;

		return MA_ERROR;
	}

	// get the result from the async operation
	ma_result deviceResult = deviceInitFuture.get();

	// check if we were abandoned during initialization (race condition protection)
	if (data->initAbandoned.load())
	{
		// we were marked as abandoned, cleanup and return error
		if (deviceResult == MA_SUCCESS)
			ma_device_uninit(&data->device);
		ma_context_uninit(&data->context);
		data->contextInitialized = false;
		return MA_ERROR;
	}

	if (deviceResult != MA_SUCCESS)
	{
		// device init failed, cleanup context
		ma_context_uninit(&data->context);
		data->contextInitialized = false;
		return deviceResult;
	}

	// success - device is initialized
	data->deviceInitialized = true;
	return MA_SUCCESS;
}

} // namespace

result miniaudio_init(SoLoud::Soloud *aSoloud, unsigned int aFlags, unsigned int aSamplerate, unsigned int aBuffer, unsigned int aChannels)
{
	auto *data = new MiniaudioData();
	aSoloud->mBackendData = data;
	data->soloudInstance = aSoloud;

	// setup logging
	data->logData.maxLogLevel = parse_log_level_from_env();
	if (data->logData.maxLogLevel > 0)
	{
		ma_result logResult = ma_log_init(nullptr, &data->log);
		if (logResult == MA_SUCCESS)
		{
			data->logInitialized = true;
			ma_log_callback callback = ma_log_callback_init(soloud_miniaudio_log_callback, &data->logData);
			ma_log_register_callback(&data->log, callback);
		}
	}

	// get list of enabled backends
	size_t enabledBackendCount = 0;
	std::array<ma_backend, MA_BACKEND_COUNT> enabledBackends{};
	ma_result result = ma_get_enabled_backends(enabledBackends.data(), MA_BACKEND_COUNT, &enabledBackendCount);
	if (result != MA_SUCCESS)
	{
		if (data->logInitialized)
			ma_log_uninit(&data->log);
		delete data;
		aSoloud->mBackendData = nullptr;
		return UNKNOWN_ERROR;
	}

	// configure device
	ma_device_config config = ma_device_config_init(ma_device_type_playback);
	config.playback.format = ma_format_unknown; // default device format, we'll do the conversion ourselves
	config.playback.channels = aChannels;
	config.dataCallback = soloud_miniaudio_audiomixer;
	config.notificationCallback = soloud_miniaudio_notification_callback;
	config.pUserData = (void *)data;
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
#if defined(_WIN32) || defined(_WIN64) || defined(_MSC_VER) || defined(__linux__)
		config.periodSizeInMilliseconds = 1; // negotiate the lowest possible period size
#endif
	}

	// backend-specific settings
	config.wasapi.noAutoConvertSRC = true; // soloud handles resampling
	config.wasapi.noDefaultQualitySRC = true;
	config.alsa.noAutoFormat = true;
	config.alsa.noAutoChannels = true;
	config.alsa.noAutoResample = true;

	// try each backend individually with timeout
	bool initialized = false;
	for (size_t i = 0; i < enabledBackendCount && !initialized; ++i)
	{
		ma_backend backend = enabledBackends[i];

		// skip null backend since we disabled it
		if (backend == ma_backend_null)
			continue;

		if (data->logData.maxLogLevel >= MA_LOG_LEVEL_INFO)
			fprintf(stderr, "[MiniAudio INFO] Trying backend: %s\n", ma_get_backend_name(backend));

		// reset abandoned flag for this attempt
		data->initAbandoned.store(false);

		result = try_backend_with_timeout(data, backend, config);

		if (result == MA_SUCCESS)
		{
			initialized = true;
			if (data->logData.maxLogLevel >= MA_LOG_LEVEL_INFO)
				fprintf(stderr, "[MiniAudio INFO] Successfully initialized with backend: %s\n", ma_get_backend_name(backend));
		}
		else
		{
			if (data->logData.maxLogLevel >= MA_LOG_LEVEL_WARNING)
			{
				if (result == MA_ERROR) // this indicates timeout
					fprintf(stderr, "[MiniAudio WARNING] Backend '%s' failed due to timeout\n", ma_get_backend_name(backend));
				else
					fprintf(stderr, "[MiniAudio WARNING] Backend '%s' failed with error: %d\n", ma_get_backend_name(backend), result);
			}
		}
	}

	if (!initialized)
	{
		if (data->logInitialized)
			ma_log_uninit(&data->log);
		delete data;
		aSoloud->mBackendData = nullptr;

		if (data->logData.maxLogLevel >= MA_LOG_LEVEL_ERROR)
			fprintf(stderr, "[MiniAudio ERROR] All backends failed to initialize\n");

		return UNKNOWN_ERROR; // this will cause soloud to try other backends like SDL3
	}

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
