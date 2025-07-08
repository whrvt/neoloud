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

#include <array>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <future>
#include <mutex>

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

	struct logCbData
	{
		ma_uint32 maxLogLevel{0};
		std::mutex logMutex; // logging messages can come from multiple threads
	} logData{};

	Soloud *soloudInstance{nullptr};

	// flag to track if initialization has timed out to prevent background thread interference
	std::atomic<bool> initAbandoned{false};

	// device management
	ma_device_info currentDeviceInfo{};
	bool hasCurrentDeviceInfo{false};
	ma_backend currentBackend{ma_backend_null};

	// cached device configuration for seamless switching
	unsigned int initFlags{0};
	unsigned int requestedSampleRate{0};
	unsigned int requestedBufferSize{0};
	unsigned int requestedChannels{0};
};

namespace // static
{

void soloud_miniaudio_log_callback(void *pUserData, ma_uint32 level, const char *pMessage)
{
	auto &data = *static_cast<MiniaudioData::logCbData *>(pUserData);
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

// helper function to create device config from cached parameters
ma_device_config create_device_config(MiniaudioData *data, const ma_device_id *pDeviceID = nullptr)
{
	ma_device_config config = ma_device_config_init(ma_device_type_playback);
	config.playback.format = ma_format_unknown; // default device format, we'll do the conversion ourselves
	config.playback.channels = data->requestedChannels;
	config.dataCallback = soloud_miniaudio_audiomixer;
	config.notificationCallback = soloud_miniaudio_notification_callback;
	config.pUserData = (void *)data;
	config.noPreSilencedOutputBuffer = true;
	config.noClip = true;
	config.performanceProfile = ma_performance_profile_low_latency;

	if (pDeviceID)
		config.playback.pDeviceID = pDeviceID;

	if (data->requestedSampleRate > 0)
		config.sampleRate = data->requestedSampleRate;
	else
		config.sampleRate = 0;

	if (data->requestedBufferSize > 0)
		config.periodSizeInFrames = data->requestedBufferSize;
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

	return config;
}

// helper function to attempt device initialization with a specific backend and timeout
ma_result try_backend_with_timeout(MiniaudioData *data, ma_backend backend, const ma_device_config &config)
{
	// create context with specific backend
	ma_context_config contextConfig = ma_context_config_init();
	contextConfig.threadPriority = ma_thread_priority_highest;
	if (data->logInitialized)
		contextConfig.pLog = &data->log;

	ma_result contextResult = ma_context_init(std::array<ma_backend, 1>{backend}.data(), 1, &contextConfig, &data->context);
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

	// wait for up to 6 seconds for device initialization
	if (deviceInitFuture.wait_for(std::chrono::seconds(6)) == std::future_status::timeout)
	{
		// timeout occurred - mark as abandoned and cleanup
		data->initAbandoned.store(true);

		if (data->logData.maxLogLevel >= MA_LOG_LEVEL_WARNING)
			fprintf(stderr, "[MiniAudio WARNING] Backend '%s' timed out after 6 seconds\n", ma_get_backend_name(backend));

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
	data->currentBackend = backend;
	return MA_SUCCESS;
}

// helper function to initialize device with specific device ID (used during device switching)
ma_result init_device_with_id(MiniaudioData *data, const ma_device_id *pDeviceID)
{
	if (!data->contextInitialized)
		return MA_ERROR;

	ma_device_config config = create_device_config(data, pDeviceID);

	// reset abandoned flag for this attempt
	data->initAbandoned.store(false);

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

	// wait for up to 6 seconds for device initialization
	if (deviceInitFuture.wait_for(std::chrono::seconds(6)) == std::future_status::timeout)
	{
		// timeout occurred - mark as abandoned and cleanup
		data->initAbandoned.store(true);

		if (data->logData.maxLogLevel >= MA_LOG_LEVEL_WARNING)
			fprintf(stderr, "[MiniAudio WARNING] Device initialization timed out after 6 seconds\n");

		return MA_ERROR;
	}

	// get the result from the async operation
	ma_result deviceResult = deviceInitFuture.get();

	// check if we were abandoned during initialization
	if (data->initAbandoned.load())
	{
		// we were marked as abandoned, cleanup and return error
		if (deviceResult == MA_SUCCESS)
			ma_device_uninit(&data->device);
		return MA_ERROR;
	}

	if (deviceResult == MA_SUCCESS)
		data->deviceInitialized = true;

	return deviceResult;
}

// convert ma_device_info to generic DeviceInfo
void convert_device_info(const ma_device_info *pMaInfo, DeviceInfo *pGenericInfo)
{
	if (!pMaInfo || !pGenericInfo)
		return;

	// copy name with bounds checking
	size_t nameLen = strlen(&pMaInfo->name[0]);
	if (nameLen >= sizeof(pGenericInfo->name))
		nameLen = sizeof(pGenericInfo->name) - 1;
	memcpy(&pGenericInfo->name[0], &pMaInfo->name[0], nameLen);
	pGenericInfo->name[nameLen] = '\0';

	// create string representation of device ID
	snprintf(&pGenericInfo->identifier[0], sizeof(pGenericInfo->identifier), "ma_%d_%d_%d_%s",
	         (int)pMaInfo->id.wasapi[0], // use first few bytes as identifier
	         (int)pMaInfo->id.wasapi[1], (int)pMaInfo->id.wasapi[2], &pMaInfo->name[0]);

	pGenericInfo->isDefault = pMaInfo->isDefault != 0;
	pGenericInfo->nativeDeviceInfo = nullptr;
}

// enumerate available playback devices
result miniaudio_enumerate_devices(SoLoud::Soloud *aSoloud)
{
	if (!aSoloud)
		return INVALID_PARAMETER;

	auto *data = static_cast<MiniaudioData *>(aSoloud->mBackendData);
	if (!data || !data->contextInitialized)
		return INVALID_PARAMETER;

	std::lock_guard<std::mutex> lock(data->deviceMutex);

	ma_device_info *maDevices = nullptr;
	ma_uint32 maDeviceCount = 0;

	ma_result result = ma_context_get_devices(&data->context, &maDevices, &maDeviceCount, nullptr, nullptr);
	if (result != MA_SUCCESS)
		return UNKNOWN_ERROR;

	// allocate internal device array
	aSoloud->mDeviceList = new DeviceInfo[maDeviceCount];
	if (!aSoloud->mDeviceList)
	{
		aSoloud->mDeviceCount = 0;
		return OUT_OF_MEMORY;
	}

	// convert ma_device_info to DeviceInfo
	for (ma_uint32 i = 0; i < maDeviceCount; i++)
	{
		convert_device_info(&maDevices[i], &aSoloud->mDeviceList[i]);
	}

	aSoloud->mDeviceCount = maDeviceCount;
	return SO_NO_ERROR;
}

// get current device information
result miniaudio_get_current_device_info(SoLoud::Soloud *aSoloud, DeviceInfo *pDeviceInfo)
{
	if (!aSoloud || !pDeviceInfo)
		return INVALID_PARAMETER;

	auto *data = static_cast<MiniaudioData *>(aSoloud->mBackendData);
	if (!data || !data->hasCurrentDeviceInfo)
		return INVALID_PARAMETER;

	convert_device_info(&data->currentDeviceInfo, pDeviceInfo);
	return SO_NO_ERROR;
}

// switch to a different playback device
result miniaudio_set_device(SoLoud::Soloud *aSoloud, const char *deviceIdentifier)
{
	if (!aSoloud)
		return INVALID_PARAMETER;

	auto *data = static_cast<MiniaudioData *>(aSoloud->mBackendData);
	if (!data || !data->contextInitialized)
		return INVALID_PARAMETER;

	ma_device_id *targetDeviceId = nullptr;

	// handle default device case
	if (!deviceIdentifier || strlen(deviceIdentifier) == 0)
	{
		targetDeviceId = nullptr; // use default device
	}
	else
	{
		// enumerate devices to find matching identifier
		ma_device_info *maDevices = nullptr;
		ma_uint32 maDeviceCount = 0;
		ma_result enumResult = ma_context_get_devices(&data->context, &maDevices, &maDeviceCount, nullptr, nullptr);
		if (enumResult != MA_SUCCESS)
			return UNKNOWN_ERROR;

		// find device with matching identifier
		for (ma_uint32 i = 0; i < maDeviceCount; i++)
		{
			DeviceInfo tempInfo{};
			convert_device_info(&maDevices[i], &tempInfo);
			if (strcmp(&tempInfo.identifier[0], deviceIdentifier) == 0)
			{
				targetDeviceId = &maDevices[i].id;
				break;
			}
		}

		if (!targetDeviceId)
			return INVALID_PARAMETER;
	}

	// lock both audio and device mutexes to ensure safe switching
	aSoloud->lockAudioMutex_internal();
	std::lock_guard<std::mutex> deviceLock(data->deviceMutex);

	// store current configuration
	unsigned int oldSampleRate = data->deviceInitialized ? data->device.sampleRate : 0;
	unsigned int oldBufferSize = data->deviceInitialized ? data->device.playback.internalPeriodSizeInFrames : 0;
	unsigned int oldChannels = data->deviceInitialized ? data->device.playback.channels : 0;

	// stop and uninitialize current device
	if (data->deviceInitialized)
	{
		ma_device_stop(&data->device);
		ma_device_uninit(&data->device);
		data->deviceInitialized = false;
		data->deviceValid.store(false);
	}

	// initialize new device
	ma_result result = init_device_with_id(data, targetDeviceId);
	if (result != MA_SUCCESS)
	{
		aSoloud->unlockAudioMutex_internal();

		if (data->logData.maxLogLevel >= MA_LOG_LEVEL_ERROR)
			fprintf(stderr, "[MiniAudio ERROR] Failed to initialize new device\n");

		return UNKNOWN_ERROR;
	}

	// update current device info
	ma_result deviceInfoResult = ma_context_get_device_info(&data->context, ma_device_type_playback, data->device.playback.pID, &data->currentDeviceInfo);
	data->hasCurrentDeviceInfo = (deviceInfoResult == MA_SUCCESS);

	// check if device configuration changed
	unsigned int newSampleRate = data->device.sampleRate;
	unsigned int newBufferSize = data->device.playback.internalPeriodSizeInFrames;
	unsigned int newChannels = data->device.playback.channels;

	bool configChanged = (oldSampleRate != newSampleRate) || (oldBufferSize != newBufferSize) || (oldChannels != newChannels);

	if (configChanged)
	{
		if (data->logData.maxLogLevel >= MA_LOG_LEVEL_INFO)
		{
			fprintf(stderr, "[MiniAudio INFO] Device configuration changed: %uHz/%uch/%uf -> %uHz/%uch/%uf\n", oldSampleRate, oldChannels, oldBufferSize,
			        newSampleRate, newChannels, newBufferSize);
		}

		// unlock audio mutex before calling postinit_internal as it may need to reinitialize internal state
		aSoloud->unlockAudioMutex_internal();

		// update SoLoud's internal configuration
		aSoloud->postinit_internal(newSampleRate, newBufferSize, data->initFlags, newChannels);

		// re-lock for device start
		aSoloud->lockAudioMutex_internal();
	}

	// start the new device
	result = ma_device_start(&data->device);
	if (result != MA_SUCCESS)
	{
		aSoloud->unlockAudioMutex_internal();

		if (data->logData.maxLogLevel >= MA_LOG_LEVEL_ERROR)
			fprintf(stderr, "[MiniAudio ERROR] Failed to start new device\n");

		return UNKNOWN_ERROR;
	}

	data->deviceValid.store(true);
	aSoloud->unlockAudioMutex_internal();

	if (data->logData.maxLogLevel >= MA_LOG_LEVEL_INFO)
		fprintf(stderr, "[MiniAudio INFO] Successfully switched to new device\n");

	return SO_NO_ERROR;
}

} // namespace

result miniaudio_init(SoLoud::Soloud *aSoloud, unsigned int aFlags, unsigned int aSamplerate, unsigned int aBuffer, unsigned int aChannels)
{
	auto *data = new MiniaudioData();
	aSoloud->mBackendData = data;
	data->soloudInstance = aSoloud;

	// cache initialization parameters for device switching
	data->initFlags = aFlags;
	data->requestedSampleRate = aSamplerate;
	data->requestedBufferSize = aBuffer;
	data->requestedChannels = aChannels;

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
	ma_device_config config = create_device_config(data);

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

	// store current device info
	ma_result deviceInfoResult = ma_context_get_device_info(&data->context, ma_device_type_playback, data->device.playback.pID, &data->currentDeviceInfo);
	if (deviceInfoResult == MA_SUCCESS)
		data->hasCurrentDeviceInfo = true;

	// use the actual device configuration that was negotiated
	unsigned int actualSampleRate = data->device.sampleRate;
	unsigned int actualBufferSize = data->device.playback.internalPeriodSizeInFrames;
	unsigned int actualChannels = data->device.playback.channels;

	aSoloud->postinit_internal(actualSampleRate, actualBufferSize, aFlags, actualChannels);

	aSoloud->mBackendCleanupFunc = soloud_miniaudio_deinit;
	aSoloud->mBackendPauseFunc = soloud_miniaudio_pause;
	aSoloud->mBackendResumeFunc = soloud_miniaudio_resume;

	// setup device management function pointers
	aSoloud->mEnumerateDevicesFunc = miniaudio_enumerate_devices;
	aSoloud->mGetCurrentDeviceFunc = miniaudio_get_current_device_info;
	aSoloud->mSetDeviceFunc = miniaudio_set_device;

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
