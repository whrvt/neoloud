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

#include "soloud_config.h"
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

#include <algorithm>
#include <array>
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <vector>

namespace SoLoud
{
using namespace detail; // SAMPLE_FORMAT

struct MiniaudioData
{
	MiniaudioData() { std::memset(&this->device, 0, sizeof(ma_device)); }
	ma_context context{};
	ma_device device;
	ma_log log{};
	ma_device_info currentDeviceInfo{};

	std::mutex deviceMutex; // protects device operations during shutdown

	// parent instance
	Soloud *soloudInstance{nullptr};

	// device management
	ma_backend currentBackend{ma_backend_null};

	// cached device configuration for seamless switching
	unsigned int initFlags{0};
	unsigned int requestedSampleRate{0};
	unsigned int requestedBufferSize{0};
	unsigned int requestedChannels{0};
	ma_uint32 maxLogLevel{0};

	// flags
	std::atomic<bool> deviceValid{true};
	std::atomic<bool> soloudInitialized{false}; // postinit_internal must be called before we use soloud->mix in the callback
	bool logInitialized{false};
	bool contextInitialized{false};
	bool deviceInitialized{false};
	bool hasCurrentDeviceInfo{false};
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

	const bool printNewline = (pMessage[strlen(pMessage) - 1] != '\n');
	SoLoud::logStdout("[MiniAudio %s] %s%s", levelStr, pMessage, printNewline ? "\n" : "");
}

// check if device supports exclusive mode by examining native data formats
bool device_supports_exclusive_mode(const ma_device_info *pDeviceInfo)
{
	if (!pDeviceInfo)
		return false;

	for (ma_uint32 i = 0; i < pDeviceInfo->nativeDataFormatCount; i++)
	{
		if (pDeviceInfo->nativeDataFormats[i].flags & MA_DATA_FORMAT_FLAG_EXCLUSIVE_MODE)
			return true;
	}
	return false;
}

// parse share mode from device identifier
ma_share_mode parse_share_mode_from_identifier(const char *identifier)
{
	if (!identifier)
		return ma_share_mode_shared;

	size_t len = strlen(identifier);
	if (len >= 2)
	{
		// check for "_e" suffix (exclusive mode)
		if (identifier[len - 2] == '_' && identifier[len - 1] == 'e')
			return ma_share_mode_exclusive;
		// check for "_s" suffix (shared mode)
		if (identifier[len - 2] == '_' && identifier[len - 1] == 's')
			return ma_share_mode_shared;
	}

	return ma_share_mode_shared; // default to shared
}

// remove share mode suffix from identifier to get base identifier for device matching
void get_base_identifier(const char *identifier, char *baseIdentifier, size_t baseIdentifierSize)
{
	if (!identifier || !baseIdentifier || baseIdentifierSize == 0)
		return;

	size_t len = strlen(identifier);

	// remove "_e" or "_s" suffix if present
	if (len >= 2 && identifier[len - 2] == '_' && (identifier[len - 1] == 'e' || identifier[len - 1] == 's'))
		len -= 2;

	if (len >= baseIdentifierSize)
		len = baseIdentifierSize - 1;

	memcpy(baseIdentifier, identifier, len);
	baseIdentifier[len] = '\0';
}

#ifdef _MSC_VER
#include <string.h>
#ifndef strncasecmp
#define strncasecmp _strnicmp
#endif
#endif

ma_uint32 parse_log_level_from_env()
{
	const char *env = getenv("SOLOUD_DEBUG");
	if (!env || !*env || *env == '0')
	{
#ifdef _DEBUG
		return MA_LOG_LEVEL_WARNING;
#else
		return 0; // none
#endif
	}

	if (strncasecmp(env, "debug", sizeof("debug") - 1) == 0)
		return MA_LOG_LEVEL_DEBUG;
	if (strncasecmp(env, "info", sizeof("info") - 1) == 0)
		return MA_LOG_LEVEL_INFO;
	if (strncasecmp(env, "warn", sizeof("warn") - 1) == 0)
		return MA_LOG_LEVEL_WARNING;
	if (strncasecmp(env, "error", sizeof("error") - 1) == 0)
		return MA_LOG_LEVEL_ERROR;
	if (strncasecmp(env, "none", sizeof("none") - 1) == 0)
		return 0;

	// default fallback for unrecognized values
#ifdef _DEBUG
	return MA_LOG_LEVEL_WARNING;
#else
	return 0;
#endif
}

// use an environment variable that's mostly compatible with SDL_AUDIO_DRIVER, so we don't have
// to do anything extra for the SDL output backend
ma_backend parse_backend_from_env()
{
	const char *env = getenv("SOLOUD_MINIAUDIO_DRIVER");
	if (!env || !*env || *env == '0')
	{
		return ma_backend_null; // none
	}

	if (strncasecmp(env, "wasapi", sizeof("wasapi") - 1) == 0)
		return ma_backend_wasapi;
	if (strncasecmp(env, "dsound", sizeof("dsound") - 1) == 0)
		return ma_backend_dsound;
	if (strncasecmp(env, "winmm", sizeof("winmm") - 1) == 0)
		return ma_backend_winmm;
	if (strncasecmp(env, "coreaudio", sizeof("coreaudio") - 1) == 0)
		return ma_backend_coreaudio;
	if (strncasecmp(env, "sndio", sizeof("sndio") - 1) == 0)
		return ma_backend_sndio;
	if (strncasecmp(env, "audio4", sizeof("audio4") - 1) == 0)
		return ma_backend_audio4;
	if (strncasecmp(env, "oss", sizeof("oss") - 1) == 0)
		return ma_backend_oss;
	if (strncasecmp(env, "pulseaudio", sizeof("pulseaudio") - 1) == 0)
		return ma_backend_pulseaudio;
	if (strncasecmp(env, "alsa", sizeof("alsa") - 1) == 0)
		return ma_backend_alsa;
	if (strncasecmp(env, "jack", sizeof("jack") - 1) == 0)
		return ma_backend_jack;
	if (strncasecmp(env, "aaudio", sizeof("aaudio") - 1) == 0)
		return ma_backend_aaudio;
	if (strncasecmp(env, "opensl", sizeof("opensl") - 1) == 0)
		return ma_backend_opensl;
	if (strncasecmp(env, "webaudio", sizeof("webaudio") - 1) == 0)
		return ma_backend_webaudio;

	// default fallback
	return ma_backend_null;
}

void soloud_miniaudio_notification_callback(const ma_device_notification *pNotification)
{
	auto *data = static_cast<MiniaudioData *>(pNotification->pDevice->pUserData);
	if (!data)
		return;

	const char *logStr = nullptr;
	switch (pNotification->type)
	{
	case ma_device_notification_type_stopped:
		logStr = "[MiniAudio INFO] Device stopped\n";
		data->deviceValid.store(false);
		break;

	case ma_device_notification_type_rerouted:
		logStr = "[MiniAudio INFO] Device rerouted\n";
		// device should continue working after reroute
		break;

	case ma_device_notification_type_interruption_began:
		logStr = "[MiniAudio INFO] Device interruption began\n";
		data->deviceValid.store(false);
		break;

	case ma_device_notification_type_interruption_ended:
		logStr = "[MiniAudio INFO] Device interruption ended\n";
		data->deviceValid.store(true);
		break;

	case ma_device_notification_type_started:
		logStr = "[MiniAudio INFO] Device started\n";
		data->deviceValid.store(true);
		break;

	case ma_device_notification_type_unlocked:
		logStr = "[MiniAudio INFO] Device unlocked\n";
		break;
	}

	if (logStr && data->maxLogLevel >= MA_LOG_LEVEL_INFO)
		SoLoud::logStdout("%s", logStr);
}

void soloud_miniaudio_audiomixer(ma_device *pDevice, void *pOutput, const void * /*pInput*/, ma_uint32 frameCount)
{
	auto *data = static_cast<MiniaudioData *>(pDevice->pUserData);
	auto *soloud = data->soloudInstance;

	// check if device is valid
	if (!(data->soloudInitialized.load(std::memory_order_acquire) && data->deviceValid.load(std::memory_order_relaxed)))
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

void soloud_miniaudio_deinit(Soloud *aSoloud)
{
	auto *data = static_cast<MiniaudioData *>(aSoloud->mBackendData);
	if (data)
	{
		// mark as invalid to prevent any ongoing background operations from interfering
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

result soloud_miniaudio_pause(Soloud *aSoloud)
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

result soloud_miniaudio_resume(Soloud *aSoloud)
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
ma_device_config create_device_config(MiniaudioData *data, ma_share_mode shareMode, const ma_device_id *pDeviceID = nullptr)
{
	ma_device_config config = ma_device_config_init(ma_device_type_playback);
	config.playback.format = ma_format_unknown; // default device format, we'll do the conversion ourselves
	config.playback.channels = data->requestedChannels;
	config.playback.shareMode = shareMode;
	config.dataCallback = soloud_miniaudio_audiomixer;
	config.notificationCallback = soloud_miniaudio_notification_callback;
	config.pUserData = (void *)data;
	config.noPreSilencedOutputBuffer = true;
	config.noFixedSizedCallback = true;
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
	config.wasapi.usage = ma_wasapi_usage_pro_audio; // sets higher priority through MMCSS
	config.alsa.noAutoFormat = true;
	config.alsa.noAutoChannels = true;
	config.alsa.noAutoResample = true;

	return config;
}

// helper function to attempt device initialization with a specific backend and timeout
ma_result try_backend_with_timeout(MiniaudioData *data, ma_backend backend, ma_share_mode shareMode)
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

	// create device config with specified share mode
	ma_device_config config = create_device_config(data, shareMode);

	// log start of device initialization for debugging potential hangs
	if (data->maxLogLevel >= MA_LOG_LEVEL_INFO)
		SoLoud::logStdout("[MiniAudio INFO] Initializing device with backend '%s' in %s mode...\n", ma_get_backend_name(backend),
		                  shareMode == ma_share_mode_exclusive ? "exclusive" : "shared");

	ma_result deviceResult = ma_device_init(&data->context, &config, &data->device);

	if (data->maxLogLevel >= MA_LOG_LEVEL_INFO)
	{
		if (deviceResult == MA_SUCCESS)
			SoLoud::logStdout("[MiniAudio INFO] Device initialization completed successfully\n");
		else
			SoLoud::logStdout("[MiniAudio INFO] Device initialization failed with result: %d\n", deviceResult);
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
ma_result init_device_with_id(MiniaudioData *data, ma_share_mode shareMode, const ma_device_id *pDeviceID)
{
	if (!data->contextInitialized)
		return MA_ERROR;

	ma_device_config config = create_device_config(data, shareMode, pDeviceID);

	// log start of device initialization
	if (data->maxLogLevel >= MA_LOG_LEVEL_INFO)
		SoLoud::logStdout("[MiniAudio INFO] Initializing device with specific ID in %s mode...\n", shareMode == ma_share_mode_exclusive ? "exclusive" : "shared");

	ma_result deviceResult = ma_device_init(&data->context, &config, &data->device);

	if (data->maxLogLevel >= MA_LOG_LEVEL_INFO)
	{
		if (deviceResult == MA_SUCCESS)
			SoLoud::logStdout("[MiniAudio INFO] Device initialization completed successfully\n");
		else
			SoLoud::logStdout("[MiniAudio INFO] Device initialization failed with result: %d\n", deviceResult);
	}

	if (deviceResult == MA_SUCCESS)
		data->deviceInitialized = true;

	return deviceResult;
}

// convert ma_device_info to generic DeviceInfo
void convert_device_info(const ma_device_info *pMaInfo, DeviceInfo *pGenericInfo, ma_share_mode shareMode, ma_backend backend)
{
	if (!pMaInfo || !pGenericInfo)
		return;

	// copy name with bounds checking
	size_t nameLen = strlen(&pMaInfo->name[0]);
	if (nameLen >= sizeof(pGenericInfo->name))
		nameLen = sizeof(pGenericInfo->name) - 1;
	memcpy(pGenericInfo->name.data(), &pMaInfo->name[0], nameLen);
	pGenericInfo->name[nameLen] = '\0';

	// append share mode suffix to name for clarity
	if (backend == ma_backend_wasapi)
	{
		const char *modeSuffix = (shareMode == ma_share_mode_exclusive) ? " (Exclusive)" : " (Shared)";
		size_t suffixLen = strlen(modeSuffix);
		if (nameLen + suffixLen < sizeof(pGenericInfo->name))
		{
			memcpy(pGenericInfo->name.data() + nameLen, modeSuffix, suffixLen);
			pGenericInfo->name[nameLen + suffixLen] = '\0';
		}
	}

	// create string representation of device ID with share mode encoded
	const char *modeTag = (shareMode == ma_share_mode_exclusive) ? "_e" : "_s";
	snprintf(pGenericInfo->identifier.data(), pGenericInfo->identifier.size(), "ma_%d_%d_%d_%s%s", (int)pMaInfo->id.wasapi[0], (int)pMaInfo->id.wasapi[1],
	         (int)pMaInfo->id.wasapi[2], &pMaInfo->name[0], modeTag);

	pGenericInfo->isDefault = (pMaInfo->isDefault != 0 && shareMode == ma_share_mode_shared); // only create 1 default device (we init in shared mode)
	pGenericInfo->isExclusive = (shareMode == ma_share_mode_exclusive);
	pGenericInfo->nativeDeviceInfo = nullptr;
}

// enumerate available playback devices
result miniaudio_enumerate_devices(Soloud *aSoloud)
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

	// build device list
	std::vector<DeviceInfo> devices;
	devices.reserve(maDeviceCount * (1ULL + (data->currentBackend == ma_backend_wasapi))); // reserve for worst case (shared + exclusive)
	bool defaultSupportsExclusive = false;

	for (ma_uint32 i = 0; i < maDeviceCount; i++)
	{
		ma_device_info detailedInfo;
		ma_result infoResult = ma_context_get_device_info(&data->context, ma_device_type_playback, &maDevices[i].id, &detailedInfo);

		bool hasDetailedInfo = (infoResult == MA_SUCCESS);
		bool supportsExclusive = false;

		if (hasDetailedInfo)
		{
			supportsExclusive = device_supports_exclusive_mode(&detailedInfo);

			// track default device's exclusive mode support for fallback logic
			if (maDevices[i].isDefault)
				defaultSupportsExclusive = supportsExclusive;
		}
		else if (infoResult == MA_BUSY && data->currentBackend == ma_backend_wasapi)
		{
			// device is busy (likely in use), assume same capabilities as default device
			supportsExclusive = defaultSupportsExclusive;
		}

		// use detailed info if available, otherwise fall back to basic info
		const ma_device_info *pInfoToUse = hasDetailedInfo ? &detailedInfo : &maDevices[i];

		// always add shared mode entry
		DeviceInfo sharedDevice;
		convert_device_info(pInfoToUse, &sharedDevice, ma_share_mode_shared, data->currentBackend);
		devices.push_back(sharedDevice);

		// add exclusive mode entry if supported
		if (supportsExclusive)
		{
			DeviceInfo exclusiveDevice;
			convert_device_info(pInfoToUse, &exclusiveDevice, ma_share_mode_exclusive, data->currentBackend);
			devices.push_back(exclusiveDevice);
		}
	}

	// move info to soloud device list
	aSoloud->mDeviceList = new DeviceInfo[devices.size()];
	std::move(devices.begin(), devices.end(), aSoloud->mDeviceList);
	aSoloud->mDeviceCount = static_cast<unsigned int>(devices.size());

	return SO_NO_ERROR;
}

// get current device information
result miniaudio_get_current_device_info(Soloud *aSoloud, DeviceInfo *pDeviceInfo)
{
	if (!aSoloud || !pDeviceInfo)
		return INVALID_PARAMETER;

	auto *data = static_cast<MiniaudioData *>(aSoloud->mBackendData);
	if (!data || !data->hasCurrentDeviceInfo)
		return INVALID_PARAMETER;

	// determine share mode from device config
	ma_share_mode currentShareMode = data->deviceInitialized ? data->device.playback.shareMode : ma_share_mode_shared;

	convert_device_info(&data->currentDeviceInfo, pDeviceInfo, currentShareMode, data->currentBackend);
	return SO_NO_ERROR;
}

// switch to a different playback device
result miniaudio_set_device(Soloud *aSoloud, const char *deviceIdentifier)
{
	if (!aSoloud)
		return INVALID_PARAMETER;

	auto *data = static_cast<MiniaudioData *>(aSoloud->mBackendData);
	if (!data || !data->contextInitialized)
		return INVALID_PARAMETER;

	ma_device_id *targetDeviceId = nullptr;
	ma_share_mode targetShareMode = ma_share_mode_shared;

	// handle default device case
	if (!deviceIdentifier || strlen(deviceIdentifier) == 0)
	{
		targetDeviceId = nullptr; // use default device
		targetShareMode = ma_share_mode_shared;
	}
	else
	{
		// parse share mode from identifier
		targetShareMode = parse_share_mode_from_identifier(deviceIdentifier);

		// get base identifier for device matching
		std::array<char, 256> baseIdentifier{};
		get_base_identifier(deviceIdentifier, baseIdentifier.data(), baseIdentifier.size());

		// enumerate devices to find matching identifier
		ma_device_info *maDevices = nullptr;
		ma_uint32 maDeviceCount = 0;
		ma_result enumResult = ma_context_get_devices(&data->context, &maDevices, &maDeviceCount, nullptr, nullptr);
		if (enumResult != MA_SUCCESS)
			return UNKNOWN_ERROR;

		// find device with matching base identifier
		for (ma_uint32 i = 0; i < maDeviceCount; i++)
		{
			DeviceInfo tempInfo{};
			convert_device_info(&maDevices[i], &tempInfo, ma_share_mode_shared, data->currentBackend);

			// get base identifier from temp info for comparison
			std::array<char, 256> tempBaseIdentifier{};
			get_base_identifier(tempInfo.identifier.data(), tempBaseIdentifier.data(), tempBaseIdentifier.size());

			if (tempBaseIdentifier == baseIdentifier)
			{
				targetDeviceId = &maDevices[i].id;
				break;
			}
		}

		if (!targetDeviceId)
			return INVALID_PARAMETER;
	}

	// lock device mutex
	// miniaudio will internally block until the data callback is done, we don't have to lock the soloud internal audio mutex
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

	// initialize new device with specified share mode
	ma_result result = init_device_with_id(data, targetShareMode, targetDeviceId);
	if (result != MA_SUCCESS)
	{
		if (data->maxLogLevel >= MA_LOG_LEVEL_ERROR)
			SoLoud::logStdout("[MiniAudio ERROR] Failed to initialize new device\n");

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
		if (data->maxLogLevel >= MA_LOG_LEVEL_INFO)
		{
			SoLoud::logStdout("[MiniAudio INFO] Device configuration changed: %uHz/%uch/%uf -> %uHz/%uch/%uf\n", oldSampleRate, oldChannels, oldBufferSize,
			                  newSampleRate, newChannels, newBufferSize);
		}
		// update SoLoud's internal configuration
		aSoloud->postinit_internal(newSampleRate, newBufferSize, data->initFlags, newChannels);
	}

	// start the new device
	result = ma_device_start(&data->device);
	if (result != MA_SUCCESS)
	{
		if (data->maxLogLevel >= MA_LOG_LEVEL_ERROR)
			SoLoud::logStdout("[MiniAudio ERROR] Failed to start new device\n");

		return UNKNOWN_ERROR;
	}

	data->deviceValid.store(true);

	if (data->maxLogLevel >= MA_LOG_LEVEL_INFO)
		SoLoud::logStdout("[MiniAudio INFO] Successfully switched to new device in %s mode\n", targetShareMode == ma_share_mode_exclusive ? "exclusive" : "shared");

	return SO_NO_ERROR;
}

} // namespace

result miniaudio_init(Soloud *aSoloud, unsigned int aFlags, unsigned int aSamplerate, unsigned int aBuffer, unsigned int aChannels)
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

	// get list of enabled backends
	size_t enabledBackendCount = 0;
	std::array<ma_backend, MA_BACKEND_COUNT + 1> enabledBackends{}; // allocate 1 more slot to hold an environment variable choice at the start
	ma_result result = ma_get_enabled_backends(enabledBackends.data() + 1, MA_BACKEND_COUNT, &enabledBackendCount);
	if (result != MA_SUCCESS || enabledBackendCount == 0)
	{
		if (data->logInitialized)
			ma_log_uninit(&data->log);
		delete data;
		aSoloud->mBackendData = nullptr;
		return UNKNOWN_ERROR;
	}

	size_t beginIndex = 1;
	const ma_backend user_backend = parse_backend_from_env();

	if (user_backend != ma_backend_null)
	{
		auto it = std::find(enabledBackends.begin(), enabledBackends.end(), user_backend);
		if (it != enabledBackends.end())
		{
			// reorder the user choice to be at the start, and set the original element to null (so we skip it later)
			beginIndex = 0;
			enabledBackends[0] = user_backend;
			*it = ma_backend_null;
		}
	}

	// try each backend individually with shared mode as default
	bool initialized = false;
	for (size_t i = beginIndex; i < (enabledBackendCount + 1 /* because we put the real elements 1 past the start */) && !initialized; ++i)
	{
		ma_backend backend = enabledBackends[i];

		// skip null backend since we disabled it
		if (backend == ma_backend_null)
			continue;

		ma_share_mode initShareMode = (backend == ma_backend_wasapi) && (data->initFlags & Soloud::INIT_EXCLUSIVE) ? ma_share_mode_exclusive : ma_share_mode_shared;
		for (size_t i = 0; i < 2; i++) // try shared if exclusive failed
		{
			const char *modeString = initShareMode == ma_share_mode_exclusive ? "exclusive mode" : "shared mode";

			if (data->maxLogLevel >= MA_LOG_LEVEL_INFO)
				SoLoud::logStdout("[MiniAudio INFO] Trying backend: %s in %s\n", ma_get_backend_name(backend), modeString);

			result = try_backend_with_timeout(data, backend, initShareMode);

			if (result == MA_SUCCESS)
			{
				initialized = true;
				if (data->maxLogLevel >= MA_LOG_LEVEL_INFO)
					SoLoud::logStdout("[MiniAudio INFO] Successfully initialized with backend: %s\n", ma_get_backend_name(backend));
				break;
			}
			else
			{
				if (data->maxLogLevel >= MA_LOG_LEVEL_WARNING)
				{
					if (result == MA_ERROR) // this indicates timeout
						SoLoud::logStdout("[MiniAudio WARNING] Backend '%s' failed in %s due to timeout\n", ma_get_backend_name(backend), modeString);
					else
						SoLoud::logStdout("[MiniAudio WARNING] Backend '%s' failed in %s with error: %d\n", ma_get_backend_name(backend), modeString, result);
				}

				if (initShareMode != ma_share_mode_exclusive)
				{
					// we failed in shared mode, move onto the next backend
					break;
				}
				else // try again in shared mode
				{
					initShareMode = ma_share_mode_shared;
				}
			}
		}
	}

	if (!initialized)
	{
		if (data->logInitialized)
			ma_log_uninit(&data->log);
		delete data;
		aSoloud->mBackendData = nullptr;

		if (data->maxLogLevel >= MA_LOG_LEVEL_ERROR)
			SoLoud::logStdout("[MiniAudio ERROR] All backends failed to initialize\n");

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
	data->soloudInitialized.store(true);

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
