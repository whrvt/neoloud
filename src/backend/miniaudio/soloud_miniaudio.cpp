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
#define MA_NO_JACK // we have pulse and alsa already, this just bloats the device struct
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
#include <condition_variable>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <thread>

namespace SoLoud
{

enum class RecoveryLevel : uint8_t
{
	None = 0,
	DeviceRestart = 1, // try stop/start device
	DeviceReinit = 2,  // reinitialize device with existing context
	ContextReinit = 3  // reinitialize entire context + device
};

struct MiniaudioData
{
	ma_context context{};
	ma_device device{};
	ma_log log{};
	bool contextInitialized{false};
	bool deviceInitialized{false};
	bool logInitialized{false};
	ma_uint32 maxLogLevel{0};

	// device state tracking
	std::atomic<bool> deviceValid{true};
	std::atomic<bool> shutdownRequested{false};
	std::atomic<bool> shutdownInProgress{false};

	// recovery thread management
	std::atomic<bool> recoveryInProgress{false};
	std::thread recoveryThread;
	std::mutex recoveryMutex;
	std::condition_variable recoveryCondition;
	std::atomic<RecoveryLevel> pendingRecoveryLevel{RecoveryLevel::None};
	std::chrono::steady_clock::time_point lastRecoveryAttempt{};

	// synchronization for shutdown
	std::mutex miniaudioMutex; // protects all miniaudio operations

	// original configuration for recovery
	ma_device_config originalConfig{};
	ma_context_config originalContextConfig{};
	Soloud *soloudInstance{nullptr};

	// initialization parameters for context recovery
	unsigned int initFlags{0};
	unsigned int initSamplerate{0};
	unsigned int initBuffer{0};
	unsigned int initChannels{0};

	MiniaudioData() = default;

	~MiniaudioData() { shutdownRecoveryThread(); }

	void shutdownRecoveryThread()
	{
		if (recoveryThread.joinable())
		{
			shutdownRequested.store(true);
			recoveryCondition.notify_all();
			recoveryThread.join();
		}
	}
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

void soloud_miniaudio_notification_callback(const ma_device_notification *pNotification)
{
	auto *data = static_cast<MiniaudioData *>(pNotification->pDevice->pUserData);
	if (!data || data->shutdownInProgress.load())
		return;

	RecoveryLevel requiredLevel = RecoveryLevel::None;

	switch (pNotification->type)
	{
	case ma_device_notification_type_stopped:
		// device stopped - might be recoverable with restart
		if (data->maxLogLevel >= MA_LOG_LEVEL_INFO)
			fprintf(stderr, "[MiniAudio INFO] Device stopped notification\n");
		requiredLevel = RecoveryLevel::DeviceRestart;
		break;

	case ma_device_notification_type_rerouted:
		// device rerouted - try device reinit first
		if (data->maxLogLevel >= MA_LOG_LEVEL_INFO)
			fprintf(stderr, "[MiniAudio INFO] Device rerouted notification\n");
		requiredLevel = RecoveryLevel::DeviceReinit;
		break;

	case ma_device_notification_type_interruption_began:
		// interruption - might indicate server death, escalate to context reinit
		if (data->maxLogLevel >= MA_LOG_LEVEL_INFO)
			fprintf(stderr, "[MiniAudio INFO] Device interruption began notification\n");
		requiredLevel = RecoveryLevel::ContextReinit;
		break;

	case ma_device_notification_type_interruption_ended:
		// interruption ended - device might be working again
		if (data->maxLogLevel >= MA_LOG_LEVEL_INFO)
			fprintf(stderr, "[MiniAudio INFO] Device interruption ended notification\n");
		// don't mark as invalid, just log
		return;

	case ma_device_notification_type_unlocked:
		if (data->maxLogLevel >= MA_LOG_LEVEL_INFO)
			fprintf(stderr, "[MiniAudio INFO] Device unlocked notification\n");
		return;

	case ma_device_notification_type_started:
		if (data->maxLogLevel >= MA_LOG_LEVEL_INFO)
			fprintf(stderr, "[MiniAudio INFO] Device started notification\n");
		data->deviceValid.store(true);
		return;
	}

	// mark device as invalid and request recovery
	data->deviceValid.store(false);

	// request recovery at the appropriate level, but only if not shutting down
	if (!data->shutdownInProgress.load())
	{
		std::lock_guard<std::mutex> lock(data->recoveryMutex);
		if (requiredLevel > data->pendingRecoveryLevel.load())
		{
			data->pendingRecoveryLevel.store(requiredLevel);
			data->recoveryCondition.notify_one();
		}
	}
}

bool attempt_device_restart(MiniaudioData *data)
{
	if (!data->deviceInitialized || data->shutdownInProgress.load())
		return false;

	if (data->maxLogLevel >= MA_LOG_LEVEL_INFO)
		fprintf(stderr, "[MiniAudio INFO] Attempting device restart\n");

	std::lock_guard<std::mutex> lock(data->miniaudioMutex);
	if (data->shutdownInProgress.load())
		return false;

	ma_result result = ma_device_stop(&data->device);
	if (result != MA_SUCCESS && data->maxLogLevel >= MA_LOG_LEVEL_WARNING)
		fprintf(stderr, "[MiniAudio WARNING] Failed to stop device during restart: %d\n", result);

	if (data->shutdownInProgress.load())
		return false;

	result = ma_device_start(&data->device);
	if (result == MA_SUCCESS)
	{
		if (data->maxLogLevel >= MA_LOG_LEVEL_INFO)
			fprintf(stderr, "[MiniAudio INFO] Device restart successful\n");
		return true;
	}
	else if (data->maxLogLevel >= MA_LOG_LEVEL_WARNING)
		fprintf(stderr, "[MiniAudio WARNING] Failed to start device during restart: %d\n", result);

	return false;
}

bool attempt_device_reinit(MiniaudioData *data)
{
	if (!data->contextInitialized || data->shutdownInProgress.load())
		return false;

	if (data->maxLogLevel >= MA_LOG_LEVEL_INFO)
		fprintf(stderr, "[MiniAudio INFO] Attempting device reinit\n");

	std::lock_guard<std::mutex> lock(data->miniaudioMutex);
	if (data->shutdownInProgress.load())
		return false;

	// uninitialize existing device
	if (data->deviceInitialized)
	{
		ma_device_uninit(&data->device);
		data->deviceInitialized = false;
	}

	if (data->shutdownInProgress.load())
		return false;

	// reinitialize device with existing context
	ma_result result = ma_device_init(&data->context, &data->originalConfig, &data->device);
	if (result == MA_SUCCESS)
	{
		data->deviceInitialized = true;

		if (data->shutdownInProgress.load())
		{
			ma_device_uninit(&data->device);
			data->deviceInitialized = false;
			return false;
		}

		result = ma_device_start(&data->device);
		if (result == MA_SUCCESS)
		{
			if (data->maxLogLevel >= MA_LOG_LEVEL_INFO)
				fprintf(stderr, "[MiniAudio INFO] Device reinit successful\n");
			return true;
		}
		else if (data->maxLogLevel >= MA_LOG_LEVEL_WARNING)
			fprintf(stderr, "[MiniAudio WARNING] Failed to start reinitialized device: %d\n", result);
	}
	else if (data->maxLogLevel >= MA_LOG_LEVEL_WARNING)
		fprintf(stderr, "[MiniAudio WARNING] Failed to reinitialize device: %d\n", result);

	return false;
}

bool attempt_context_reinit(MiniaudioData *data)
{
	if (data->shutdownInProgress.load())
		return false;

	if (data->maxLogLevel >= MA_LOG_LEVEL_INFO)
		fprintf(stderr, "[MiniAudio INFO] Attempting context reinit\n");

	std::lock_guard<std::mutex> lock(data->miniaudioMutex);
	if (data->shutdownInProgress.load())
		return false;

	// uninitialize existing device
	if (data->deviceInitialized)
	{
		ma_device_uninit(&data->device);
		data->deviceInitialized = false;
	}

	if (data->shutdownInProgress.load())
		return false;

	// uninitialize existing context
	if (data->contextInitialized)
	{
		ma_context_uninit(&data->context);
		data->contextInitialized = false;
	}

	if (data->shutdownInProgress.load())
		return false;

	// reinitialize context
	ma_result result = ma_context_init(nullptr, 0, &data->originalContextConfig, &data->context);
	if (result != MA_SUCCESS)
	{
		if (data->maxLogLevel >= MA_LOG_LEVEL_WARNING)
			fprintf(stderr, "[MiniAudio WARNING] Failed to reinitialize context: %d\n", result);
		return false;
	}
	data->contextInitialized = true;

	if (data->shutdownInProgress.load())
		return false;

	// reinitialize device
	result = ma_device_init(&data->context, &data->originalConfig, &data->device);
	if (result != MA_SUCCESS)
	{
		if (data->maxLogLevel >= MA_LOG_LEVEL_WARNING)
			fprintf(stderr, "[MiniAudio WARNING] Failed to reinitialize device after context reinit: %d\n", result);
		return false;
	}
	data->deviceInitialized = true;

	if (data->shutdownInProgress.load())
		return false;

	// start device
	result = ma_device_start(&data->device);
	if (result == MA_SUCCESS)
	{
		if (data->maxLogLevel >= MA_LOG_LEVEL_INFO)
			fprintf(stderr, "[MiniAudio INFO] Context reinit successful\n");

		// update soloud with potentially new device parameters
		if (data->soloudInstance)
		{
			unsigned int actualSampleRate = data->device.sampleRate;
			unsigned int actualBufferSize = data->device.playback.internalPeriodSizeInFrames;
			unsigned int actualChannels = data->device.playback.channels;
			data->soloudInstance->postinit_internal(actualSampleRate, actualBufferSize, data->initFlags, actualChannels);
		}

		return true;
	}
	else if (data->maxLogLevel >= MA_LOG_LEVEL_WARNING)
		fprintf(stderr, "[MiniAudio WARNING] Failed to start device after context reinit: %d\n", result);

	return false;
}

void recovery_thread_function(MiniaudioData *data)
{
	if (data->maxLogLevel >= MA_LOG_LEVEL_INFO)
		fprintf(stderr, "[MiniAudio INFO] Recovery thread started\n");

	while (!data->shutdownRequested.load())
	{
		RecoveryLevel recoveryLevel = RecoveryLevel::None;

		// wait for recovery request or shutdown
		{
			std::unique_lock<std::mutex> lock(data->recoveryMutex);
			data->recoveryCondition.wait(lock, [data]() { return data->shutdownRequested.load() || data->pendingRecoveryLevel.load() != RecoveryLevel::None; });

			if (data->shutdownRequested.load())
				break;

			recoveryLevel = data->pendingRecoveryLevel.exchange(RecoveryLevel::None);
		}

		if (recoveryLevel == RecoveryLevel::None || data->shutdownInProgress.load())
			continue;

		// check if we should throttle recovery attempts
		auto now = std::chrono::steady_clock::now();
		if (now - data->lastRecoveryAttempt < std::chrono::seconds(2))
		{
			std::this_thread::sleep_for(std::chrono::seconds(1));
			continue;
		}
		data->lastRecoveryAttempt = now;

		data->recoveryInProgress.store(true);

		bool success = false;
		switch (recoveryLevel)
		{
		case RecoveryLevel::DeviceRestart:
			success = attempt_device_restart(data);
			if (!success && data->maxLogLevel >= MA_LOG_LEVEL_INFO && !data->shutdownInProgress.load())
				fprintf(stderr, "[MiniAudio INFO] Device restart failed, escalating to device reinit\n");
			// fall through to device reinit if restart failed
			if (!success && !data->shutdownInProgress.load())
				success = attempt_device_reinit(data);
			break;

		case RecoveryLevel::DeviceReinit:
			success = attempt_device_reinit(data);
			if (!success && data->maxLogLevel >= MA_LOG_LEVEL_INFO && !data->shutdownInProgress.load())
				fprintf(stderr, "[MiniAudio INFO] Device reinit failed, escalating to context reinit\n");
			// fall through to context reinit if device reinit failed
			if (!success && !data->shutdownInProgress.load())
				success = attempt_context_reinit(data);
			break;

		case RecoveryLevel::ContextReinit:
			success = attempt_context_reinit(data);
			break;

		default:
			break;
		}

		if (success)
		{
			data->deviceValid.store(true);
		}
		else if (data->maxLogLevel >= MA_LOG_LEVEL_WARNING && !data->shutdownInProgress.load())
		{
			fprintf(stderr, "[MiniAudio WARNING] All recovery attempts failed\n");
		}

		data->recoveryInProgress.store(false);
	}

	if (data->maxLogLevel >= MA_LOG_LEVEL_INFO)
		fprintf(stderr, "[MiniAudio INFO] Recovery thread stopped\n");
}

void soloud_miniaudio_audiomixer(ma_device *pDevice, void *pOutput, const void * /*pInput*/, ma_uint32 frameCount)
{
	auto *data = static_cast<MiniaudioData *>(pDevice->pUserData);
	auto *soloud = data->soloudInstance;

	// check if device is valid or shutting down
	if (!data->deviceValid.load(std::memory_order_relaxed) || data->shutdownInProgress.load(std::memory_order_relaxed))
	{
		// output silence while device is invalid or shutting down
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
	switch (pDevice->playback.internalFormat)
	{
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
	case ma_format_f32:
	default: // fallback to float if format is unknown or unsupported
		soloud->mix(static_cast<float *>(pOutput), frameCount);
		break;
	}
}

void soloud_miniaudio_deinit(SoLoud::Soloud *aSoloud)
{
	auto *data = static_cast<MiniaudioData *>(aSoloud->mBackendData);
	if (data)
	{
		// signal shutdown to prevent new recovery operations
		data->shutdownInProgress.store(true);
		data->deviceValid.store(false);

		// stop device first to stop audio callback
		{
			std::lock_guard<std::mutex> lock(data->miniaudioMutex);
			if (data->deviceInitialized)
			{
				ma_device_stop(&data->device);
			}
		}

		// now signal recovery thread shutdown and wait for it
		data->shutdownRequested.store(true);
		data->recoveryCondition.notify_all();
		data->shutdownRecoveryThread();

		// now safely uninitialize miniaudio objects
		{
			std::lock_guard<std::mutex> lock(data->miniaudioMutex);

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
	if (data && data->deviceInitialized && data->deviceValid.load() && !data->shutdownInProgress.load())
	{
		std::lock_guard<std::mutex> lock(data->miniaudioMutex);
		if (data->deviceInitialized && !data->shutdownInProgress.load())
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
	if (data && data->deviceInitialized && data->deviceValid.load() && !data->shutdownInProgress.load())
	{
		std::lock_guard<std::mutex> lock(data->miniaudioMutex);
		if (data->deviceInitialized && !data->shutdownInProgress.load())
		{
			if (ma_device_start(&data->device) != MA_SUCCESS)
				return UNKNOWN_ERROR;
		}
	}
	return SO_NO_ERROR;
}

} // namespace

result miniaudio_init(SoLoud::Soloud *aSoloud, unsigned int aFlags, unsigned int aSamplerate, unsigned int aBuffer, unsigned int aChannels)
{
	auto *data = new MiniaudioData();
	aSoloud->mBackendData = data;
	data->soloudInstance = aSoloud;
	data->initFlags = aFlags;
	data->initSamplerate = aSamplerate;
	data->initBuffer = aBuffer;
	data->initChannels = aChannels;

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

	// store original context config for recovery
	data->originalContextConfig = contextConfig;

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

	// store original config for recovery
	data->originalConfig = config;

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
		case MA_DEVICE_ALREADY_INITIALIZED:
		case MA_DEVICE_NOT_STARTED:
		case MA_DEVICE_NOT_STOPPED:
		default:
			return UNKNOWN_ERROR;
		}
	}
	data->deviceInitialized = true;

	// start recovery thread before starting device
	data->recoveryThread = std::thread(recovery_thread_function, data);

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
