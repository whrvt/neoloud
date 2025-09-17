/*
SoLoud audio engine
Copyright (c) 2013-2015 Jari Komppa

SDL3 static entry point
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

#include "soloud_internal.h"

#if !defined(WITH_SDL3)

namespace SoLoud
{
result sdl3_init(SoLoud::Soloud *aSoloud, unsigned int aFlags, unsigned int aSamplerate, unsigned int aBuffer, unsigned int aChannels)
{
	return NOT_IMPLEMENTED;
}
} // namespace SoLoud

#else

#include <SDL3/SDL_audio.h>
#include <SDL3/SDL_hints.h>
#include <SDL3/SDL_init.h>
#include <SDL3/SDL_log.h>
#include <SDL3/SDL_version.h>

#include <array>
#include <atomic>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <string>
#include <vector>

#ifdef _MSC_VER
#ifndef strcasecmp
#define strcasecmp _stricmp
#endif
#endif

namespace SoLoud
{
using namespace detail; // SAMPLE_FORMAT

struct SDL3Data
{
	SDL_AudioStream *audioStream{nullptr};
	SDL_AudioSpec streamSpec{}; // stream input format (what we provide)
	SDL_AudioDeviceID deviceID{0};
	SDL_AudioSpec deviceSpec{}; // actual device format

	struct
	{
		unsigned int flags{Soloud::CLIP_ROUNDOFF};
		unsigned int samplerate{Soloud::AUTO};
		unsigned int bufferSize{Soloud::AUTO};
		unsigned int channels{2};
	} initParams; // parameters we were first created with

	// device management
	SDL_AudioDeviceID currentPhysicalDeviceID{SDL_AUDIO_DEVICE_DEFAULT_PLAYBACK}; // which physical device we opened
	DeviceInfo currentDeviceInfo{};
	bool hasCurrentDeviceInfo{false};

	// for init/shutdown
	std::mutex deviceMutex;

	// mix buffer to avoid allocation in audio callback
#if SDL_VERSION_ATLEAST(3, 3, 0)
	// switch between a few separate buffers so we can use SDL_PutAudioStreamNoCopy
	// in practice, SDL only keeps up to 2 buffers alive at once,
	// but using some extras to be safe (the documentation is unclear on if there is a hard limit or not)
	static constexpr const size_t MAX_MIXBUFFERS{4};

	std::array<std::vector<uint8_t>, MAX_MIXBUFFERS> mixBuffers;
	std::array<bool, MAX_MIXBUFFERS> mixBufNumBusy{}; // mostly for debugging, keep track of which buffer SDL doesn't need anymore
	uint8_t cb{0};                                    // current mix buffer
#else
	std::vector<uint8_t> mixBuffer;
#endif

	// parent instance
	Soloud *soloudInstance{nullptr};

	std::atomic<bool> deviceValid{true};
	std::atomic<bool> soloudInitialized{false}; // postinit_internal must be called before we use soloud->mix in the callback
	bool weInitSDLAudio{false};
	bool streamInitialized{false};
	bool deviceInitialized{false};
};

namespace // static
{

#if SDL_VERSION_ATLEAST(3, 3, 0)
void nocopy_complete_callback(void *userdata, const void * /*buf*/, int /*buflen*/)
{
	if (userdata)
	{
		bool *busy = static_cast<bool *>(userdata);
		*busy = false;
	}
}
#endif

void soloud_sdl3_stream_callback(void *pUserData, SDL_AudioStream *stream, int additional_amount, int total_amount)
{
	auto *data = static_cast<SDL3Data *>(pUserData);
	auto *soloud = data->soloudInstance;

#if SDL_VERSION_ATLEAST(3, 3, 0)
	auto &currentBuf = data->mixBuffers[data->cb];
#else
	auto &currentBuf = data->mixBuffer;
#endif

	// use additional_amount since that's what's actually needed right now
	// total_amount might include data already queued in the stream
	int bytesToProvide = (additional_amount > 0) ? additional_amount : total_amount;
	if (bytesToProvide <= 0)
		return;

#if SDL_VERSION_ATLEAST(3, 3, 0)
	if (data->mixBufNumBusy[data->cb])
		SDL_LogError(SDL_LOG_CATEGORY_AUDIO, "Mix buffer %u was in use by SDL before PutAudioStreamDataNoCopy!", data->cb);

	// mark this buffer as in use
	data->mixBufNumBusy[data->cb] = true;
#endif

	// calculate how many frames we need
	int frameSize = SDL_AUDIO_FRAMESIZE(data->streamSpec);
	int requestedFrames = bytesToProvide / frameSize;

	// ensure our mix buffer is large enough
	if (currentBuf.size() < static_cast<size_t>(bytesToProvide))
		currentBuf.resize(static_cast<size_t>(bytesToProvide));

	// determine output format for soloud mixer to match stream input format exactly
	SAMPLE_FORMAT outputFormat = SAMPLE_FLOAT32;
	if (data->streamSpec.format == SDL_AUDIO_S16)
		outputFormat = SAMPLE_SIGNED16;
	else if (data->streamSpec.format == SDL_AUDIO_S32)
		outputFormat = SAMPLE_SIGNED32;
	else if (data->streamSpec.format == SDL_AUDIO_U8)
		outputFormat = SAMPLE_UNSIGNED8;

	if (data->soloudInitialized.load(std::memory_order_acquire) && data->deviceValid.load(std::memory_order_relaxed))
	{
		// mix audio directly into our buffer in the correct format
		soloud->mix(currentBuf.data(), requestedFrames, outputFormat);
	}
	else
	{
		// device not valid, fill with silence
		int silenceValue = SDL_GetSilenceValueForFormat(data->streamSpec.format);
		memset(currentBuf.data(), silenceValue, bytesToProvide);
	}

	// provide data to the stream
#if SDL_VERSION_ATLEAST(3, 3, 0)
	SDL_PutAudioStreamDataNoCopy(stream, currentBuf.data(), bytesToProvide, nocopy_complete_callback, &data->mixBufNumBusy[data->cb]);

	// switch buffers
	data->cb++;
	data->cb = data->cb % data->MAX_MIXBUFFERS;
#else
	SDL_PutAudioStreamData(stream, currentBuf.data(), bytesToProvide);
#endif
}

void soloud_sdl3_deinit(SoLoud::Soloud *aSoloud)
{
	auto *data = static_cast<SDL3Data *>(aSoloud->mBackendData);
	if (data)
	{
		data->deviceValid.store(false);

		{
			std::lock_guard<std::mutex> lock(data->deviceMutex);

			if (data->streamInitialized && data->audioStream)
			{
				SDL_DestroyAudioStream(data->audioStream);
				data->audioStream = nullptr;
				data->streamInitialized = false;
			}

			if (data->deviceInitialized && data->deviceID)
			{
				SDL_CloseAudioDevice(data->deviceID);
				data->deviceID = 0;
				data->deviceInitialized = false;
			}
		}

		if (data->weInitSDLAudio && SDL_WasInit(SDL_INIT_AUDIO))
		{
			SDL_QuitSubSystem(SDL_INIT_AUDIO);
			data->weInitSDLAudio = false;
		}

		delete data;
		aSoloud->mBackendData = nullptr;
	}
}

result soloud_sdl3_pause(SoLoud::Soloud *aSoloud)
{
	auto *data = static_cast<SDL3Data *>(aSoloud->mBackendData);
	if (data && data->deviceInitialized)
	{
		std::lock_guard<std::mutex> lock(data->deviceMutex);
		if (data->deviceInitialized && data->deviceID)
		{
			if (!SDL_PauseAudioDevice(data->deviceID))
			{
				SDL_LogError(SDL_LOG_CATEGORY_AUDIO, "Failed to pause audio device: %s", SDL_GetError());
				return UNKNOWN_ERROR;
			}
			data->deviceValid.store(false);
		}
	}
	return SO_NO_ERROR;
}

result soloud_sdl3_resume(SoLoud::Soloud *aSoloud)
{
	auto *data = static_cast<SDL3Data *>(aSoloud->mBackendData);
	if (data && data->deviceInitialized)
	{
		std::lock_guard<std::mutex> lock(data->deviceMutex);
		if (data->deviceInitialized && data->deviceID)
		{
			if (!SDL_ResumeAudioDevice(data->deviceID))
			{
				SDL_LogError(SDL_LOG_CATEGORY_AUDIO, "Failed to resume audio device: %s", SDL_GetError());
				return UNKNOWN_ERROR;
			}
			data->deviceValid.store(true);
		}
	}
	return SO_NO_ERROR;
}

SDL_LogPriority parse_log_level_from_env()
{
	const char *env = getenv("SOLOUD_DEBUG");
	if (!env)
	{
#ifdef _DEBUG
		return SDL_LOG_PRIORITY_WARN;
#else
		return static_cast<SDL_LogPriority>(SDL_LOG_PRIORITY_CRITICAL + 1); // disable
#endif
	}

	if (strcmp(env, "debug") == 0)
		return SDL_LOG_PRIORITY_TRACE;
	if (strcmp(env, "info") == 0)
		return SDL_LOG_PRIORITY_INFO;
	if (strcmp(env, "warn") == 0)
		return SDL_LOG_PRIORITY_WARN;
	if (strcmp(env, "error") == 0)
		return SDL_LOG_PRIORITY_ERROR;
	if (strcmp(env, "none") == 0)
		return static_cast<SDL_LogPriority>(SDL_LOG_PRIORITY_CRITICAL + 1);

	// default fallback for unrecognized values
#ifdef _DEBUG
	return SDL_LOG_PRIORITY_WARN;
#else
	return static_cast<SDL_LogPriority>(SDL_LOG_PRIORITY_CRITICAL + 1);
#endif
}

// try to get an optimal set of parameters/device spec for a given device
void params_to_optimal_devconfig(SDL3Data *instance, SDL_AudioSpec *outSpec, unsigned int *outBufsize, SDL_AudioDeviceID devid = SDL_AUDIO_DEVICE_DEFAULT_PLAYBACK)
{
	SDL_AudioSpec deviceSpec;
	int deviceFrames = 0;
	bool hasDeviceInfo = SDL_GetAudioDeviceFormat(devid, &deviceSpec, &deviceFrames);

	// configure device specification - prefer device native when using AUTO settings
	if (instance->initParams.samplerate == Soloud::AUTO && hasDeviceInfo)
		outSpec->freq = deviceSpec.freq;
	else
		outSpec->freq = (instance->initParams.samplerate > 0) ? instance->initParams.samplerate : 44100;

	if (instance->initParams.channels == Soloud::AUTO && hasDeviceInfo)
		outSpec->channels = deviceSpec.channels;
	else
		outSpec->channels = (instance->initParams.channels > 0) ? instance->initParams.channels : 2;

	// prefer device native format to avoid unnecessary conversion
	if (hasDeviceInfo &&
	    (deviceSpec.format == SDL_AUDIO_F32 || deviceSpec.format == SDL_AUDIO_S32 || deviceSpec.format == SDL_AUDIO_S16 || deviceSpec.format == SDL_AUDIO_U8))
		outSpec->format = deviceSpec.format;
	else // fallback to float32
		outSpec->format = SDL_AUDIO_F32;

	// determine buffer size and set hint
	if (instance->initParams.bufferSize == Soloud::AUTO) // aim for ~1ms at target sample rate for low latency (ceil'd)
		*outBufsize = (outSpec->freq + 999) / 1000;
	else if (hasDeviceInfo && deviceFrames > 0) // use device buffer size (preferring explicit request if smaller)
		*outBufsize = deviceFrames > instance->initParams.bufferSize ? instance->initParams.bufferSize : deviceFrames;
	else
		*outBufsize = instance->initParams.bufferSize;

	const char *driver = SDL_GetCurrentAudioDriver();

#ifdef WINDOWS_VERSION
	// unless it was explicitly requested to be lower, clamp to min 10ms on Windows if not using WASAPI (empirically, directsound blows up with anything lower)
	if (strcasecmp(driver, "wasapi") != 0)
		*outBufsize = std::max(((outSpec->freq * 10) + 999) / 1000u, *outBufsize);
#elif defined(__linux__)
	if (strcasecmp(driver, "pulseaudio") == 0) // SDL seems to always return a buffer size 2x of what was requested (but only for pulseaudio)... why?
		*outBufsize = (*outBufsize + 1) / 2;   // curse integer truncation
#else // other platforms/drivers, just clamp to min 10ms (like windows non-wasapi)
	*outBufsize = std::max(((outSpec->freq * 10) + 999) / 1000u, *outBufsize);
#endif

	return;
}

// device identifier conversion helpers
std::string device_id_to_identifier(SDL_AudioDeviceID deviceID)
{
	if (deviceID == SDL_AUDIO_DEVICE_DEFAULT_PLAYBACK)
		return "sdl3_default_playback";
	return "sdl3_" + std::to_string(deviceID);
}

SDL_AudioDeviceID identifier_to_device_id(const char *identifier)
{
	if (!identifier)
		return SDL_AUDIO_DEVICE_DEFAULT_PLAYBACK;

	if (strcmp(identifier, "sdl3_default_playback") == 0)
		return SDL_AUDIO_DEVICE_DEFAULT_PLAYBACK;

	if (strncmp(identifier, "sdl3_", 5) == 0)
	{
		char *endptr;
		unsigned long id = strtoul(identifier + 5, &endptr, 10);
		if (*endptr == '\0' && id <= UINT32_MAX)
			return static_cast<SDL_AudioDeviceID>(id);
	}

	return 0; // invalid
}

// update current device info tracking
void update_current_device_info(SDL3Data *data)
{
	const char *deviceName = nullptr;

	if (data->currentPhysicalDeviceID == SDL_AUDIO_DEVICE_DEFAULT_PLAYBACK)
	{
		deviceName = "Default Playback Device";
		data->currentDeviceInfo.isDefault = true;
	}
	else
	{
		deviceName = SDL_GetAudioDeviceName(data->currentPhysicalDeviceID);
		if (!deviceName)
			deviceName = "Unknown Device";
		data->currentDeviceInfo.isDefault = false;
	}

	// copy name with bounds checking
	size_t nameLen = strlen(deviceName);
	if (nameLen >= sizeof(data->currentDeviceInfo.name))
		nameLen = sizeof(data->currentDeviceInfo.name) - 1;
	memcpy(data->currentDeviceInfo.name, deviceName, nameLen);
	data->currentDeviceInfo.name[nameLen] = '\0';

	// set identifier
	std::string identifier = device_id_to_identifier(data->currentPhysicalDeviceID);
	size_t identifierLen = identifier.length();
	if (identifierLen >= sizeof(data->currentDeviceInfo.identifier))
		identifierLen = sizeof(data->currentDeviceInfo.identifier) - 1;
	memcpy(data->currentDeviceInfo.identifier, identifier.c_str(), identifierLen);
	data->currentDeviceInfo.identifier[identifierLen] = '\0';

	data->currentDeviceInfo.nativeDeviceInfo = nullptr;
	data->hasCurrentDeviceInfo = true;
}

// common device setup logic shared between init and device switching
result setup_device_and_stream(SDL3Data *data, SDL_AudioDeviceID targetDeviceID, SDL_AudioSpec *outDeviceSpec, int *outDeviceFrames)
{
	// get optimal configuration for target device
	SDL_AudioSpec requestSpec;
	unsigned int requestBufsize;
	params_to_optimal_devconfig(data, &requestSpec, &requestBufsize, targetDeviceID);

	// set buffer size hint
	std::string bufferStr = std::to_string(requestBufsize);
	SDL_SetHintWithPriority(SDL_HINT_AUDIO_DEVICE_SAMPLE_FRAMES, bufferStr.c_str(), SDL_HINT_NORMAL);

	SDL_LogDebug(SDL_LOG_CATEGORY_AUDIO, "Opening device %s: %dHz, %d channels, format %#x, buffer %u frames", device_id_to_identifier(targetDeviceID).c_str(),
	             requestSpec.freq, requestSpec.channels, requestSpec.format, requestBufsize);

	// open device
	data->deviceID = SDL_OpenAudioDevice(targetDeviceID, &requestSpec);
	if (!data->deviceID)
	{
		SDL_LogError(SDL_LOG_CATEGORY_AUDIO, "Failed to open audio device: %s", SDL_GetError());
		return UNKNOWN_ERROR;
	}
	data->deviceInitialized = true;

	// get actual device format that was negotiated
	int actualDeviceFrames = 0;
	if (!SDL_GetAudioDeviceFormat(data->deviceID, &data->deviceSpec, &actualDeviceFrames))
	{
		SDL_LogError(SDL_LOG_CATEGORY_AUDIO, "Failed to get device format: %s", SDL_GetError());
		return UNKNOWN_ERROR;
	}

	SDL_LogDebug(SDL_LOG_CATEGORY_AUDIO, "Device opened with: %dHz, %d channels, format %#x, buffer %d frames", data->deviceSpec.freq, data->deviceSpec.channels,
	             data->deviceSpec.format, actualDeviceFrames);

	// update tracking
	data->currentPhysicalDeviceID = targetDeviceID;
	update_current_device_info(data);

	// use same format as device to avoid conversion
	data->streamSpec = data->deviceSpec;

	// create stream
	data->audioStream = SDL_CreateAudioStream(&data->streamSpec, &data->deviceSpec);
	if (!data->audioStream)
	{
		SDL_LogError(SDL_LOG_CATEGORY_AUDIO, "Failed to create audio stream: %s", SDL_GetError());
		return UNKNOWN_ERROR;
	}
	data->streamInitialized = true;

	// bind stream to device
	if (!SDL_BindAudioStream(data->deviceID, data->audioStream))
	{
		SDL_LogError(SDL_LOG_CATEGORY_AUDIO, "Failed to bind stream to device: %s", SDL_GetError());
		return UNKNOWN_ERROR;
	}

	// set up callback
	if (!SDL_SetAudioStreamGetCallback(data->audioStream, soloud_sdl3_stream_callback, data))
	{
		SDL_LogError(SDL_LOG_CATEGORY_AUDIO, "Failed to set stream callback: %s", SDL_GetError());
		return UNKNOWN_ERROR;
	}

	// resize mix buffers if needed
	size_t maxBufferBytes = static_cast<size_t>(actualDeviceFrames) * data->streamSpec.channels * SDL_AUDIO_BYTESIZE(data->streamSpec.format);
#if SDL_VERSION_ATLEAST(3, 3, 0)
	for (auto &buf : data->mixBuffers)
	{
		if (buf.size() < maxBufferBytes)
			buf.resize(maxBufferBytes);
	}
#else
	if (data->mixBuffer.size() < maxBufferBytes)
		data->mixBuffer.resize(maxBufferBytes);
#endif

	// return actual configuration for caller
	if (outDeviceSpec)
		*outDeviceSpec = data->deviceSpec;
	if (outDeviceFrames)
		*outDeviceFrames = actualDeviceFrames;

	return SO_NO_ERROR;
}

// device management functions
result sdl3_enumerate_devices(Soloud *aSoloud)
{
	if (!aSoloud)
		return INVALID_PARAMETER;

	auto *data = static_cast<SDL3Data *>(aSoloud->mBackendData);
	if (!data)
		return INVALID_PARAMETER;

	std::lock_guard<std::mutex> lock(data->deviceMutex);

	int physicalDeviceCount = 0;
	SDL_AudioDeviceID *physicalDevices = SDL_GetAudioPlaybackDevices(&physicalDeviceCount);
	if (!physicalDevices)
	{
		SDL_LogError(SDL_LOG_CATEGORY_AUDIO, "Failed to enumerate audio devices: %s", SDL_GetError());
		return UNKNOWN_ERROR;
	}

	// allocate device list (including default device entry)
	aSoloud->mDeviceList = new DeviceInfo[physicalDeviceCount + 1];
	if (!aSoloud->mDeviceList)
	{
		SDL_free(physicalDevices);
		aSoloud->mDeviceCount = 0;
		return OUT_OF_MEMORY;
	}

	// add default device first
	aSoloud->mDeviceList[0] = {};
	strncpy(aSoloud->mDeviceList[0].name, "Default Playback Device", sizeof(aSoloud->mDeviceList[0].name) - 1);
	strncpy(aSoloud->mDeviceList[0].identifier, "sdl3_default_playback", sizeof(aSoloud->mDeviceList[0].identifier) - 1);
	aSoloud->mDeviceList[0].isDefault = true;
	aSoloud->mDeviceList[0].nativeDeviceInfo = nullptr;

	// add physical devices
	for (int i = 0; i < physicalDeviceCount; i++)
	{
		const char *deviceName = SDL_GetAudioDeviceName(physicalDevices[i]);
		if (!deviceName)
			deviceName = "Unknown Device";

		aSoloud->mDeviceList[i + 1] = {};

		// copy name with bounds checking
		size_t nameLen = strlen(deviceName);
		if (nameLen >= sizeof(aSoloud->mDeviceList[i + 1].name))
			nameLen = sizeof(aSoloud->mDeviceList[i + 1].name) - 1;
		memcpy(aSoloud->mDeviceList[i + 1].name, deviceName, nameLen);
		aSoloud->mDeviceList[i + 1].name[nameLen] = '\0';

		// set identifier
		std::string identifier = device_id_to_identifier(physicalDevices[i]);
		size_t identifierLen = identifier.length();
		if (identifierLen >= sizeof(aSoloud->mDeviceList[i + 1].identifier))
			identifierLen = sizeof(aSoloud->mDeviceList[i + 1].identifier) - 1;
		memcpy(aSoloud->mDeviceList[i + 1].identifier, identifier.c_str(), identifierLen);
		aSoloud->mDeviceList[i + 1].identifier[identifierLen] = '\0';

		aSoloud->mDeviceList[i + 1].isDefault = false;
		aSoloud->mDeviceList[i + 1].nativeDeviceInfo = nullptr;
	}

	aSoloud->mDeviceCount = physicalDeviceCount + 1;
	SDL_free(physicalDevices);
	return SO_NO_ERROR;
}

result sdl3_get_current_device(Soloud *aSoloud, DeviceInfo *pDeviceInfo)
{
	if (!aSoloud || !pDeviceInfo)
		return INVALID_PARAMETER;

	auto *data = static_cast<SDL3Data *>(aSoloud->mBackendData);
	if (!data || !data->hasCurrentDeviceInfo)
		return INVALID_PARAMETER;

	*pDeviceInfo = data->currentDeviceInfo;
	return SO_NO_ERROR;
}

result sdl3_set_device(Soloud *aSoloud, const char *deviceIdentifier)
{
	if (!aSoloud)
		return INVALID_PARAMETER;

	auto *data = static_cast<SDL3Data *>(aSoloud->mBackendData);
	if (!data)
		return INVALID_PARAMETER;

	// determine target device ID
	SDL_AudioDeviceID targetDeviceID;
	if (!deviceIdentifier || strlen(deviceIdentifier) == 0)
	{
		targetDeviceID = SDL_AUDIO_DEVICE_DEFAULT_PLAYBACK;
	}
	else
	{
		targetDeviceID = identifier_to_device_id(deviceIdentifier);
		if (targetDeviceID == 0)
		{
			SDL_LogError(SDL_LOG_CATEGORY_AUDIO, "Invalid device identifier: %s", deviceIdentifier);
			return INVALID_PARAMETER;
		}
	}

	// if already using this device, no change needed
	if (data->currentPhysicalDeviceID == targetDeviceID)
	{
		SDL_LogDebug(SDL_LOG_CATEGORY_AUDIO, "Already using requested device");
		return SO_NO_ERROR;
	}

	// lock both audio and device mutexes to ensure safe switching
	aSoloud->lockAudioMutex_internal();
	std::lock_guard<std::mutex> deviceLock(data->deviceMutex);

	// store current configuration for comparison
	unsigned int oldSampleRate = data->deviceSpec.freq;
	unsigned int oldChannels = data->deviceSpec.channels;
	unsigned int oldBufferSize = 0;
	if (data->deviceInitialized)
	{
		int deviceFrames = 0;
		SDL_GetAudioDeviceFormat(data->deviceID, nullptr, &deviceFrames);
		oldBufferSize = deviceFrames;
	}

	// mark device as invalid during transition
	data->deviceValid.store(false);

	// destroy current stream and close device
	if (data->streamInitialized && data->audioStream)
	{
		SDL_DestroyAudioStream(data->audioStream);
		data->audioStream = nullptr;
		data->streamInitialized = false;
	}

	if (data->deviceInitialized && data->deviceID)
	{
		SDL_CloseAudioDevice(data->deviceID);
		data->deviceID = 0;
		data->deviceInitialized = false;
	}

	// set up new device and stream
	SDL_AudioSpec newDeviceSpec;
	int newDeviceFrames = 0;
	result setupResult = setup_device_and_stream(data, targetDeviceID, &newDeviceSpec, &newDeviceFrames);
	if (setupResult != SO_NO_ERROR)
	{
		aSoloud->unlockAudioMutex_internal();
		return setupResult;
	}

	// check if device configuration changed
	unsigned int newSampleRate = newDeviceSpec.freq;
	unsigned int newChannels = newDeviceSpec.channels;
	unsigned int newBufferSize = newDeviceFrames;

	bool configChanged = (oldSampleRate != newSampleRate) || (oldChannels != newChannels) || (oldBufferSize != newBufferSize);

	if (configChanged)
	{
		SDL_LogDebug(SDL_LOG_CATEGORY_AUDIO, "Device configuration changed: %uHz/%uch/%uf -> %uHz/%uch/%uf", oldSampleRate, oldChannels, oldBufferSize,
		             newSampleRate, newChannels, newBufferSize);

		// unlock audio mutex before calling postinit_internal as it may need to reinitialize internal state
		aSoloud->unlockAudioMutex_internal();

		// update SoLoud's internal configuration
		aSoloud->postinit_internal(newSampleRate, newBufferSize, data->initParams.flags, newChannels);

		// re-lock for device start
		aSoloud->lockAudioMutex_internal();
	}

	// start the new device
	if (!SDL_ResumeAudioDevice(data->deviceID))
	{
		aSoloud->unlockAudioMutex_internal();
		SDL_LogError(SDL_LOG_CATEGORY_AUDIO, "Failed to start new audio device: %s", SDL_GetError());
		return UNKNOWN_ERROR;
	}

	data->deviceValid.store(true);
	aSoloud->unlockAudioMutex_internal();

	SDL_LogInfo(SDL_LOG_CATEGORY_AUDIO, "Successfully switched to device: %s", data->currentDeviceInfo.name);
	return SO_NO_ERROR;
}

} // namespace

result sdl3_init(SoLoud::Soloud *aSoloud, unsigned int aFlags /*Soloud::CLIP_ROUNDOFF*/, unsigned int aSamplerate /*Soloud::AUTO (0)*/,
                 unsigned int aBuffer /*Soloud::AUTO (0)*/, unsigned int aChannels /*2*/)
{
	if (!aSoloud)
		return INVALID_PARAMETER;

	auto *data = new SDL3Data();
	aSoloud->mBackendData = data;
	data->soloudInstance = aSoloud;

	// save starting parameters (for later implementing device switching)
	data->initParams = {aFlags, aSamplerate, aBuffer, aChannels};

	// setup logging based on SOLOUD_DEBUG envvar
	SDL_LogPriority logLevel = parse_log_level_from_env();
	if (logLevel <= SDL_LOG_PRIORITY_CRITICAL)
		SDL_SetLogPriority(SDL_LOG_CATEGORY_AUDIO, logLevel);

	SDL_LogDebug(SDL_LOG_CATEGORY_AUDIO, "Initializing SDL3 audio backend");

	const bool sdlEventsWereAlreadyInit = SDL_WasInit(SDL_INIT_EVENTS);

	// initialize SDL audio subsystem if needed
	if (!SDL_WasInit(SDL_INIT_AUDIO))
	{
		if (!SDL_InitSubSystem(SDL_INIT_AUDIO))
		{
			SDL_LogError(SDL_LOG_CATEGORY_AUDIO, "Failed to initialize SDL audio subsystem: %s", SDL_GetError());
			delete data;
			aSoloud->mBackendData = nullptr;
			return UNKNOWN_ERROR;
		}

		// we don't want SDL event handling, which is automatically initialized by SDL_INIT_AUDIO
		if (!sdlEventsWereAlreadyInit)
			SDL_QuitSubSystem(SDL_INIT_EVENTS);

		// used to keep track when exiting if we should quit the audio subsystem (which might have been externally initialized)
		data->weInitSDLAudio = true;
	}

	// set up device and stream
	int actualDeviceFrames = 0;
	result setupResult = setup_device_and_stream(data, SDL_AUDIO_DEVICE_DEFAULT_PLAYBACK, nullptr, &actualDeviceFrames);
	if (setupResult != SO_NO_ERROR)
	{
		soloud_sdl3_deinit(aSoloud);
		return setupResult;
	}

	// log final configuration
	{
		const char *formatName = "Unknown";
		switch (data->streamSpec.format)
		{
		case SDL_AUDIO_U8:
			formatName = "U8";
			break;
		case SDL_AUDIO_S16:
			formatName = "S16";
			break;
		case SDL_AUDIO_S32:
			formatName = "S32";
			break;
		case SDL_AUDIO_F32:
			formatName = "F32";
			break;
		default:
			formatName = "Other";
			break;
		}

		SDL_LogInfo(SDL_LOG_CATEGORY_AUDIO, "Stream created: %dHz, %d channels, %s format, %d frame buffer", data->streamSpec.freq, data->streamSpec.channels,
		            formatName, actualDeviceFrames);
	}

	// initialize SoLoud with our configuration
	aSoloud->postinit_internal(data->streamSpec.freq, actualDeviceFrames, aFlags, data->streamSpec.channels);

	data->soloudInitialized.store(true);

	// set up control functions
	aSoloud->mBackendCleanupFunc = soloud_sdl3_deinit;
	aSoloud->mBackendPauseFunc = soloud_sdl3_pause;
	aSoloud->mBackendResumeFunc = soloud_sdl3_resume;

	// set up device management function pointers
	aSoloud->mEnumerateDevicesFunc = sdl3_enumerate_devices;
	aSoloud->mGetCurrentDeviceFunc = sdl3_get_current_device;
	aSoloud->mSetDeviceFunc = sdl3_set_device;

	aSoloud->mBackendString = "SDL3";

	// start audio playback
	if (!SDL_ResumeAudioDevice(data->deviceID))
	{
		SDL_LogError(SDL_LOG_CATEGORY_AUDIO, "Failed to start audio device: %s", SDL_GetError());
		soloud_sdl3_deinit(aSoloud);
		return UNKNOWN_ERROR;
	}

	SDL_LogInfo(SDL_LOG_CATEGORY_AUDIO, "SDL3 audio backend initialized");

	return SO_NO_ERROR;
}

} // namespace SoLoud

#endif
