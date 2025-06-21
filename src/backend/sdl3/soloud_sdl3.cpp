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

#include <atomic>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <vector>

namespace SoLoud
{
using namespace detail; // SAMPLE_FORMAT

struct SDL3Data
{
	SDL_AudioStream *audioStream{nullptr};
	SDL_AudioDeviceID deviceID{0};
	SDL_AudioSpec deviceSpec{}; // actual device format
	SDL_AudioSpec streamSpec{}; // stream input format (what we provide)
	bool sdlInitialized{false};
	bool streamInitialized{false};
	bool deviceInitialized{false};

	std::atomic<bool> deviceValid{true};
	std::mutex deviceMutex;

	// mix buffer to avoid allocation in audio callback
	std::vector<uint8_t> mixBuffer;
	unsigned int bufferFrames{0};

	Soloud *soloudInstance{nullptr};
};

namespace // static
{

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

void soloud_sdl3_stream_callback(void *pUserData, SDL_AudioStream *stream, int additional_amount, int total_amount)
{
	auto *data = static_cast<SDL3Data *>(pUserData);
	auto *soloud = data->soloudInstance;

	// use additional_amount since that's what's actually needed right now
	// total_amount might include data already queued in the stream
	int bytesToProvide = (additional_amount > 0) ? additional_amount : total_amount;
	if (bytesToProvide <= 0)
		return;

	// calculate how many frames we need
	int frameSize = SDL_AUDIO_FRAMESIZE(data->streamSpec);
	int requestedFrames = bytesToProvide / frameSize;

	// ensure our mix buffer is large enough
	if (data->mixBuffer.size() < static_cast<size_t>(bytesToProvide))
		data->mixBuffer.resize(static_cast<size_t>(bytesToProvide));

	// determine output format for soloud mixer to match stream input format exactly
	SAMPLE_FORMAT outputFormat = SAMPLE_FLOAT32;
	if (data->streamSpec.format == SDL_AUDIO_S16)
		outputFormat = SAMPLE_SIGNED16;
	else if (data->streamSpec.format == SDL_AUDIO_S32)
		outputFormat = SAMPLE_SIGNED32;
	else if (data->streamSpec.format == SDL_AUDIO_U8)
		outputFormat = SAMPLE_UNSIGNED8;

	if (data->deviceValid.load(std::memory_order_relaxed))
	{
		// mix audio directly into our buffer in the correct format
		soloud->mix(data->mixBuffer.data(), requestedFrames, outputFormat);
	}
	else
	{
		// device not valid, fill with silence
		int silenceValue = SDL_GetSilenceValueForFormat(data->streamSpec.format);
		memset(data->mixBuffer.data(), silenceValue, bytesToProvide);
	}

	// provide data to the stream
	SDL_PutAudioStreamData(stream, data->mixBuffer.data(), bytesToProvide);
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

		if (data->sdlInitialized && SDL_WasInit(SDL_INIT_AUDIO))
		{
			SDL_QuitSubSystem(SDL_INIT_AUDIO);
			data->sdlInitialized = false;
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

} // namespace

result sdl3_init(SoLoud::Soloud *aSoloud, unsigned int aFlags, unsigned int aSamplerate, unsigned int aBuffer, unsigned int aChannels)
{
	if (!aSoloud)
		return INVALID_PARAMETER;

	auto *data = new SDL3Data();
	aSoloud->mBackendData = data;
	data->soloudInstance = aSoloud;

	// setup logging based on SOLOUD_DEBUG envvar
	SDL_LogPriority logLevel = parse_log_level_from_env();
	if (logLevel <= SDL_LOG_PRIORITY_CRITICAL)
		SDL_SetLogPriority(SDL_LOG_CATEGORY_AUDIO, logLevel);

	SDL_LogDebug(SDL_LOG_CATEGORY_AUDIO, "Initializing SDL3 audio backend");

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
		data->sdlInitialized = true;
	}

	// get device capabilities to choose optimal format
	SDL_AudioSpec deviceSpec;
	int deviceFrames = 0;
	bool hasDeviceInfo = SDL_GetAudioDeviceFormat(SDL_AUDIO_DEVICE_DEFAULT_PLAYBACK, &deviceSpec, &deviceFrames);

	// configure device specification - prefer device native when using AUTO settings
	SDL_AudioSpec deviceRequestSpec = {};
	if (aSamplerate == Soloud::AUTO && hasDeviceInfo)
		deviceRequestSpec.freq = deviceSpec.freq;
	else
		deviceRequestSpec.freq = (aSamplerate > 0) ? aSamplerate : 44100;

	if (aChannels == Soloud::AUTO && hasDeviceInfo)
		deviceRequestSpec.channels = deviceSpec.channels;
	else
		deviceRequestSpec.channels = (aChannels > 0) ? aChannels : 2;

	// prefer device native format to avoid unnecessary conversion
	if (hasDeviceInfo &&
	    (deviceSpec.format == SDL_AUDIO_F32 || deviceSpec.format == SDL_AUDIO_S32 || deviceSpec.format == SDL_AUDIO_S16 || deviceSpec.format == SDL_AUDIO_U8))
		deviceRequestSpec.format = deviceSpec.format;
	else // fallback to float32
		deviceRequestSpec.format = SDL_AUDIO_F32;

	// determine buffer size and set hint
	unsigned int requestedBufferSize = 0;
	if (aBuffer == Soloud::AUTO) // aim for ~1ms at target sample rate for low latency (ceil'd)
		requestedBufferSize = (deviceRequestSpec.freq + 999) / 1000;
	else if (hasDeviceInfo && deviceFrames > 0) // use device buffer size (preferring explicit request if smaller)
		requestedBufferSize = deviceFrames > aBuffer ? aBuffer : deviceFrames;
	else
		requestedBufferSize = aBuffer;

	// NOTE: SDL seems to always return a buffer size 2x of what was requested... why?
	// This doesn't match up with `pw-mon` (or Windows buffer sizes), so just request half the amount we want here.
	requestedBufferSize = (requestedBufferSize + 1) / 2; // curse integer truncation

#ifdef WINDOWS_VERSION
	// unless it was explicitly requested to be lower, clamp to min 10ms on Windows if not using WASAPI (empirically, directsound blows up with anything lower)
	if (strcasecmp(SDL_GetCurrentAudioDriver(), "wasapi") != 0)
		requestedBufferSize = std::max(((deviceRequestSpec.freq * 10) + 999) / 1000u, requestedBufferSize);
#endif

	char bufferStr[32];
	snprintf(&bufferStr[0], sizeof(bufferStr), "%u", requestedBufferSize);
	SDL_SetHintWithPriority(SDL_HINT_AUDIO_DEVICE_SAMPLE_FRAMES, &bufferStr[0], SDL_HINT_NORMAL);
	SDL_LogDebug(SDL_LOG_CATEGORY_AUDIO, "Setting buffer size hint: %u frames", requestedBufferSize);

	SDL_LogDebug(SDL_LOG_CATEGORY_AUDIO, "Requesting device: %dHz, %d channels, format %#x, buffer %u frames (driver: %s)", deviceRequestSpec.freq, deviceRequestSpec.channels,
	             deviceRequestSpec.format, requestedBufferSize, SDL_GetCurrentAudioDriver());

	// open device first
	data->deviceID = SDL_OpenAudioDevice(SDL_AUDIO_DEVICE_DEFAULT_PLAYBACK, &deviceRequestSpec);
	if (!data->deviceID)
	{
		SDL_LogError(SDL_LOG_CATEGORY_AUDIO, "Failed to open audio device: %s", SDL_GetError());
		soloud_sdl3_deinit(aSoloud);
		return UNKNOWN_ERROR;
	}
	data->deviceInitialized = true;

	// get actual device format that was negotiated
	int actualDeviceFrames = 0;
	if (!SDL_GetAudioDeviceFormat(data->deviceID, &data->deviceSpec, &actualDeviceFrames))
	{
		SDL_LogError(SDL_LOG_CATEGORY_AUDIO, "Failed to get device format: %s", SDL_GetError());
		soloud_sdl3_deinit(aSoloud);
		return UNKNOWN_ERROR;
	}

	SDL_LogDebug(SDL_LOG_CATEGORY_AUDIO, "Device opened with: %dHz, %d channels, format %#x, buffer %d frames", data->deviceSpec.freq, data->deviceSpec.channels,
	             data->deviceSpec.format, actualDeviceFrames);

	// for our stream input, use the same format as the device to avoid conversion
	data->streamSpec = data->deviceSpec;

	// create stream that matches device format exactly (no conversion)
	data->audioStream = SDL_CreateAudioStream(&data->streamSpec, &data->deviceSpec);
	if (!data->audioStream)
	{
		SDL_LogError(SDL_LOG_CATEGORY_AUDIO, "Failed to create audio stream: %s", SDL_GetError());
		soloud_sdl3_deinit(aSoloud);
		return UNKNOWN_ERROR;
	}
	data->streamInitialized = true;

	// bind stream to device
	if (!SDL_BindAudioStream(data->deviceID, data->audioStream))
	{
		SDL_LogError(SDL_LOG_CATEGORY_AUDIO, "Failed to bind stream to device: %s", SDL_GetError());
		soloud_sdl3_deinit(aSoloud);
		return UNKNOWN_ERROR;
	}

	// set up callback
	if (!SDL_SetAudioStreamGetCallback(data->audioStream, soloud_sdl3_stream_callback, data))
	{
		SDL_LogError(SDL_LOG_CATEGORY_AUDIO, "Failed to set stream callback: %s", SDL_GetError());
		soloud_sdl3_deinit(aSoloud);
		return UNKNOWN_ERROR;
	}

	data->bufferFrames = actualDeviceFrames;

	// pre-allocate mix buffer for the audio callback
	size_t maxBufferBytes = static_cast<size_t>(actualDeviceFrames) * data->streamSpec.channels * SDL_AUDIO_BYTESIZE(data->streamSpec.format);
	data->mixBuffer.resize(maxBufferBytes);

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

	// set up control functions
	aSoloud->mBackendCleanupFunc = soloud_sdl3_deinit;
	aSoloud->mBackendPauseFunc = soloud_sdl3_pause;
	aSoloud->mBackendResumeFunc = soloud_sdl3_resume;
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
