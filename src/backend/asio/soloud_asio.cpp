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

#if !defined(WITH_ASIO)

namespace SoLoud
{
result asio_init(Soloud * /*aSoloud*/, unsigned int /*aFlags*/, unsigned int /*aSamplerate*/, unsigned int /*aBufferSize*/, unsigned int /*aChannels*/,
                 const char * /*aDeviceIdentifier*/)
{
	return NOT_IMPLEMENTED;
}

result asio_enumerate_devices(Soloud * /*aSoloud*/)
{
	return NOT_IMPLEMENTED;
}
} // namespace SoLoud

#else

#include "soloud_ll_mixing.h"

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <objbase.h>

// steinberg asio sdk, interface declarations only: the driver is a plain com object and is called directly
#include "asiosys.h"
#include "asio.h"
#include "iasiodrv.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

namespace SoLoud
{
using namespace mixing; // SAMPLE_FORMAT

namespace
{
constexpr const wchar_t *ASIO_REGISTRY_PATH = L"SOFTWARE\\ASIO";
constexpr const char *IDENTIFIER_PREFIX = "asio_";
constexpr long ASIO_HOST_VERSION = 2; // asio 2.x host
constexpr double FALLBACK_SAMPLE_RATE = 44100.0;

enum LogLevel : unsigned int
{
	LOG_NONE = 0,
	LOG_ERROR,
	LOG_WARNING,
	LOG_INFO,
	LOG_DEBUG
};

struct AsioDriverEntry
{
	CLSID clsid{};
	DeviceInfo info{};
};

struct AsioData
{
	Soloud *soloud{nullptr};
	IASIO *driver{nullptr};
	DeviceInfo currentDevice{};

	// cached init parameters for device switching
	unsigned int initFlags{0};
	unsigned int requestedSampleRate{0};
	unsigned int requestedBufferSize{0};
	unsigned int requestedChannels{0};

	// negotiated stream configuration
	std::vector<ASIOBufferInfo> bufferInfos; // output channels only
	ASIOCallbacks callbacks{};
	long bufferSize{0}; // frames
	long numChannels{0};
	double sampleRate{0};
	SAMPLE_FORMAT format{SAMPLE_FLOAT32};
	unsigned int bytesPerSample{0};
	int sampleShift{0}; // right shift for 32-bit containers that hold fewer significant bits (ASIOSTInt32LSB16 etc.)
	bool outputReadySupported{false};
	std::atomic<bool> running{false}; // set around start()/stop(), checked by the buffer switch callback
	bool comInitialized{false};
	std::vector<unsigned char> interleaved; // mix() output, deinterleaved into the driver's per-channel buffers

	std::atomic<bool> deviceLost{false};
	unsigned int logLevel{LOG_NONE};
};

// asio callbacks carry no user data, so the instance has to be reachable through a global (which also limits asio to one instance per process)
std::atomic<AsioData *> gInstance{nullptr};

unsigned int parse_log_level_from_env()
{
	const char *env = getenv("SOLOUD_DEBUG");
	if (!env || !*env || *env == '0')
	{
#ifdef _DEBUG
		return LOG_WARNING;
#else
		return LOG_NONE;
#endif
	}

	if (_strnicmp(env, "debug", sizeof("debug") - 1) == 0)
		return LOG_DEBUG;
	if (_strnicmp(env, "info", sizeof("info") - 1) == 0)
		return LOG_INFO;
	if (_strnicmp(env, "warn", sizeof("warn") - 1) == 0)
		return LOG_WARNING;
	if (_strnicmp(env, "error", sizeof("error") - 1) == 0)
		return LOG_ERROR;

#ifdef _DEBUG
	return LOG_WARNING;
#else
	return LOG_NONE;
#endif
}

template <size_t N>
void wide_to_utf8(const wchar_t *aWide, std::array<char, N> &aOut)
{
	int written = WideCharToMultiByte(CP_UTF8, 0, aWide, static_cast<int>(wcslen(aWide)), aOut.data(), static_cast<int>(aOut.size()) - 1, nullptr, nullptr);
	if (written < 0)
		written = 0;
	aOut[written] = '\0';
}

template <size_t N>
bool read_registry_string(HKEY aKey, const wchar_t *aSubKey, const wchar_t *aValueName, std::array<wchar_t, N> &aOut)
{
	DWORD size = static_cast<DWORD>(aOut.size() * sizeof(wchar_t));
	if (RegGetValueW(aKey, aSubKey, aValueName, RRF_RT_REG_SZ, nullptr, aOut.data(), &size) != ERROR_SUCCESS)
		return false;
	aOut[N - 1] = L'\0';
	return aOut[0] != L'\0';
}

// lists the installed drivers from the registry without loading any of them
void enumerate_drivers(std::vector<AsioDriverEntry> &aDrivers)
{
	HKEY asioKey = nullptr;
	if (RegOpenKeyExW(HKEY_LOCAL_MACHINE, ASIO_REGISTRY_PATH, 0, KEY_READ, &asioKey) != ERROR_SUCCESS)
		return;

	for (DWORD index = 0;; index++)
	{
		std::array<wchar_t, 256> keyName{};
		DWORD keyNameLength = static_cast<DWORD>(keyName.size());
		if (RegEnumKeyExW(asioKey, index, keyName.data(), &keyNameLength, nullptr, nullptr, nullptr, nullptr) != ERROR_SUCCESS)
			break;

		std::array<wchar_t, 64> clsidString{};
		AsioDriverEntry entry{};
		if (!read_registry_string(asioKey, keyName.data(), L"CLSID", clsidString) || FAILED(CLSIDFromString(clsidString.data(), &entry.clsid)))
			continue;

		// the description is the user-facing name, the key name is what the sdk's own host code shows
		std::array<wchar_t, 256> description{};
		if (read_registry_string(asioKey, keyName.data(), L"Description", description))
			wide_to_utf8(description.data(), entry.info.name);
		else
			wide_to_utf8(keyName.data(), entry.info.name);

		std::array<char, 64> clsidUtf8{};
		wide_to_utf8(clsidString.data(), clsidUtf8);
		snprintf(entry.info.identifier.data(), entry.info.identifier.size(), "%s%s", IDENTIFIER_PREFIX, clsidUtf8.data());

		// asio has no notion of a default device, the first driver stands in for it
		entry.info.isDefault = aDrivers.empty();
		entry.info.isExclusive = true;
		entry.info.nativeDeviceInfo = nullptr;
		aDrivers.push_back(entry);
	}

	RegCloseKey(asioKey);
}

// resolves an identifier from enumerate_drivers to its driver, NULL or empty picks the default
bool find_driver(const char *aIdentifier, AsioDriverEntry &aEntry)
{
	std::vector<AsioDriverEntry> drivers;
	enumerate_drivers(drivers);
	for (const AsioDriverEntry &entry : drivers)
	{
		if (!aIdentifier || !*aIdentifier || strcmp(entry.info.identifier.data(), aIdentifier) == 0)
		{
			aEntry = entry;
			return true;
		}
	}
	return false;
}

// maps an asio sample type to the soloud output format that can be written to it directly
bool map_sample_type(ASIOSampleType aType, SAMPLE_FORMAT &aFormat, unsigned int &aBytesPerSample, int &aShift)
{
	aShift = 0;
	switch (aType)
	{
	case ASIOSTInt16LSB:
		aFormat = SAMPLE_SIGNED16;
		aBytesPerSample = 2;
		return true;
	case ASIOSTInt24LSB:
		aFormat = SAMPLE_SIGNED24;
		aBytesPerSample = 3;
		return true;
	case ASIOSTInt32LSB:
		aFormat = SAMPLE_SIGNED32;
		aBytesPerSample = 4;
		return true;
	case ASIOSTFloat32LSB:
		aFormat = SAMPLE_FLOAT32;
		aBytesPerSample = 4;
		return true;
	// 32-bit containers with the significant bits right-aligned
	case ASIOSTInt32LSB16:
	case ASIOSTInt32LSB18:
	case ASIOSTInt32LSB20:
	case ASIOSTInt32LSB24:
		aFormat = SAMPLE_SIGNED32;
		aBytesPerSample = 4;
		aShift = aType == ASIOSTInt32LSB16 ? 16 : aType == ASIOSTInt32LSB18 ? 14 : aType == ASIOSTInt32LSB20 ? 12 : 8;
		return true;
	default:
		// big-endian, 64-bit float and dsd formats
		return false;
	}
}

// clamps a requested buffer size to the sizes the driver reports as valid, AUTO picks the driver's preferred size
long clamp_buffer_size(long aMinSize, long aMaxSize, long aPreferredSize, long aGranularity, unsigned int aRequested)
{
	if (aRequested == Soloud::AUTO)
		return aPreferredSize;

	long size = static_cast<long>(aRequested);
	if (size <= aMinSize)
		return aMinSize;
	if (size >= aMaxSize)
		return aMaxSize;

	if (aGranularity > 0)
		return aMinSize + ((size - aMinSize) / aGranularity) * aGranularity;

	if (aGranularity == -1)
	{
		// powers of two counted from the minimum
		long candidate = aMinSize > 0 ? aMinSize : 1;
		while (candidate * 2 <= size)
			candidate *= 2;
		return candidate;
	}

	return size;
}

void deinterleave_channel(const AsioData *data, long aChannel, void *aDest)
{
	const unsigned char *src = data->interleaved.data() + static_cast<size_t>(aChannel) * data->bytesPerSample;
	const size_t stride = static_cast<size_t>(data->numChannels) * data->bytesPerSample;

	if (data->sampleShift != 0)
	{
		int *dst = static_cast<int *>(aDest);
		for (long i = 0; i < data->bufferSize; i++)
		{
			int sample;
			memcpy(&sample, src + i * stride, sizeof(sample));
			dst[i] = sample >> data->sampleShift;
		}
		return;
	}

	unsigned char *dst = static_cast<unsigned char *>(aDest);
	for (long i = 0; i < data->bufferSize; i++)
		memcpy(dst + i * data->bytesPerSample, src + i * stride, data->bytesPerSample);
}

void asio_buffer_switch(long aBufferIndex, ASIOBool /*aDirectProcess*/)
{
	AsioData *data = gInstance.load(std::memory_order_acquire);
	if (!data || !data->running.load(std::memory_order_acquire))
		return;

	data->soloud->mix(data->interleaved.data(), static_cast<unsigned int>(data->bufferSize), data->format);
	for (long channel = 0; channel < data->numChannels; channel++)
		deinterleave_channel(data, channel, data->bufferInfos[channel].buffers[aBufferIndex]);

	if (data->outputReadySupported)
		data->driver->outputReady();
}

ASIOTime *asio_buffer_switch_time_info(ASIOTime * /*aTimeInfo*/, long aBufferIndex, ASIOBool aDirectProcess)
{
	asio_buffer_switch(aBufferIndex, aDirectProcess);
	return nullptr;
}

void asio_sample_rate_did_change(ASIOSampleRate aSampleRate)
{
	// the stream has to be rebuilt around the new rate, so let the application reopen the device
	AsioData *data = gInstance.load(std::memory_order_acquire);
	if (data && static_cast<double>(aSampleRate) != data->sampleRate)
		data->deviceLost.store(true);
}

long asio_message(long aSelector, long aValue, void * /*aMessage*/, double * /*aOpt*/)
{
	AsioData *data = gInstance.load(std::memory_order_acquire);
	switch (aSelector)
	{
	case kAsioSelectorSupported:
		return (aValue == kAsioEngineVersion || aValue == kAsioResetRequest || aValue == kAsioResyncRequest || aValue == kAsioLatenciesChanged ||
		        aValue == kAsioSupportsTimeInfo)
		           ? 1
		           : 0;
	case kAsioEngineVersion:
		return ASIO_HOST_VERSION;
	case kAsioResetRequest:
		// the driver wants to be closed and reopened, e.g. after a buffer size change in its control panel
		if (data)
			data->deviceLost.store(true);
		return 1;
	case kAsioResyncRequest:
	case kAsioLatenciesChanged:
		// nothing is cached that depends on timestamps or latencies
		return 1;
	case kAsioSupportsTimeInfo:
		return 1;
	default:
		return 0;
	}
}

void stop_driver(AsioData *data)
{
	if (data->running.load())
	{
		data->running.store(false);
		data->driver->stop(); // no callbacks arrive after this returns
	}
}

void close_driver(AsioData *data)
{
	if (!data->driver)
		return;

	stop_driver(data);
	if (!data->bufferInfos.empty())
	{
		data->driver->disposeBuffers();
		data->bufferInfos.clear();
	}
	data->driver->Release();
	data->driver = nullptr;
}

result start_driver(AsioData *data)
{
	// whatever was left in the buffers from before a stop would otherwise play for one period
	for (const ASIOBufferInfo &info : data->bufferInfos)
		for (void *half : info.buffers)
			memset(half, 0, static_cast<size_t>(data->bufferSize) * data->bytesPerSample);

	data->running.store(true, std::memory_order_release);
	if (data->driver->start() != ASE_OK)
	{
		data->running.store(false);
		if (data->logLevel >= LOG_ERROR)
			SoLoud::logStdout("[ASIO ERROR] Failed to start '%s'\n", data->currentDevice.name.data());
		return UNKNOWN_ERROR;
	}

	return SO_NO_ERROR;
}

// loads the driver and negotiates the stream configuration from the cached init parameters, but doesn't start it
result open_driver(AsioData *data, const AsioDriverEntry &aEntry)
{
	const char *name = aEntry.info.name.data();

	IASIO *driver = nullptr;
	HRESULT hr = CoCreateInstance(aEntry.clsid, nullptr, CLSCTX_INPROC_SERVER, aEntry.clsid, reinterpret_cast<void **>(&driver));
	if (FAILED(hr) || !driver)
	{
		if (data->logLevel >= LOG_ERROR)
			SoLoud::logStdout("[ASIO ERROR] Failed to load driver '%s' (0x%08lx)\n", name, static_cast<unsigned long>(hr));
		return UNKNOWN_ERROR;
	}
	data->driver = driver;

	// the window handle only serves as the parent of the driver's control panel
	if (!driver->init(GetDesktopWindow()))
	{
		std::array<char, 128> message{}; // ASIODriverInfo::errorMessage size
		driver->getErrorMessage(message.data());
		if (data->logLevel >= LOG_ERROR)
			SoLoud::logStdout("[ASIO ERROR] Driver '%s' failed to initialize: %s\n", name, message.data());
		close_driver(data);
		return UNKNOWN_ERROR;
	}

	long numInputs = 0, numOutputs = 0;
	if (driver->getChannels(&numInputs, &numOutputs) != ASE_OK || numOutputs < 1)
	{
		if (data->logLevel >= LOG_ERROR)
			SoLoud::logStdout("[ASIO ERROR] Driver '%s' has no output channels\n", name);
		close_driver(data);
		return UNKNOWN_ERROR;
	}
	data->numChannels = std::min(static_cast<long>(data->requestedChannels), numOutputs);

	long minSize = 0, maxSize = 0, preferredSize = 0, granularity = 0;
	if (driver->getBufferSize(&minSize, &maxSize, &preferredSize, &granularity) != ASE_OK)
	{
		if (data->logLevel >= LOG_ERROR)
			SoLoud::logStdout("[ASIO ERROR] Driver '%s' didn't report its buffer sizes\n", name);
		close_driver(data);
		return UNKNOWN_ERROR;
	}
	data->bufferSize = clamp_buffer_size(minSize, maxSize, preferredSize, granularity, data->requestedBufferSize);

	// only touch the driver's clock when a specific rate was asked for and the driver accepts it, otherwise keep what the user configured
	if (data->requestedSampleRate != Soloud::AUTO && driver->canSampleRate(static_cast<ASIOSampleRate>(data->requestedSampleRate)) == ASE_OK)
		driver->setSampleRate(static_cast<ASIOSampleRate>(data->requestedSampleRate));

	ASIOSampleRate currentRate = 0;
	if (driver->getSampleRate(&currentRate) != ASE_OK || currentRate <= 0)
	{
		// some drivers report no rate until one has been set
		driver->setSampleRate(FALLBACK_SAMPLE_RATE);
		if (driver->getSampleRate(&currentRate) != ASE_OK || currentRate <= 0)
		{
			if (data->logLevel >= LOG_ERROR)
				SoLoud::logStdout("[ASIO ERROR] Driver '%s' didn't report a sample rate\n", name);
			close_driver(data);
			return UNKNOWN_ERROR;
		}
	}
	data->sampleRate = currentRate;

	// every output we use has to share one sample type
	for (long channel = 0; channel < data->numChannels; channel++)
	{
		ASIOChannelInfo channelInfo{};
		channelInfo.channel = channel;
		channelInfo.isInput = ASIOFalse;

		SAMPLE_FORMAT format = SAMPLE_FLOAT32;
		unsigned int bytesPerSample = 0;
		int shift = 0;
		if (driver->getChannelInfo(&channelInfo) != ASE_OK || !map_sample_type(channelInfo.type, format, bytesPerSample, shift) ||
		    (channel > 0 && (format != data->format || shift != data->sampleShift)))
		{
			if (data->logLevel >= LOG_ERROR)
				SoLoud::logStdout("[ASIO ERROR] Driver '%s' uses an unsupported sample type (%ld) on output %ld\n", name, static_cast<long>(channelInfo.type),
				                  channel);
			close_driver(data);
			return UNKNOWN_ERROR;
		}
		data->format = format;
		data->bytesPerSample = bytesPerSample;
		data->sampleShift = shift;
	}

	data->callbacks.bufferSwitch = asio_buffer_switch;
	data->callbacks.sampleRateDidChange = asio_sample_rate_did_change;
	data->callbacks.asioMessage = asio_message;
	data->callbacks.bufferSwitchTimeInfo = asio_buffer_switch_time_info;

	data->bufferInfos.assign(static_cast<size_t>(data->numChannels), ASIOBufferInfo{});
	for (long channel = 0; channel < data->numChannels; channel++)
	{
		data->bufferInfos[channel].isInput = ASIOFalse;
		data->bufferInfos[channel].channelNum = channel;
	}

	ASIOError err = driver->createBuffers(data->bufferInfos.data(), data->numChannels, data->bufferSize, &data->callbacks);
	if (err != ASE_OK && data->bufferSize != preferredSize)
	{
		// some drivers only accept the size configured in their control panel
		data->bufferSize = preferredSize;
		err = driver->createBuffers(data->bufferInfos.data(), data->numChannels, data->bufferSize, &data->callbacks);
	}
	if (err != ASE_OK)
	{
		data->bufferInfos.clear();
		if (data->logLevel >= LOG_ERROR)
			SoLoud::logStdout("[ASIO ERROR] Driver '%s' failed to create buffers (%ld)\n", name, static_cast<long>(err));
		close_driver(data);
		return UNKNOWN_ERROR;
	}

	// drivers that support it wait for outputReady() after each buffer switch instead of adding a period of latency
	data->outputReadySupported = (driver->outputReady() == ASE_OK);
	data->interleaved.assign(static_cast<size_t>(data->bufferSize) * data->numChannels * data->bytesPerSample, 0);

	data->currentDevice = aEntry.info;
	data->deviceLost.store(false);

	if (data->logLevel >= LOG_INFO)
	{
		long inputLatency = 0, outputLatency = 0;
		driver->getLatencies(&inputLatency, &outputLatency);
		SoLoud::logStdout("[ASIO INFO] Opened '%s': %u Hz, %ld/%ld outputs, %ld frames (%ld-%ld, preferred %ld), %ld frames output latency\n", name,
		                  static_cast<unsigned int>(data->sampleRate), data->numChannels, numOutputs, data->bufferSize, minSize, maxSize, preferredSize,
		                  outputLatency);
	}

	return SO_NO_ERROR;
}

void asio_deinit(Soloud *aSoloud)
{
	AsioData *data = static_cast<AsioData *>(aSoloud->mBackendData);
	if (!data)
		return;

	close_driver(data);
	if (data->comInitialized)
		CoUninitialize();

	gInstance.store(nullptr);
	delete data;
	aSoloud->mBackendData = nullptr;
}

result asio_pause(Soloud *aSoloud)
{
	AsioData *data = static_cast<AsioData *>(aSoloud->mBackendData);
	if (!data || !data->driver)
		return INVALID_PARAMETER;

	stop_driver(data);
	return SO_NO_ERROR;
}

result asio_resume(Soloud *aSoloud)
{
	AsioData *data = static_cast<AsioData *>(aSoloud->mBackendData);
	if (!data || !data->driver)
		return INVALID_PARAMETER;

	if (data->running.load())
		return SO_NO_ERROR;
	return start_driver(data);
}

result asio_get_current_device(Soloud *aSoloud, DeviceInfo *pDeviceInfo)
{
	AsioData *data = static_cast<AsioData *>(aSoloud->mBackendData);
	if (!data || !data->driver)
		return INVALID_PARAMETER;

	*pDeviceInfo = data->currentDevice;
	return SO_NO_ERROR;
}

result asio_set_device(Soloud *aSoloud, const char *aDeviceIdentifier)
{
	AsioData *data = static_cast<AsioData *>(aSoloud->mBackendData);
	if (!data)
		return INVALID_PARAMETER;

	AsioDriverEntry entry{};
	if (!find_driver(aDeviceIdentifier, entry))
		return INVALID_PARAMETER;

	// already playing on that driver (a lost device is reopened instead, that's the recovery path)
	if (data->driver && !data->deviceLost.load() && strcmp(entry.info.identifier.data(), data->currentDevice.identifier.data()) == 0)
		return SO_NO_ERROR;

	const double oldSampleRate = data->sampleRate;
	const long oldBufferSize = data->bufferSize;
	const long oldChannels = data->numChannels;

	close_driver(data);
	result res = open_driver(data, entry);
	if (res != SO_NO_ERROR)
	{
		// the previous device is gone at this point, the application has to pick another one
		data->deviceLost.store(true);
		return res;
	}

	if (data->sampleRate != oldSampleRate || data->bufferSize != oldBufferSize || data->numChannels != oldChannels)
	{
		aSoloud->postinit_internal(static_cast<unsigned int>(data->sampleRate), static_cast<unsigned int>(data->bufferSize), data->initFlags,
		                           static_cast<unsigned int>(data->numChannels));
	}

	return start_driver(data);
}

result asio_get_device_latency(Soloud *aSoloud, unsigned int *pLatencyFrames)
{
	AsioData *data = static_cast<AsioData *>(aSoloud->mBackendData);
	if (!data || !data->driver)
		return INVALID_PARAMETER;

	long inputLatency = 0, outputLatency = 0;
	if (data->driver->getLatencies(&inputLatency, &outputLatency) != ASE_OK)
		return UNKNOWN_ERROR;

	*pLatencyFrames = outputLatency > 0 ? static_cast<unsigned int>(outputLatency) : 0;
	return SO_NO_ERROR;
}

result asio_get_buffer_size_limits(Soloud *aSoloud, unsigned int *pMinSize, unsigned int *pMaxSize, unsigned int *pPreferredSize, int *pGranularity)
{
	AsioData *data = static_cast<AsioData *>(aSoloud->mBackendData);
	if (!data || !data->driver)
		return INVALID_PARAMETER;

	// queried live, since the control panel can change these at any time
	long minSize = 0, maxSize = 0, preferredSize = 0, granularity = 0;
	if (data->driver->getBufferSize(&minSize, &maxSize, &preferredSize, &granularity) != ASE_OK)
		return UNKNOWN_ERROR;

	*pMinSize = minSize > 0 ? static_cast<unsigned int>(minSize) : 0;
	*pMaxSize = maxSize > 0 ? static_cast<unsigned int>(maxSize) : 0;
	*pPreferredSize = preferredSize > 0 ? static_cast<unsigned int>(preferredSize) : 0;
	*pGranularity = static_cast<int>(granularity);
	return SO_NO_ERROR;
}

result asio_open_control_panel(Soloud *aSoloud)
{
	AsioData *data = static_cast<AsioData *>(aSoloud->mBackendData);
	if (!data || !data->driver)
		return INVALID_PARAMETER;

	ASIOError err = data->driver->controlPanel();
	if (err == ASE_NotPresent)
		return NOT_IMPLEMENTED;
	return err == ASE_OK ? SO_NO_ERROR : UNKNOWN_ERROR;
}

bool asio_is_device_lost(Soloud *aSoloud)
{
	AsioData *data = static_cast<AsioData *>(aSoloud->mBackendData);
	return data && data->deviceLost.load();
}
} // namespace

result asio_enumerate_devices(Soloud *aSoloud)
{
	if (!aSoloud)
		return INVALID_PARAMETER;

	std::vector<AsioDriverEntry> drivers;
	enumerate_drivers(drivers);

	aSoloud->mDeviceList = new DeviceInfo[drivers.size()];
	aSoloud->mDeviceCount = static_cast<unsigned int>(drivers.size());
	for (size_t i = 0; i < drivers.size(); i++)
		aSoloud->mDeviceList[i] = drivers[i].info;

	return SO_NO_ERROR;
}

result asio_init(Soloud *aSoloud, unsigned int aFlags, unsigned int aSamplerate, unsigned int aBufferSize, unsigned int aChannels, const char *aDeviceIdentifier)
{
	AsioData *data = new AsioData();
	data->logLevel = parse_log_level_from_env();

	AsioData *expected = nullptr;
	if (!gInstance.compare_exchange_strong(expected, data))
	{
		if (data->logLevel >= LOG_ERROR)
			SoLoud::logStdout("[ASIO ERROR] Another SoLoud instance is already using ASIO\n");
		delete data;
		return UNKNOWN_ERROR;
	}

	aSoloud->mBackendData = data;
	data->soloud = aSoloud;

	// cache initialization parameters for device switching
	data->initFlags = aFlags;
	data->requestedSampleRate = aSamplerate;
	data->requestedBufferSize = aBufferSize;
	data->requestedChannels = aChannels;

	// asio drivers are apartment-threaded com objects, every call into them has to come from the thread that loaded them (this one)
	HRESULT hr = CoInitializeEx(nullptr, COINIT_APARTMENTTHREADED);
	data->comInitialized = SUCCEEDED(hr);
	if (hr == RPC_E_CHANGED_MODE && data->logLevel >= LOG_WARNING)
		SoLoud::logStdout("[ASIO WARNING] Calling thread already uses multithreaded COM, drivers may misbehave\n");

	AsioDriverEntry entry{};
	if (!find_driver(aDeviceIdentifier, entry))
	{
		if (data->logLevel >= LOG_ERROR)
			SoLoud::logStdout(aDeviceIdentifier && *aDeviceIdentifier ? "[ASIO ERROR] No driver matches '%s'\n" : "[ASIO ERROR] No drivers installed\n",
			                  aDeviceIdentifier);
		asio_deinit(aSoloud);
		return aDeviceIdentifier && *aDeviceIdentifier ? INVALID_PARAMETER : UNKNOWN_ERROR;
	}

	result res = open_driver(data, entry);
	if (res != SO_NO_ERROR)
	{
		asio_deinit(aSoloud);
		return res;
	}

	aSoloud->postinit_internal(static_cast<unsigned int>(data->sampleRate), static_cast<unsigned int>(data->bufferSize), aFlags,
	                           static_cast<unsigned int>(data->numChannels));

	res = start_driver(data);
	if (res != SO_NO_ERROR)
	{
		asio_deinit(aSoloud);
		return res;
	}

	aSoloud->mBackendCleanupFunc = asio_deinit;
	aSoloud->mBackendPauseFunc = asio_pause;
	aSoloud->mBackendResumeFunc = asio_resume;

	aSoloud->mEnumerateDevicesFunc = asio_enumerate_devices;
	aSoloud->mGetCurrentDeviceFunc = asio_get_current_device;
	aSoloud->mSetDeviceFunc = asio_set_device;
	aSoloud->mGetDeviceLatencyFunc = asio_get_device_latency;
	aSoloud->mGetBufferSizeLimitsFunc = asio_get_buffer_size_limits;
	aSoloud->mOpenControlPanelFunc = asio_open_control_panel;
	aSoloud->mIsDeviceLostFunc = asio_is_device_lost;

	aSoloud->mBackendString = "ASIO";
	return SO_NO_ERROR;
}
} // namespace SoLoud

#endif // WITH_ASIO
