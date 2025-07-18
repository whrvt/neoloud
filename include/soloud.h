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

#ifndef SOLOUD_H
#define SOLOUD_H

#include <math.h>   // sin
#include <stdlib.h> // rand

#include "soloud_config.h"

#include "soloud_audiosource3d.h"
#include "soloud_fader.h"
#include "soloud_filter.h"
#include "soloud_intrin.h"

#include "soloud_error.h"

namespace SoLoud
{
	class AudioSource;
	class AudioSourceInstance;
	class Bus;
	class BusInstance;
	class Queue;
	class QueueInstance;
	class AudioSourceInstance3dData;

	// Generic device information structure for cross-backend compatibility
	struct DeviceInfo
	{
		char name[256];         // Human-readable device name
		char identifier[128];   // Backend-specific device identifier
		bool isDefault;         // Whether this is the default device
		void *nativeDeviceInfo; // Backend-specific device info (optional)
	};

	// Device management function pointer types
	typedef result (*enumerateDevicesFunc)(Soloud *aSoloud);
	typedef result (*getCurrentDeviceFunc)(Soloud *aSoloud, DeviceInfo *pDeviceInfo);
	typedef result (*setDeviceFunc)(Soloud *aSoloud, const char *deviceIdentifier);

	// Soloud core class.
	class Soloud
	{
	public:
		// Back-end data; content is up to the back-end implementation.
		void *mBackendData;
		// Pointer for the audio thread mutex.
		void *mAudioThreadMutex;
		// Flag for when we're inside the mutex, used for debugging.
		bool mInsideAudioThreadMutex;
		// Called by SoLoud to shut down the back-end. If NULL, not called. Should be set by back-end.
		soloudCallFunction mBackendCleanupFunc;

		// Some backends like CoreAudio on iOS must be paused/resumed in some cases. On incoming call as instance.
		soloudResultFunction mBackendPauseFunc;
		soloudResultFunction mBackendResumeFunc;

		// Device management function pointers (set by backend during initialization)
		enumerateDevicesFunc mEnumerateDevicesFunc;
		getCurrentDeviceFunc mGetCurrentDeviceFunc;
		setDeviceFunc mSetDeviceFunc;

		// Internal device list storage
		DeviceInfo *mDeviceList;
		unsigned int mDeviceCount;

		// CTor
		Soloud();
		// DTor
		~Soloud();

		enum BACKENDS
		{
			AUTO = 0,
			MINIAUDIO,
			SDL3,
			NOSOUND,
			NULLDRIVER,
			BACKEND_MAX,
		};

		enum FLAGS
		{
			// Use round-off clipper
			CLIP_ROUNDOFF = 1,
			ENABLE_VISUALIZATION = 2,
			LEFT_HANDED_3D = 4,
			NO_FPU_REGISTER_CHANGE = 8
		};

		enum WAVEFORM
		{
			WAVE_SQUARE = 0,
			WAVE_SAW,
			WAVE_SIN,
			WAVE_TRIANGLE,
			WAVE_BOUNCE,
			WAVE_JAWS,
			WAVE_HUMPS,
			WAVE_FSQUARE,
			WAVE_FSAW
		};

		enum RESAMPLER
		{
			RESAMPLER_POINT,
			RESAMPLER_LINEAR,
			RESAMPLER_CATMULLROM
		};

		// Initialize SoLoud. Must be called before SoLoud can be used.
		result init(unsigned int aFlags = Soloud::CLIP_ROUNDOFF, unsigned int aBackend = Soloud::AUTO, unsigned int aSamplerate = Soloud::AUTO,
		            unsigned int aBufferSize = Soloud::AUTO, unsigned int aChannels = 2);

		result pause();
		result resume();

		// Deinitialize SoLoud. Must be called before shutting down.
		void deinit();

		// Query SoLoud version number (should equal to SOLOUD_VERSION macro)
		unsigned int getVersion() const;

		// Translate error number to an asciiz string
		const char *getErrorString(result aErrorCode) const;

		// Returns current backend ID (BACKENDS enum)
		unsigned int getBackendId();
		// Returns current backend string. May be NULL.
		const char *getBackendString();
		// Returns current backend channel count (1 mono, 2 stereo, etc)
		unsigned int getBackendChannels();
		// Returns current backend sample rate
		unsigned int getBackendSamplerate();
		// Returns current backend buffer size
		unsigned int getBackendBufferSize();

		// Device management API
		/*
		Enumerates available playback devices.

		Parameters
		----------
		ppDevices (out)
		    Pointer to array of DeviceInfo structures. Points to internal storage that remains
		    valid until the next call to enumerateDevices() or until SoLoud is destroyed.

		pDeviceCount (out)
		    Number of devices returned in ppDevices array.

		Return Value
		------------
		SO_NO_ERROR if successful; error code otherwise.
		NOT_IMPLEMENTED if current backend doesn't support device enumeration.

		Thread Safety
		-------------
		Unsafe. Do not call from multiple threads simultaneously.

		Remarks
		-------
		The returned device list points to internal storage and becomes invalid on the next
		call to enumerateDevices(). Device availability may change over time (devices can
		be plugged/unplugged), so re-enumerate when needed.
		*/
		result enumerateDevices(DeviceInfo * *ppDevices, unsigned int *pDeviceCount);

		/*
		Gets information about the currently active playback device.

		Parameters
		----------
		pDeviceInfo (out)
		    DeviceInfo structure to receive current device information.

		Return Value
		------------
		SO_NO_ERROR if successful; error code otherwise.
		NOT_IMPLEMENTED if current backend doesn't support device enumeration.

		Thread Safety
		-------------
		Safe.
		*/
		result getCurrentDevice(DeviceInfo * pDeviceInfo);

		/*
		Switches to a different playback device.

		Parameters
		----------
		deviceIdentifier (in)
		    Device identifier string from DeviceInfo.identifier, or NULL for default device.

		Return Value
		------------
		SO_NO_ERROR if successful; error code otherwise.
		NOT_IMPLEMENTED if current backend doesn't support device switching.
		INVALID_PARAMETER if device identifier is invalid.

		Thread Safety
		-------------
		Safe.

		Remarks
		-------
		Audio playback will continue seamlessly on the new device. If the new device has
		different capabilities (sample rate, channels), internal buffers will be reconfigured
		automatically. Currently playing sounds will continue from their current position.
		*/
		result setDevice(const char *deviceIdentifier);

		// Set speaker position in 3d space
		result setSpeakerPosition(unsigned int aChannel, float aX, float aY, float aZ);
		// Get speaker position in 3d space
		result getSpeakerPosition(unsigned int aChannel, float &aX, float &aY, float &aZ);

		// Start playing a sound. Returns voice handle, which can be ignored or used to alter the playing sound's parameters. Negative volume means to use default.
		handle play(AudioSource & aSound, float aVolume = -1.0f, float aPan = 0.0f, bool aPaused = 0, unsigned int aBus = 0);
		// Start playing a sound delayed in relation to other sounds called via this function. Negative volume means to use default.
		handle playClocked(time aSoundTime, AudioSource & aSound, float aVolume = -1.0f, float aPan = 0.0f, unsigned int aBus = 0);
		// Start playing a 3d audio source
		handle play3d(AudioSource & aSound, float aPosX, float aPosY, float aPosZ, float aVelX = 0.0f, float aVelY = 0.0f, float aVelZ = 0.0f, float aVolume = 1.0f,
		              bool aPaused = 0, unsigned int aBus = 0);
		// Start playing a 3d audio source, delayed in relation to other sounds called via this function.
		handle play3dClocked(time aSoundTime, AudioSource & aSound, float aPosX, float aPosY, float aPosZ, float aVelX = 0.0f, float aVelY = 0.0f,
		                     float aVelZ = 0.0f, float aVolume = 1.0f, unsigned int aBus = 0);
		// Start playing a sound without any panning. It will be played at full volume.
		handle playBackground(AudioSource & aSound, float aVolume = -1.0f, bool aPaused = 0, unsigned int aBus = 0);

		// Seek the audio stream to certain point in time. Some streams can't seek backwards. Relative play speed affects time.
		result seek(handle aVoiceHandle, time aSeconds);
		// Stop the sound.
		void stop(handle aVoiceHandle);
		// Stop all voices.
		void stopAll();
		// Stop all voices that play this sound source
		void stopAudioSource(AudioSource & aSound);
		// Count voices that play this audio source
		int countAudioSource(AudioSource & aSound);

		// Set a live filter parameter. Use 0 for the global filters.
		void setFilterParameter(handle aVoiceHandle, unsigned int aFilterId, unsigned int aAttributeId, float aValue);
		// Get a live filter parameter. Use 0 for the global filters.
		float getFilterParameter(handle aVoiceHandle, unsigned int aFilterId, unsigned int aAttributeId);
		// Fade a live filter parameter. Use 0 for the global filters.
		void fadeFilterParameter(handle aVoiceHandle, unsigned int aFilterId, unsigned int aAttributeId, float aTo, time aTime);
		// Oscillate a live filter parameter. Use 0 for the global filters.
		void oscillateFilterParameter(handle aVoiceHandle, unsigned int aFilterId, unsigned int aAttributeId, float aFrom, float aTo, time aTime);

		// Get current play time, in seconds.
		time getStreamTime(handle aVoiceHandle);
		// Get current sample position, in seconds.
		time getStreamPosition(handle aVoiceHandle);
		// Get current pause state.
		bool getPause(handle aVoiceHandle);
		// Get current volume.
		float getVolume(handle aVoiceHandle);
		// Get current overall volume (set volume * 3d volume)
		float getOverallVolume(handle aVoiceHandle);
		// Get current pan.
		float getPan(handle aVoiceHandle);
		// Get current sample rate.
		float getSamplerate(handle aVoiceHandle);
		// Get current voice protection state.
		bool getProtectVoice(handle aVoiceHandle);
		// Get the current number of busy voices.
		unsigned int getActiveVoiceCount();
		// Get the current number of voices in SoLoud
		unsigned int getVoiceCount();
		// Check if the handle is still valid, or if the sound has stopped.
		bool isValidVoiceHandle(handle aVoiceHandle);
		// Get current relative play speed.
		float getRelativePlaySpeed(handle aVoiceHandle);
		// Get current post-clip scaler value.
		float getPostClipScaler() const;
		// Get the current main resampler
		unsigned int getMainResampler() const;
		// Get current global volume
		float getGlobalVolume() const;
		// Get current maximum active voice setting
		unsigned int getMaxActiveVoiceCount() const;
		// Query whether a voice is set to loop.
		bool getLooping(handle aVoiceHandle);
		// Query whether a voice is set to auto-stop when it ends.
		bool getAutoStop(handle aVoiceHandle);
		// Get voice loop point value
		time getLoopPoint(handle aVoiceHandle);

		// Set voice loop point value
		void setLoopPoint(handle aVoiceHandle, time aLoopPoint);
		// Set voice's loop state
		void setLooping(handle aVoiceHandle, bool aLooping);
		// Set whether sound should auto-stop when it ends
		void setAutoStop(handle aVoiceHandle, bool aAutoStop);
		// Set current maximum active voice setting
		result setMaxActiveVoiceCount(unsigned int aVoiceCount);
		// Set behavior for inaudible sounds
		void setInaudibleBehavior(handle aVoiceHandle, bool aMustTick, bool aKill);
		// Set the global volume
		void setGlobalVolume(float aVolume);
		// Set the post clip scaler value
		void setPostClipScaler(float aScaler);
		// Set the main resampler
		void setMainResampler(unsigned int aResampler);
		// Set the pause state
		void setPause(handle aVoiceHandle, bool aPause);
		// Pause all voices
		void setPauseAll(bool aPause);
		// Set the relative play speed
		result setRelativePlaySpeed(handle aVoiceHandle, float aSpeed);
		// Set the voice protection state
		void setProtectVoice(handle aVoiceHandle, bool aProtect);
		// Set the sample rate
		void setSamplerate(handle aVoiceHandle, float aSamplerate);
		// Set panning value; -1 is left, 0 is center, 1 is right
		void setPan(handle aVoiceHandle, float aPan);
		// Set absolute left/right volumes
		void setPanAbsolute(handle aVoiceHandle, float aLVolume, float aRVolume);
		// Set channel volume (volume for a specific speaker)
		void setChannelVolume(handle aVoiceHandle, unsigned int aChannel, float aVolume);
		// Set overall volume
		void setVolume(handle aVoiceHandle, float aVolume);
		// Set delay, in samples, before starting to play samples. Calling this on a live sound will cause glitches.
		void setDelaySamples(handle aVoiceHandle, unsigned int aSamples);

		// Set up volume fader
		void fadeVolume(handle aVoiceHandle, float aTo, time aTime);
		// Set up panning fader
		void fadePan(handle aVoiceHandle, float aTo, time aTime);
		// Set up relative play speed fader
		void fadeRelativePlaySpeed(handle aVoiceHandle, float aTo, time aTime);
		// Set up global volume fader
		void fadeGlobalVolume(float aTo, time aTime);
		// Schedule a stream to pause
		void schedulePause(handle aVoiceHandle, time aTime);
		// Schedule a stream to stop
		void scheduleStop(handle aVoiceHandle, time aTime);

		// Set up volume oscillator
		void oscillateVolume(handle aVoiceHandle, float aFrom, float aTo, time aTime);
		// Set up panning oscillator
		void oscillatePan(handle aVoiceHandle, float aFrom, float aTo, time aTime);
		// Set up relative play speed oscillator
		void oscillateRelativePlaySpeed(handle aVoiceHandle, float aFrom, float aTo, time aTime);
		// Set up global volume oscillator
		void oscillateGlobalVolume(float aFrom, float aTo, time aTime);

		// Set global filters. Set to NULL to clear the filter.
		void setGlobalFilter(unsigned int aFilterId, Filter *aFilter);

		// Enable or disable visualization data gathering
		void setVisualizationEnable(bool aEnable);

		// Calculate and get 256 floats of FFT data for visualization. Visualization has to be enabled before use.
		float *calcFFT();

		// Get 256 floats of wave data for visualization. Visualization has to be enabled before use.
		float *getWave();

		// Get approximate output volume for a channel for visualization. Visualization has to be enabled before use.
		float getApproximateVolume(unsigned int aChannel);

		// Get current loop count. Returns 0 if handle is not valid. (All audio sources may not update loop count)
		unsigned int getLoopCount(handle aVoiceHandle);

		// Get audiosource-specific information from a voice.
		float getInfo(handle aVoiceHandle, unsigned int aInfoKey);

		// Create a voice group. Returns 0 if unable (out of voice groups / out of memory)
		handle createVoiceGroup();
		// Destroy a voice group.
		result destroyVoiceGroup(handle aVoiceGroupHandle);
		// Add a voice handle to a voice group
		result addVoiceToGroup(handle aVoiceGroupHandle, handle aVoiceHandle);
		// Is this handle a valid voice group?
		bool isVoiceGroup(handle aVoiceGroupHandle);
		// Is this voice group empty?
		bool isVoiceGroupEmpty(handle aVoiceGroupHandle);

		// Perform 3d audio parameter update
		void update3dAudio();

		// Set the speed of sound constant for doppler
		result set3dSoundSpeed(float aSpeed);
		// Get the current speed of sound constant for doppler
		float get3dSoundSpeed();
		// Set 3d listener parameters
		void set3dListenerParameters(float aPosX, float aPosY, float aPosZ, float aAtX, float aAtY, float aAtZ, float aUpX, float aUpY, float aUpZ,
		                             float aVelocityX = 0.0f, float aVelocityY = 0.0f, float aVelocityZ = 0.0f);
		// Set 3d listener position
		void set3dListenerPosition(float aPosX, float aPosY, float aPosZ);
		// Set 3d listener "at" vector
		void set3dListenerAt(float aAtX, float aAtY, float aAtZ);
		// set 3d listener "up" vector
		void set3dListenerUp(float aUpX, float aUpY, float aUpZ);
		// Set 3d listener velocity
		void set3dListenerVelocity(float aVelocityX, float aVelocityY, float aVelocityZ);

		// Set 3d audio source parameters
		void set3dSourceParameters(handle aVoiceHandle, float aPosX, float aPosY, float aPosZ, float aVelocityX = 0.0f, float aVelocityY = 0.0f,
		                           float aVelocityZ = 0.0f);
		// Set 3d audio source position
		void set3dSourcePosition(handle aVoiceHandle, float aPosX, float aPosY, float aPosZ);
		// Set 3d audio source velocity
		void set3dSourceVelocity(handle aVoiceHandle, float aVelocityX, float aVelocityY, float aVelocityZ);
		// Set 3d audio source min/max distance (distance < min means max volume)
		void set3dSourceMinMaxDistance(handle aVoiceHandle, float aMinDistance, float aMaxDistance);
		// Set 3d audio source attenuation parameters
		void set3dSourceAttenuation(handle aVoiceHandle, unsigned int aAttenuationModel, float aAttenuationRolloffFactor);
		// Set 3d audio source doppler factor to reduce or enhance doppler effect. Default = 1.0
		void set3dSourceDopplerFactor(handle aVoiceHandle, float aDopplerFactor);

		// Rest of the stuff is used internally.

		/**
		 * Apply volume scaling and clipping to audio buffer
		 *
		 * @param aSoloud       SoLoud instance (for configuration flags)
		 * @param aBuffer       Input buffer with source samples
		 * @param aDestBuffer   Output buffer for processed samples
		 * @param aSamples      Number of samples to process
		 * @param aVolume0      Starting volume level
		 * @param aVolume1      Ending volume level (for smooth ramping)
		 *
		 * Supports two clipping modes:
		 * - Hard clipping: Simple [-1,1] bounds
		 * - Roundoff clipping: Smooth saturation curve that approaches limits asymptotically
		 *
		 * Uses AVX2 intrinsics when available for 8x performance improvement,
		 * falls back to SSE intrinsics for 4x performance improvement,
		 * or scalar fallback for compatibility
		 */
		void clip_internal(AlignedFloatBuffer & aBuffer, AlignedFloatBuffer & aDestBuffer, unsigned int aSamples, float aVolume0, float aVolume1) const;

		// Returns mixed float samples in buffer. Called by the back-end, or user with null driver.
		void mix(void *aBuffer, unsigned int aSamples, detail::SAMPLE_FORMAT aFormat = detail::SAMPLE_FLOAT32);

	public:
		// Mix N samples * M channels. Called by other mix_ functions.
		void mix_internal(unsigned int aSamples, unsigned int aStride);

		// Handle rest of initialization (called from backend)
		void postinit_internal(unsigned int aSamplerate, unsigned int aBufferSize, unsigned int aFlags, unsigned int aChannels);

		// Update list of active voices
		void calcActiveVoices_internal();
		unsigned int resampleVoicePrecise_internal(AudioSourceInstance * voice, float *outputBuffer, unsigned int outputSamples, unsigned int outputStride,
		                                           double outputSampleRate, unsigned int resampler, float *scratchBuffer, unsigned int scratchSize);
		// Perform mixing for a specific bus
		void mixBus_internal(float *aBuffer, unsigned int aSamplesToRead, unsigned int aBufferSize, float *aScratch, unsigned int aBus, float aSamplerate,
		                     unsigned int aChannels, unsigned int aResampler);
		// Find a free voice, stopping the oldest if no free voice is found.
		int findFreeVoice_internal();
		// Converts handle to voice, if the handle is valid. Returns -1 if not.
		int getVoiceFromHandle_internal(handle aVoiceHandle) const;
		// Converts voice + playindex into handle
		handle getHandleFromVoice_internal(unsigned int aVoice) const;
		// Stop voice (not handle).
		void stopVoice_internal(unsigned int aVoice);
		// Set voice (not handle) pan.
		void setVoicePan_internal(unsigned int aVoice, float aPan);
		// Set voice (not handle) relative play speed.
		result setVoiceRelativePlaySpeed_internal(unsigned int aVoice, float aSpeed);
		// Set voice (not handle) volume.
		void setVoiceVolume_internal(unsigned int aVoice, float aVolume);
		// Set voice (not handle) pause state.
		void setVoicePause_internal(unsigned int aVoice, int aPause);
		// Update overall volume from set and 3d volumes
		void updateVoiceVolume_internal(unsigned int aVoice);
		// Update overall relative play speed from set and 3d speeds
		void updateVoiceRelativePlaySpeed_internal(unsigned int aVoice);
		// Perform 3d audio calculation for array of voices
		void update3dVoices_internal(unsigned int *aVoiceList, unsigned int aVoiceCount);
		// Remove all non-active voices from group
		void trimVoiceGroup_internal(handle aVoiceGroupHandle);
		// Get pointer to the zero-terminated array of voice handles in a voice group
		handle *voiceGroupHandleToArray_internal(handle aVoiceGroupHandle) const;

		// Lock audio thread mutex.
		void lockAudioMutex_internal();
		// Unlock audio thread mutex.
		void unlockAudioMutex_internal();

		// Max. number of active voices. Busses and tickable inaudibles also count against this.
		unsigned int mMaxActiveVoices;
		// Highest voice in use so far
		unsigned int mHighestVoice;
		// Scratch buffer, used for resampling.
		AlignedFloatBuffer mScratch;
		// Current size of the scratch, in samples.
		unsigned int mScratchSize;
		// Output scratch buffer, used in mix_().
		AlignedFloatBuffer mOutputScratch;
		// Audio voices.
		AudioSourceInstance *mVoice[VOICE_COUNT];
		// Resampler for the main bus
		unsigned int mResampler;
		// Output sample rate (not float)
		unsigned int mSamplerate;
		// Output channel count
		unsigned int mChannels;
		// Current backend ID
		unsigned int mBackendID;
		// Current backend string
		const char *mBackendString;
		// Maximum size of output buffer; used to calculate needed scratch.
		unsigned int mBufferSize;
		// Flags; see Soloud::FLAGS
		unsigned int mFlags;
		// Global volume. Applied before clipping.
		float mGlobalVolume;
		// Post-clip scaler. Applied after clipping.
		float mPostClipScaler;
		// Current play index. Used to create audio handles.
		unsigned int mPlayIndex;
		// Current sound source index. Used to create sound source IDs.
		unsigned int mAudioSourceID;
		// Fader for the global volume.
		Fader mGlobalVolumeFader;
		// Global stream time, for the global volume fader.
		time mStreamTime;
		// Last time seen by the playClocked call
		time mLastClockedTime;
		// Global filter
		Filter *mFilter[FILTERS_PER_STREAM];
		// Global filter instance
		FilterInstance *mFilterInstance[FILTERS_PER_STREAM];

		// Approximate volume for channels.
		float mVisualizationChannelVolume[MAX_CHANNELS];
		// Mono-mixed wave data for visualization and for visualization FFT input
		float mVisualizationWaveData[256];
		// FFT output data
		float mFFTData[256];
		// Snapshot of wave data for visualization
		float mWaveData[256];

		// 3d listener position
		float m3dPosition[3];
		// 3d listener look-at
		float m3dAt[3];
		// 3d listener up
		float m3dUp[3];
		// 3d listener velocity
		float m3dVelocity[3];
		// 3d speed of sound (for doppler)
		float m3dSoundSpeed;

		// 3d position of speakers
		float m3dSpeakerPosition[3 * MAX_CHANNELS];

		// Data related to 3d processing, separate from AudioSource so we can do 3d calculations without audio mutex.
		AudioSourceInstance3dData m3dData[VOICE_COUNT];

		// For each voice group, first int is number of ints alocated.
		unsigned int **mVoiceGroup;
		unsigned int mVoiceGroupCount;

		// List of currently active voices
		unsigned int mActiveVoice[VOICE_COUNT];
		// Number of currently active voices
		unsigned int mActiveVoiceCount;
		// Active voices list needs to be recalculated
		bool mActiveVoiceDirty;
	};
}; // namespace SoLoud

#endif
