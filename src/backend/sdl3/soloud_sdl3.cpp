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
#include <cstdlib>

#include "soloud_internal.h"

#if !defined(WITH_SDL3)

namespace SoLoud
{
    result sdl3_init(SoLoud::Soloud *aSoloud, unsigned int aFlags, unsigned int aSamplerate, unsigned int aBuffer, unsigned int aChannels)
    {
        return NOT_IMPLEMENTED;
    }
}

#else

#include <SDL3/SDL_init.h>
#include <SDL3/SDL_audio.h>
#include <math.h>

namespace SoLoud
{
    static bool gWeInitSDL = false;
    static SDL_AudioDeviceID gAudioDeviceID = 0;
    static float* gFloatBuffer = nullptr;
    static short* gShortBuffer = nullptr;
    static int gBufferSize = 0;
    static SDL_AudioFormat gDeviceFormat = SDL_AUDIO_UNKNOWN;

    static void soloud_sdl3_postmix(void *userdata, const SDL_AudioSpec *spec, float *buffer, int buflen)
    {
        SoLoud::Soloud *soloud = (SoLoud::Soloud *)userdata;
        if (!soloud) return;

        int sampleFrames = buflen / (spec->channels * sizeof(float));
        soloud->mix(buffer, sampleFrames);
    }

    static void soloud_sdl3_deinit(SoLoud::Soloud *aSoloud)
    {
        if (gFloatBuffer) {
            SDL_free(gFloatBuffer);
            gFloatBuffer = nullptr;
        }

        if (gShortBuffer) {
            SDL_free(gShortBuffer);
            gShortBuffer = nullptr;
        }

        gBufferSize = 0;

        if (gAudioDeviceID) {
            SDL_SetAudioPostmixCallback(gAudioDeviceID, NULL, NULL);
            SDL_CloseAudioDevice(gAudioDeviceID);
            gAudioDeviceID = 0;
        }

        if (gWeInitSDL && SDL_WasInit(SDL_INIT_AUDIO)) {
            SDL_QuitSubSystem(SDL_INIT_AUDIO);
            gWeInitSDL = false;
        }
    }

    result sdl3_init(SoLoud::Soloud *aSoloud, unsigned int aFlags, unsigned int aSamplerate, unsigned int aBuffer, unsigned int aChannels)
    {
        if (!aSoloud) return INVALID_PARAMETER;

        if (gAudioDeviceID) {
            soloud_sdl3_deinit(aSoloud);
        }

        if (!SDL_WasInit(SDL_INIT_AUDIO)) {
            if (!SDL_InitSubSystem(SDL_INIT_AUDIO)) {
                return UNKNOWN_ERROR;
            }
            gWeInitSDL = true;
        }

        SDL_AudioSpec preferredSpec;
        bool hasPreferredFormat = SDL_GetAudioDeviceFormat(SDL_AUDIO_DEVICE_DEFAULT_PLAYBACK,
                                                          &preferredSpec,
                                                          NULL);

        aSamplerate = hasPreferredFormat ? preferredSpec.freq : aSamplerate;
        aChannels = hasPreferredFormat ? preferredSpec.channels : aChannels;

        SDL_AudioSpec spec;
        spec.freq = aSamplerate;
        spec.channels = aChannels;
        spec.format = SDL_AUDIO_F32;

        gAudioDeviceID = SDL_OpenAudioDevice(SDL_AUDIO_DEVICE_DEFAULT_PLAYBACK, &spec);
        if (gAudioDeviceID == 0) {
            spec.format = SDL_AUDIO_S16;
            gAudioDeviceID = SDL_OpenAudioDevice(SDL_AUDIO_DEVICE_DEFAULT_PLAYBACK, &spec);

            if (gAudioDeviceID == 0) {
                gAudioDeviceID = SDL_OpenAudioDevice(SDL_AUDIO_DEVICE_DEFAULT_PLAYBACK, NULL);
                if (gAudioDeviceID == 0) {
                    return UNKNOWN_ERROR;
                }
            }
        }

        SDL_AudioSpec actualSpec;
        int deviceSampleFrames;
        if (!SDL_GetAudioDeviceFormat(gAudioDeviceID, &actualSpec, &deviceSampleFrames)) {
            SDL_CloseAudioDevice(gAudioDeviceID);
            gAudioDeviceID = 0;
            return UNKNOWN_ERROR;
        }

        gDeviceFormat = actualSpec.format;

        unsigned int soloudBufferSize = aBuffer;
        if (deviceSampleFrames > 0 && deviceSampleFrames < (int)soloudBufferSize) {
            // 2048 as in default/not set, otherwise respect the set choice if the queried buffer size is way too small
            soloudBufferSize = (((deviceSampleFrames >= 16 ? deviceSampleFrames : (soloudBufferSize == 2048 ? 16 : soloudBufferSize)) + 1) / 2) * 2;
        }

        if (!SDL_SetAudioPostmixCallback(gAudioDeviceID, soloud_sdl3_postmix, aSoloud)) {
            SDL_CloseAudioDevice(gAudioDeviceID);
            gAudioDeviceID = 0;
            return UNKNOWN_ERROR;
        }

        aSoloud->postinit_internal(actualSpec.freq, soloudBufferSize, aFlags, actualSpec.channels);
        aSoloud->mBackendCleanupFunc = soloud_sdl3_deinit;
        aSoloud->mBackendString = "SDL3 (static)";

        return 0;
    }
};
#endif
