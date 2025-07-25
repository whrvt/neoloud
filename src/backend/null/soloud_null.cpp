/*
SoLoud audio engine
Copyright (c) 2013-2015 Jari Komppa

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

namespace SoLoud
{
    static void nullCleanup(Soloud * /*aSoloud*/)
    {
    }

    result null_init(Soloud *aSoloud, unsigned int aFlags, unsigned int aSamplerate, unsigned int aBuffer, unsigned int aChannels)
    {
		if (aChannels == 0 || aChannels == 3 || aChannels == 5 || aChannels == 7 || aChannels > MAX_CHANNELS || aBuffer < SAMPLE_GRANULARITY)
			return INVALID_PARAMETER;
        aSoloud->mBackendData = nullptr;
        aSoloud->mBackendCleanupFunc = nullCleanup;

        aSoloud->postinit_internal(aSamplerate, aBuffer, aFlags, aChannels);
        aSoloud->mBackendString = "null driver";
        return SO_NO_ERROR;
    }
};
