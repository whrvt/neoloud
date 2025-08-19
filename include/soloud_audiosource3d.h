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

#ifndef SOLOUD_AUDIOSOURCE3D_H
#define SOLOUD_AUDIOSOURCE3D_H

#include "soloud_config.h"

namespace SoLoud
{
class AudioSource;
class AudioCollider;
class AudioAttenuator;

class AudioSourceInstance3dData
{
public:
	// ctor
	AudioSourceInstance3dData();
	// Set settings from audiosource
	void init(AudioSource &aSource);
	// 3d position
	float m3dPosition[3];
	// 3d velocity
	float m3dVelocity[3];
	// 3d cone direction
	/*
	float m3dConeDirection[3];
	// 3d cone inner angle
	float m3dConeInnerAngle;
	// 3d cone outer angle
	float m3dConeOuterAngle;
	// 3d cone outer volume multiplier
	float m3dConeOuterVolume;
	*/
	// 3d min distance
	float m3dMinDistance;
	// 3d max distance
	float m3dMaxDistance;
	// 3d attenuation rolloff factor
	float m3dAttenuationRolloff;
	// 3d attenuation model
	unsigned int m3dAttenuationModel;
	// 3d doppler factor
	float m3dDopplerFactor;
	// Pointer to a custom audio collider object
	AudioCollider *mCollider;
	// Pointer to a custom audio attenuator object
	AudioAttenuator *mAttenuator;
	// User data related to audio collider
	int mColliderData;

	// Doppler sample rate multiplier
	float mDopplerValue;
	// Overall 3d volume
	float m3dVolume;
	// Channel volume
	float mChannelVolume[MAX_CHANNELS];
	// Copy of flags
	unsigned int mFlags;
	// Latest handle for this voice
	handle mHandle;
};
} // namespace SoLoud

#endif
