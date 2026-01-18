/*
SoLoud audio engine - low-level mixing (internal)
Copyright (c) 2013-2020 Jari Komppa
Copyright (c) 2025-2026 William Horvath

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

#ifndef SOLOUD_MIXING_INTERNAL_H
#define SOLOUD_MIXING_INTERNAL_H

#include "soloud_cpu.h"
#if defined(SOLOUD_IS_X86) && !defined(DISABLE_SIMD)
#include "soloud_ll_mixing.h"
#endif

namespace SoLoud::mixing
{
// Resampling algorithm constants
namespace ResamplingConstants
{
// Lookahead samples required for each resampling algorithm
inline constexpr unsigned int POINT_LOOKAHEAD = 1;      // Point sampling: just current sample
inline constexpr unsigned int LINEAR_LOOKAHEAD = 2;     // Linear: current + next sample
inline constexpr unsigned int CATMULLROM_LOOKAHEAD = 4; // Catmull-Rom: 4-point cubic interpolation

// Safety margin for high sample rate ratios to prevent buffer underruns
inline constexpr unsigned int LOOKAHEAD_SAFETY_MARGIN = 8;

// Catmull-Rom cubic interpolation coefficients
// Formula: 0.5 * ((2*p1) + (-p0+p2)*t + (2*p0-5*p1+4*p2-p3)*t^2 + (-p0+3*p1-3*p2+p3)*t^3)
inline constexpr float CATMULLROM_SCALE = 0.5f;
inline constexpr float CATMULLROM_LINEAR_COEFF = 2.0f;
inline constexpr float CATMULLROM_QUAD_COEFFS[4]{2.0f, -5.0f, 4.0f, -1.0f};  // p0, p1, p2, p3 coefficients
inline constexpr float CATMULLROM_CUBIC_COEFFS[4]{-1.0f, 3.0f, -3.0f, 1.0f}; // p0, p1, p2, p3 coefficients
} // namespace ResamplingConstants

// Roundoff clipping algorithm constants
// These values provide smooth saturation instead of hard clipping
namespace ClippingConstants
{
// Input threshold bounds for roundoff clipping
inline constexpr float ROUNDOFF_NEG_THRESHOLD = -1.65f;
inline constexpr float ROUNDOFF_POS_THRESHOLD = 1.65f;

// Output saturation limits for roundoff clipping
inline constexpr float ROUNDOFF_NEG_WALL = -0.9862875f;
inline constexpr float ROUNDOFF_POS_WALL = 0.9862875f;

// Roundoff curve coefficients: output = LINEAR_SCALE * input + CUBIC_SCALE * input^3
// This creates a smooth S-curve that approaches the walls asymptotically
inline constexpr float ROUNDOFF_LINEAR_SCALE = 0.87f;
inline constexpr float ROUNDOFF_CUBIC_SCALE = -0.1f;

// Hard clipping bounds
inline constexpr float HARD_CLIP_MIN = -1.0f;
inline constexpr float HARD_CLIP_MAX = 1.0f;
} // namespace ClippingConstants

// Channel mixing coefficients for downmixing operations
namespace ChannelMixingConstants
{
// Coefficients for mixing multi-channel content to fewer channels
// These values are chosen to preserve perceived loudness while avoiding clipping

inline constexpr float MIX_8_TO_2_SCALE = 0.2f; // 8->2: Mix 5 channels per output
inline constexpr float MIX_6_TO_2_SCALE = 0.3f; // 6->2: Mix 4 channels per output
inline constexpr float MIX_4_TO_2_SCALE = 0.5f; // 4->2: Mix 2 channels per output

inline constexpr float CENTER_SUB_MIX_SCALE = 0.7f; // Center + sub mixing level
inline constexpr float SURROUND_MIX_SCALE = 0.5f;   // Surround channel mixing level
inline constexpr float QUAD_MIX_SCALE = 0.25f;      // 4-channel average mixing
} // namespace ChannelMixingConstants

// Declare AVX/SSE-optimized mixer implementations here, as well.
#if defined(SOLOUD_IS_X86) && !defined(DISABLE_SIMD)
class MixerAVX final : public Mixer
{
	friend Mixer;
	// Created by Mixer::createMixer()
	MixerAVX() = default;

public:
	~MixerAVX() override = default;

	MixerAVX(const MixerAVX &) = delete;
	MixerAVX &operator=(const MixerAVX &) = delete;
	MixerAVX(MixerAVX &&) = delete;
	MixerAVX &operator=(MixerAVX &&) = delete;

	void resample_channels(float **srcChannels, float *outputBuffer, unsigned int outputStride, unsigned int numChannels, unsigned int outputSamples,
	                       double srcPosition, double stepSize, unsigned int availableInput, unsigned int lookahead) override;
	void clip_samples(const float *const aBuffer, float *aDestBuffer, unsigned int aSamples, unsigned int aChannels, float aVolume0, float aVolume1, float aScaler,
	                  bool aRoundoff) override;
	void interlace_samples(void *outputBuffer, const float *const rawBuffer, unsigned int aSamples, unsigned int stride, unsigned int aChannels,
	                       SAMPLE_FORMAT aFormat = SAMPLE_FLOAT32) override;
	void panAndExpand(AudioSourceInstance *aVoice, float *aBuffer, unsigned int aSamplesToRead, unsigned int aBufferSize, float *aScratch,
	                  unsigned int aChannels) override;
};

class MixerSSE final : public Mixer
{
	friend Mixer;
	// Created by Mixer::createMixer()
	MixerSSE() = default;

public:
	~MixerSSE() override = default;

	MixerSSE(const MixerSSE &) = delete;
	MixerSSE &operator=(const MixerSSE &) = delete;
	MixerSSE(MixerSSE &&) = delete;
	MixerSSE &operator=(MixerSSE &&) = delete;

	void resample_channels(float **srcChannels, float *outputBuffer, unsigned int outputStride, unsigned int numChannels, unsigned int outputSamples,
	                       double srcPosition, double stepSize, unsigned int availableInput, unsigned int lookahead) override;
	void clip_samples(const float *const aBuffer, float *aDestBuffer, unsigned int aSamples, unsigned int aChannels, float aVolume0, float aVolume1, float aScaler,
	                  bool aRoundoff) override;
	void interlace_samples(void *outputBuffer, const float *const rawBuffer, unsigned int aSamples, unsigned int stride, unsigned int aChannels,
	                       SAMPLE_FORMAT aFormat = SAMPLE_FLOAT32) override;
	void panAndExpand(AudioSourceInstance *aVoice, float *aBuffer, unsigned int aSamplesToRead, unsigned int aBufferSize, float *aScratch,
	                  unsigned int aChannels) override;
};

#endif

} // namespace SoLoud::mixing

#endif
