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

#include "soloud.h"
#include "soloud_audiosource.h"
#include "soloud_cpu.h"
#include "soloud_ll_mixing.h"
#include "soloud_mixing_internal.h"

#include <climits> // _MAX
#include <cstring>

// NOLINTBEGIN(cppcoreguidelines-avoid-c-arrays, modernize-avoid-c-arrays, cppcoreguidelines-pro-type-member-init, hicpp-member-init)

namespace SoLoud::mixing
{
// Factory
Mixer *Mixer::createMixer()
{
#if !defined(SOLOUD_DISABLE_SIMD)
	const unsigned int CPUType = detectCPUextensions();
#if defined(SOLOUD_SUPPORT_AVX2)
	if (CPUType & CPUFEATURE_AVX2)
	{
		return new MixerAVX();
	}
	else
#endif
#if defined(SOLOUD_SUPPORT_SSE2)
	    if (CPUType & CPUFEATURE_SSE2)
	{
		return new MixerSSE();
	}
	else
#endif
#endif
	{
		// Base implementation.
		return new Mixer();
	}
}

namespace
{
constexpr size_t SCALAR_OPTIMAL_CHUNK_SAMPLES = 4; // Process 4 samples at a time even for scalar
}

void Mixer::interlace_samples(void *outputBuffer, const float *const rawBuffer, unsigned int aSamples, unsigned int stride, unsigned int aChannels,
                              SAMPLE_FORMAT aFormat)
{
	unsigned int j = 0, c = 0;
	// scalar implementation
	int i = 0;
	switch (aFormat)
	{
	case SAMPLE_FLOAT32: {
		float *buffer = static_cast<float *>(outputBuffer);
		for (j = 0; j < aChannels; j++)
		{
			c = j * stride;
			for (i = j; i < aSamples * aChannels; i += aChannels)
			{
				buffer[i] = rawBuffer[c];
				c++;
			}
		}
	}
	break;

	case SAMPLE_UNSIGNED8: {
		unsigned char *buffer = static_cast<unsigned char *>(outputBuffer);
		for (j = 0; j < aChannels; j++)
		{
			c = j * stride;
			for (i = j; i < aSamples * aChannels; i += aChannels)
			{
				int sample = (int)(rawBuffer[c] * 127.0f + 128.0f);
				if (sample < 0)
					sample = 0;
				if (sample > 255)
					sample = 255;
				buffer[i] = (unsigned char)sample;
				c++;
			}
		}
	}
	break;

	case SAMPLE_SIGNED16: {
		short *buffer = static_cast<short *>(outputBuffer);
		for (j = 0; j < aChannels; j++)
		{
			c = j * stride;
			for (i = j; i < aSamples * aChannels; i += aChannels)
			{
				buffer[i] = (short)(rawBuffer[c] * 0x7fff);
				c++;
			}
		}
	}
	break;

	case SAMPLE_SIGNED24: {
		unsigned char *buffer = static_cast<unsigned char *>(outputBuffer);
		for (j = 0; j < aChannels; j++)
		{
			c = j * stride;
			for (i = j; i < aSamples * aChannels; i += aChannels)
			{
				int sample = (int)(rawBuffer[c] * (float)(INT_MAX >> 8));
				if (sample < (INT_MIN >> 8))
					sample = (INT_MIN >> 8);
				if (sample > (INT_MAX >> 8))
					sample = (INT_MAX >> 8);

				unsigned int destIdx = i * 3;
				buffer[destIdx] = (unsigned char)(sample & 0xff);
				buffer[destIdx + 1] = (unsigned char)((sample >> 8) & 0xff);
				buffer[destIdx + 2] = (unsigned char)((sample >> 16) & 0xff);
				c++;
			}
		}
	}
	break;

	case SAMPLE_SIGNED32: {
		int *buffer = static_cast<int *>(outputBuffer);
		for (j = 0; j < aChannels; j++)
		{
			c = j * stride;
			for (i = j; i < aSamples * aChannels; i += aChannels)
			{
				double sample = (double)rawBuffer[c] * (double)INT_MAX;
				if (sample < (double)INT_MIN)
					sample = (double)INT_MIN;
				if (sample > (double)INT_MAX)
					sample = (double)INT_MAX;
				buffer[i] = (int)sample;
				c++;
			}
		}
	}
	break;
	}
}

// Scalar implementation
void Mixer::clip_samples(const float *const aBuffer, float *aDestBuffer, unsigned int aSamples, unsigned int aChannels, float aVolume0, float aVolume1,
                         float aScaler, bool aRoundoff)
{
	using namespace ClippingConstants;

	float volumeDelta = (aVolume1 - aVolume0) / aSamples;
	float volume = aVolume0;
	unsigned int sampleQuads = (aSamples + 3) / 4; // Process in groups of 4 for consistency

	if (aRoundoff)
	{
		unsigned int srcIdx = 0, dstIdx = 0;

		for (unsigned int channel = 0; channel < aChannels; channel++)
		{
			volume = aVolume0;

			for (unsigned int quad = 0; quad < sampleQuads; quad++)
			{
				// Process 4 samples at a time
				for (unsigned int sample = 0; sample < SCALAR_OPTIMAL_CHUNK_SAMPLES; sample++)
				{
					float input = aBuffer[srcIdx++] * volume;
					volume += volumeDelta;

					// Apply roundoff clipping curve
					float output;
					if (input <= ROUNDOFF_NEG_THRESHOLD)
					{
						output = ROUNDOFF_NEG_WALL;
					}
					else if (input >= ROUNDOFF_POS_THRESHOLD)
					{
						output = ROUNDOFF_POS_WALL;
					}
					else
					{
						// Smooth curve: linear + cubic term
						output = ROUNDOFF_LINEAR_SCALE * input + ROUNDOFF_CUBIC_SCALE * input * input * input;
					}

					aDestBuffer[dstIdx++] = output * aScaler;
				}
			}
		}
	}
	else
	{
		unsigned int srcIdx = 0, dstIdx = 0;

		for (unsigned int channel = 0; channel < aChannels; channel++)
		{
			volume = aVolume0;

			for (unsigned int quad = 0; quad < sampleQuads; quad++)
			{
				// Process 4 samples at a time
				for (unsigned int sample = 0; sample < SCALAR_OPTIMAL_CHUNK_SAMPLES; sample++)
				{
					float input = aBuffer[srcIdx++] * volume;
					volume += volumeDelta;

					// Hard clipping to [-1, 1] range
					float output = (input <= HARD_CLIP_MIN) ? HARD_CLIP_MIN : (input >= HARD_CLIP_MAX) ? HARD_CLIP_MAX : input;

					aDestBuffer[dstIdx++] = output * aScaler;
				}
			}
		}
	}
}

void Mixer::panAndExpand(AudioSourceInstance *aVoice, float *aBuffer, unsigned int aSamplesToRead, unsigned int aBufferSize, float *aScratch, unsigned int aChannels)
{
	using namespace ChannelMixingConstants;
	using namespace ClippingConstants;

	// Calculate volume ramping for smooth transitions
	float currentPan[MAX_CHANNELS]; // Current speaker volumes
	float targetPan[MAX_CHANNELS];  // Target speaker volumes
	float panDelta[MAX_CHANNELS];   // Volume change per sample

	for (unsigned int ch = 0; ch < aChannels; ch++)
	{
		currentPan[ch] = aVoice->mCurrentChannelVolume[ch];
		targetPan[ch] = aVoice->mChannelVolume[ch] * aVoice->mOverallVolume;
		panDelta[ch] = (targetPan[ch] - currentPan[ch]) / aSamplesToRead;
	}

	// Channel mapping and mixing logic
	switch (aChannels)
	{
	case 1: // Mono output: sum all input channels
	{
		for (unsigned int srcCh = 0; srcCh < aVoice->mChannels; srcCh++)
		{
			unsigned int srcOffset = srcCh * aBufferSize;
			currentPan[0] = aVoice->mCurrentChannelVolume[0];

			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				currentPan[0] += panDelta[0];
				aBuffer[sample] += aScratch[srcOffset + sample] * currentPan[0];
			}
		}
	}
	break;

	case 2: // Stereo output
		switch (aVoice->mChannels)
		{
		case 8: // 8.0 -> 2.0: mix left/right groups with center and sub contribution
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				currentPan[0] += panDelta[0];
				currentPan[1] += panDelta[1];

				// Load all 8 channels
				float ch[8];
				for (int i = 0; i < 8; i++)
				{
					ch[i] = aScratch[aBufferSize * i + sample];
				}

				// Mix: FL + C + LFE + SL + BL for left, FR + C + LFE + SR + BR for right
				aBuffer[sample] += MIX_8_TO_2_SCALE * (ch[0] + ch[2] + ch[3] + ch[4] + ch[6]) * currentPan[0];
				aBuffer[aBufferSize + sample] += MIX_8_TO_2_SCALE * (ch[1] + ch[2] + ch[3] + ch[5] + ch[7]) * currentPan[1];
			}
			break;

		case 6: // 5.1 -> 2.0: mix left/right groups with center and sub contribution
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				currentPan[0] += panDelta[0];
				currentPan[1] += panDelta[1];

				// Load all 6 channels: FL, FR, C, LFE, SL, SR
				float ch[6];
				for (int i = 0; i < 6; i++)
				{
					ch[i] = aScratch[aBufferSize * i + sample];
				}

				// Mix: FL + C + LFE + SL for left, FR + C + LFE + SR for right
				aBuffer[sample] += MIX_6_TO_2_SCALE * (ch[0] + ch[2] + ch[3] + ch[4]) * currentPan[0];
				aBuffer[aBufferSize + sample] += MIX_6_TO_2_SCALE * (ch[1] + ch[2] + ch[3] + ch[5]) * currentPan[1];
			}
			break;

		case 4: // 4.0 -> 2.0: sum front/back pairs
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				currentPan[0] += panDelta[0];
				currentPan[1] += panDelta[1];

				float frontLeft = aScratch[sample];
				float frontRight = aScratch[aBufferSize + sample];
				float backLeft = aScratch[aBufferSize * 2 + sample];
				float backRight = aScratch[aBufferSize * 3 + sample];

				aBuffer[sample] += MIX_4_TO_2_SCALE * (frontLeft + backLeft) * currentPan[0];
				aBuffer[aBufferSize + sample] += MIX_4_TO_2_SCALE * (frontRight + backRight) * currentPan[1];
			}
			break;

		case 2: // 2.0 -> 2.0: direct mapping with volume ramping
		        // Scalar implementation
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				currentPan[0] += panDelta[0];
				currentPan[1] += panDelta[1];

				float leftSrc = aScratch[sample];
				float rightSrc = aScratch[aBufferSize + sample];

				aBuffer[sample] += leftSrc * currentPan[0];
				aBuffer[aBufferSize + sample] += rightSrc * currentPan[1];
			}
			break;

		case 1: // 1.0 -> 2.0: distribute mono to stereo with panning
		        // Scalar implementation
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				currentPan[0] += panDelta[0];
				currentPan[1] += panDelta[1];

				float monoSrc = aScratch[sample];

				aBuffer[sample] += monoSrc * currentPan[0];
				aBuffer[aBufferSize + sample] += monoSrc * currentPan[1];
			}
			break;
		}
		break;

	case 4: // Quadraphonic output
		switch (aVoice->mChannels)
		{
		case 8: // 8.0 -> 4.0: group channels and add center/sub
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				for (int ch = 0; ch < 4; ch++)
				{
					currentPan[ch] += panDelta[ch];
				}

				float ch[8];
				for (int i = 0; i < 8; i++)
				{
					ch[i] = aScratch[aBufferSize * i + sample];
				}

				float centerSub = (ch[2] + ch[3]) * CENTER_SUB_MIX_SCALE;

				aBuffer[sample] += ch[0] * currentPan[0] + centerSub;
				aBuffer[aBufferSize + sample] += ch[1] * currentPan[1] + centerSub;
				aBuffer[aBufferSize * 2 + sample] += SURROUND_MIX_SCALE * (ch[4] + ch[6]) * currentPan[2];
				aBuffer[aBufferSize * 3 + sample] += SURROUND_MIX_SCALE * (ch[5] + ch[7]) * currentPan[3];
			}
			break;

		case 6: // 5.1 -> 4.0: map channels and add center/sub
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				for (int ch = 0; ch < 4; ch++)
				{
					currentPan[ch] += panDelta[ch];
				}

				float ch[6];
				for (int i = 0; i < 6; i++)
				{
					ch[i] = aScratch[aBufferSize * i + sample];
				}

				float centerSub = (ch[2] + ch[3]) * CENTER_SUB_MIX_SCALE;

				aBuffer[sample] += ch[0] * currentPan[0] + centerSub;
				aBuffer[aBufferSize + sample] += ch[1] * currentPan[1] + centerSub;
				aBuffer[aBufferSize * 2 + sample] += ch[4] * currentPan[2];
				aBuffer[aBufferSize * 3 + sample] += ch[5] * currentPan[3];
			}
			break;

		case 4: // 4.0 -> 4.0: direct mapping
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				for (int ch = 0; ch < 4; ch++)
				{
					currentPan[ch] += panDelta[ch];
					float src = aScratch[aBufferSize * ch + sample];
					aBuffer[aBufferSize * ch + sample] += src * currentPan[ch];
				}
			}
			break;

		case 2: // 2.0 -> 4.0: duplicate to front and back
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				for (int ch = 0; ch < 4; ch++)
				{
					currentPan[ch] += panDelta[ch];
				}

				float left = aScratch[sample];
				float right = aScratch[aBufferSize + sample];

				aBuffer[sample] += left * currentPan[0];
				aBuffer[aBufferSize + sample] += right * currentPan[1];
				aBuffer[aBufferSize * 2 + sample] += left * currentPan[2];
				aBuffer[aBufferSize * 3 + sample] += right * currentPan[3];
			}
			break;

		case 1: // 1.0 -> 4.0: distribute mono to all channels
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				for (int ch = 0; ch < 4; ch++)
				{
					currentPan[ch] += panDelta[ch];
				}

				float monoSrc = aScratch[sample];

				for (int ch = 0; ch < 4; ch++)
				{
					aBuffer[aBufferSize * ch + sample] += monoSrc * currentPan[ch];
				}
			}
			break;
		}
		break;

	case 6: // 5.1 Surround output
		switch (aVoice->mChannels)
		{
		case 8: // 8.0 -> 5.1: combine side and back channels
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				for (int ch = 0; ch < 6; ch++)
				{
					currentPan[ch] += panDelta[ch];
				}

				float ch[8];
				for (int i = 0; i < 8; i++)
				{
					ch[i] = aScratch[aBufferSize * i + sample];
				}

				// FL, FR, C, LFE, SL+BL, SR+BR
				aBuffer[sample] += ch[0] * currentPan[0];
				aBuffer[aBufferSize + sample] += ch[1] * currentPan[1];
				aBuffer[aBufferSize * 2 + sample] += ch[2] * currentPan[2];
				aBuffer[aBufferSize * 3 + sample] += ch[3] * currentPan[3];
				aBuffer[aBufferSize * 4 + sample] += SURROUND_MIX_SCALE * (ch[4] + ch[6]) * currentPan[4];
				aBuffer[aBufferSize * 5 + sample] += SURROUND_MIX_SCALE * (ch[5] + ch[7]) * currentPan[5];
			}
			break;

		case 6: // 5.1 -> 5.1: direct mapping
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				for (int ch = 0; ch < 6; ch++)
				{
					currentPan[ch] += panDelta[ch];
					float src = aScratch[aBufferSize * ch + sample];
					aBuffer[aBufferSize * ch + sample] += src * currentPan[ch];
				}
			}
			break;

		case 4: // 4.0 -> 5.1: create center/sub from front channels
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				for (int ch = 0; ch < 6; ch++)
				{
					currentPan[ch] += panDelta[ch];
				}

				float ch[4];
				for (int i = 0; i < 4; i++)
				{
					ch[i] = aScratch[aBufferSize * i + sample];
				}

				// FL, FR, C (mix), LFE (mix), SL, SR
				aBuffer[sample] += ch[0] * currentPan[0];
				aBuffer[aBufferSize + sample] += ch[1] * currentPan[1];
				aBuffer[aBufferSize * 2 + sample] += SURROUND_MIX_SCALE * (ch[0] + ch[1]) * currentPan[2];
				aBuffer[aBufferSize * 3 + sample] += QUAD_MIX_SCALE * (ch[0] + ch[1] + ch[2] + ch[3]) * currentPan[3];
				aBuffer[aBufferSize * 4 + sample] += ch[2] * currentPan[4];
				aBuffer[aBufferSize * 5 + sample] += ch[3] * currentPan[5];
			}
			break;

		case 2: // 2.0 -> 5.1: create surround content
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				for (int ch = 0; ch < 6; ch++)
				{
					currentPan[ch] += panDelta[ch];
				}

				float left = aScratch[sample];
				float right = aScratch[aBufferSize + sample];

				// FL, FR, C (mix), LFE (mix), SL, SR
				aBuffer[sample] += left * currentPan[0];
				aBuffer[aBufferSize + sample] += right * currentPan[1];
				aBuffer[aBufferSize * 2 + sample] += SURROUND_MIX_SCALE * (left + right) * currentPan[2];
				aBuffer[aBufferSize * 3 + sample] += SURROUND_MIX_SCALE * (left + right) * currentPan[3];
				aBuffer[aBufferSize * 4 + sample] += left * currentPan[4];
				aBuffer[aBufferSize * 5 + sample] += right * currentPan[5];
			}
			break;

		case 1: // 1.0 -> 5.1: distribute mono to all channels
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				for (int ch = 0; ch < 6; ch++)
				{
					currentPan[ch] += panDelta[ch];
				}

				float monoSrc = aScratch[sample];

				for (int ch = 0; ch < 6; ch++)
				{
					aBuffer[aBufferSize * ch + sample] += monoSrc * currentPan[ch];
				}
			}
			break;
		}
		break;

	case 8: // 7.1 Surround output
		switch (aVoice->mChannels)
		{
		case 8: // 8.0 -> 7.1: direct mapping
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				for (int ch = 0; ch < 8; ch++)
				{
					currentPan[ch] += panDelta[ch];
					float src = aScratch[aBufferSize * ch + sample];
					aBuffer[aBufferSize * ch + sample] += src * currentPan[ch];
				}
			}
			break;

		case 6: // 5.1 -> 7.1: create side channels from surrounds
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				for (int ch = 0; ch < 8; ch++)
				{
					currentPan[ch] += panDelta[ch];
				}

				float ch[6];
				for (int i = 0; i < 6; i++)
				{
					ch[i] = aScratch[aBufferSize * i + sample];
				}

				// FL, FR, C, LFE, SL+blend, SR+blend, SL, SR
				aBuffer[sample] += ch[0] * currentPan[0];
				aBuffer[aBufferSize + sample] += ch[1] * currentPan[1];
				aBuffer[aBufferSize * 2 + sample] += ch[2] * currentPan[2];
				aBuffer[aBufferSize * 3 + sample] += ch[3] * currentPan[3];
				aBuffer[aBufferSize * 4 + sample] += SURROUND_MIX_SCALE * (ch[4] + ch[0]) * currentPan[4];
				aBuffer[aBufferSize * 5 + sample] += SURROUND_MIX_SCALE * (ch[5] + ch[1]) * currentPan[5];
				aBuffer[aBufferSize * 6 + sample] += ch[4] * currentPan[6];
				aBuffer[aBufferSize * 7 + sample] += ch[5] * currentPan[7];
			}
			break;

		case 4: // 4.0 -> 7.1: create center/sub and side channels
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				for (int ch = 0; ch < 8; ch++)
				{
					currentPan[ch] += panDelta[ch];
				}

				float ch[4];
				for (int i = 0; i < 4; i++)
				{
					ch[i] = aScratch[aBufferSize * i + sample];
				}

				// FL, FR, C (mix), LFE (mix), SL+FL, SR+FR, SL, SR
				aBuffer[sample] += ch[0] * currentPan[0];
				aBuffer[aBufferSize + sample] += ch[1] * currentPan[1];
				aBuffer[aBufferSize * 2 + sample] += SURROUND_MIX_SCALE * (ch[0] + ch[1]) * currentPan[2];
				aBuffer[aBufferSize * 3 + sample] += QUAD_MIX_SCALE * (ch[0] + ch[1] + ch[2] + ch[3]) * currentPan[3];
				aBuffer[aBufferSize * 4 + sample] += SURROUND_MIX_SCALE * (ch[0] + ch[2]) * currentPan[4];
				aBuffer[aBufferSize * 5 + sample] += SURROUND_MIX_SCALE * (ch[1] + ch[3]) * currentPan[5];
				aBuffer[aBufferSize * 6 + sample] += ch[2] * currentPan[4];
				aBuffer[aBufferSize * 7 + sample] += ch[3] * currentPan[5];
			}
			break;

		case 2: // 2.0 -> 7.1: create full surround field
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				for (int ch = 0; ch < 8; ch++)
				{
					currentPan[ch] += panDelta[ch];
				}

				float left = aScratch[sample];
				float right = aScratch[aBufferSize + sample];

				// FL, FR, C (mix), LFE (mix), SL, SR, BL, BR
				aBuffer[sample] += left * currentPan[0];
				aBuffer[aBufferSize + sample] += right * currentPan[1];
				aBuffer[aBufferSize * 2 + sample] += SURROUND_MIX_SCALE * (left + right) * currentPan[2];
				aBuffer[aBufferSize * 3 + sample] += SURROUND_MIX_SCALE * (left + right) * currentPan[3];
				aBuffer[aBufferSize * 4 + sample] += left * currentPan[4];
				aBuffer[aBufferSize * 5 + sample] += right * currentPan[5];
				aBuffer[aBufferSize * 6 + sample] += left * currentPan[6];
				aBuffer[aBufferSize * 7 + sample] += right * currentPan[7];
			}
			break;

		case 1: // 1.0 -> 7.1: distribute mono to all channels
			for (unsigned int sample = 0; sample < aSamplesToRead; sample++)
			{
				for (int ch = 0; ch < 8; ch++)
				{
					currentPan[ch] += panDelta[ch];
				}

				float monoSrc = aScratch[sample];

				for (int ch = 0; ch < 8; ch++)
				{
					aBuffer[aBufferSize * ch + sample] += monoSrc * currentPan[ch];
				}
			}
			break;
		}
		break;
	}

	// Update voice state with final channel volumes
	for (unsigned int ch = 0; ch < aChannels; ch++)
	{
		aVoice->mCurrentChannelVolume[ch] = targetPan[ch];
	}
}

void Mixer::handleResampleRemainder(unsigned int i, float *src, float *dst, unsigned int outputSamples, double srcPosition, double stepSize,
                                    unsigned int availableInput, unsigned int lookahead)
{
	using namespace ResamplingConstants;

	double srcPos = srcPosition + i * stepSize;
	for (; i < outputSamples; i++)
	{
		unsigned int intPos = (unsigned int)srcPos;
		float frac = (float)(srcPos - intPos);

		float sample;
		if (lookahead == POINT_LOOKAHEAD)
		{
			sample = src[intPos];
		}
		else if (lookahead == LINEAR_LOOKAHEAD)
		{
			float s0 = src[intPos];
			float s1 = (intPos + 1 < availableInput) ? src[intPos + 1] : s0;
			sample = s0 + (s1 - s0) * frac;
		}
		else // CATMULLROM_LOOKAHEAD
		{
			float p0 = (intPos >= 1) ? src[intPos - 1] : src[intPos];
			float p1 = src[intPos];
			float p2 = (intPos + 1 < availableInput) ? src[intPos + 1] : p1;
			float p3 = (intPos + 2 < availableInput) ? src[intPos + 2] : p1;

			float t = frac;
			float t2 = t * t;
			float t3 = t2 * t;

			sample = CATMULLROM_SCALE *
			         ((CATMULLROM_LINEAR_COEFF * p1) + (-p0 + p2) * t +
			          (CATMULLROM_QUAD_COEFFS[0] * p0 + CATMULLROM_QUAD_COEFFS[1] * p1 + CATMULLROM_QUAD_COEFFS[2] * p2 + CATMULLROM_QUAD_COEFFS[3] * p3) * t2 +
			          (CATMULLROM_CUBIC_COEFFS[0] * p0 + CATMULLROM_CUBIC_COEFFS[1] * p1 + CATMULLROM_CUBIC_COEFFS[2] * p2 + CATMULLROM_CUBIC_COEFFS[3] * p3) * t3);
		}

		dst[i] = sample;
		srcPos += stepSize;
	}
}

void Mixer::resample_channels(float **srcChannels, float *outputBuffer, unsigned int outputStride, unsigned int numChannels, unsigned int outputSamples,
                              double srcPosition, double stepSize, unsigned int availableInput, unsigned int lookahead)
{
	using namespace ResamplingConstants;

	SOLOUD_ASSERT(srcChannels != nullptr);
	SOLOUD_ASSERT(outputBuffer != nullptr);
	SOLOUD_ASSERT(numChannels > 0);
	SOLOUD_ASSERT(availableInput >= lookahead);

	for (unsigned int ch = 0; ch < numChannels; ch++)
	{
		float *src = srcChannels[ch];
		float *dst = outputBuffer + ch * outputStride;

		Mixer::handleResampleRemainder(0, src, dst, outputSamples, srcPosition, stepSize, availableInput, lookahead);
	}

	return;
}
} // namespace SoLoud::mixing

// NOLINTEND(cppcoreguidelines-avoid-c-arrays, modernize-avoid-c-arrays, cppcoreguidelines-pro-type-member-init, hicpp-member-init)
