/*
SoLoud audio engine - low-level mixing interface (**mostly** internal)
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

#ifndef SOLOUD_LL_MIXING_H
#define SOLOUD_LL_MIXING_H

#include <cstddef>
namespace SoLoud
{

class AudioSourceInstance;
class Soloud;

// For use by backends to specify which format they'd like from the mixer.
namespace mixing
{
enum SAMPLE_FORMAT : unsigned char
{
	SAMPLE_FLOAT32,
	SAMPLE_UNSIGNED8,
	SAMPLE_SIGNED16,
	SAMPLE_SIGNED24,
	SAMPLE_SIGNED32
};

class Mixer
{
private:
	friend Soloud;
	static Mixer *createMixer();

protected:
	Mixer() = default;
	// For deduplication between implementations
	static void handleResampleRemainder(unsigned int i, float *src, float *dst, unsigned int outputSamples, double srcPosition, double stepSize,
	                                    unsigned int availableInput, unsigned int lookahead);

public:
	virtual ~Mixer() = default;

	Mixer(const Mixer &) = delete;
	Mixer &operator=(const Mixer &) = delete;
	Mixer(Mixer &&) = delete;
	Mixer &operator=(Mixer &&) = delete;

	/**
	 * Resample audio channels from source to destination sample rate
	 *
	 * @param srcChannels       Array of pointers to per-channel source buffers
	 * @param outputBuffer      Output buffer, channels separated by outputStride
	 * @param outputStride      Stride between channels in output buffer (samples)
	 * @param numChannels       Number of audio channels to process
	 * @param outputSamples     Number of samples to process and output
	 * @param srcPosition       Starting fractional position in source stream
	 * @param stepSize          Source position increment per output sample (srcRate/dstRate)
	 * @param availableInput    Number of available samples in source buffers
	 * @param lookahead         Lookahead samples required by the resampler algorithm
	 * @return                  outputSamples (always processes exactly the requested amount)
	 *
	 * Caller must ensure: availableInput >= floor(srcPosition + outputSamples * stepSize) + lookahead
	 *
	 * Uses AVX2 gather instructions when available, SSE with manual gathering as fallback,
	 * or scalar for compatibility.
	 */
	virtual void resample_channels(float **srcChannels, float *outputBuffer, unsigned int outputStride, unsigned int numChannels, unsigned int outputSamples,
	                               double srcPosition, double stepSize, unsigned int availableInput, unsigned int lookahead);

	/**
	 * Apply volume scaling and clipping to audio buffer
	 *
	 * @param aBuffer       Input buffer with source samples
	 * @param aDestBuffer   Output buffer for processed samples
	 * @param aSamples      Number of samples to process
	 * @param aChannels     Number of channels
	 * @param aVolume0      Starting volume level
	 * @param aVolume1      Ending volume level (for smooth ramping)
	 * @param aScaler       Post-clip scale factor
	 * @param aRoundoff     Whether to do roundoff clipping
	 *
	 * Supports two clipping modes:
	 * - Hard clipping: Simple [-1,1] bounds
	 * - Roundoff clipping: Smooth saturation curve that approaches limits asymptotically
	 *
	 * Uses AVX2 intrinsics when available for 8x performance improvement,
	 * falls back to SSE intrinsics for 4x performance improvement,
	 * or scalar fallback for compatibility
	 */
	virtual void clip_samples(const float *const aBuffer, float *aDestBuffer, unsigned int aSamples, unsigned int aChannels, float aVolume0, float aVolume1,
	                          float aScaler, bool aRoundoff);

	/**
	 * Convert channel-separated audio data to interleaved format with sample format conversion
	 *
	 * @param outputBuffer   Destination buffer for interleaved samples
	 * @param rawBuffer      Source buffer with channel-separated float data
	 * @param aSamples       Number of samples per channel to process
	 * @param stride         Size of each channel buffer (including padding for alignment)
	 * @param aChannels      Number of audio channels
	 * @param aFormat        Target sample format for conversion
	 *
	 * Converts from SoLoud's internal channel-separated float format to interleaved output
	 * in the specified sample format. Handles scaling, clamping, and type conversion for:
	 * - SAMPLE_FLOAT32: Direct copy (32-bit float)
	 * - SAMPLE_UNSIGNED8: Scale to [0,255] range (8-bit unsigned)
	 * - SAMPLE_SIGNED16: Scale to [-32767,32767] range (16-bit signed)
	 * - SAMPLE_SIGNED24: Scale to 24-bit range, store as 3-byte packed
	 * - SAMPLE_SIGNED32: Scale to full 32-bit signed integer range
	 *
	 */
	virtual void interlace_samples(void *outputBuffer, const float *const rawBuffer, unsigned int aSamples, unsigned int stride, unsigned int aChannels,
	                               SAMPLE_FORMAT aFormat = SAMPLE_FLOAT32);

	/**
	 * Pan and expand audio from source channel count to output channel count
	 *
	 * @param aVoice         Voice instance containing channel volume settings
	 * @param aBuffer        Output buffer (accumulative mixing)
	 * @param aSamplesToRead Number of samples to process
	 * @param aBufferSize    Size of each channel buffer (for stride calculation)
	 * @param aScratch       Input source data (separated by channel)
	 * @param aChannels      Number of output channels
	 *
	 * Handles all common channel configurations:
	 * - 1.0 (mono) -> 1.0, 2.0, 4.0, 5.1, 7.1
	 * - 2.0 (stereo) -> 1.0, 2.0, 4.0, 5.1, 7.1
	 * - 4.0 (quad) -> 1.0, 2.0, 4.0, 5.1, 7.1
	 * - 5.1 (surround) -> 1.0, 2.0, 4.0, 5.1, 7.1
	 * - 7.1 (surround) -> 1.0, 2.0, 4.0, 5.1, 7.1
	 *
	 * Provides smooth volume ramping to avoid audio artifacts.
	 * Uses intelligent downmixing coefficients to preserve perceived loudness.
	 * Optimized with AVX2 intrinsics for stereo output paths when available,
	 * falls back to SSE intrinsics or scalar implementation.
	 */
	virtual void panAndExpand(AudioSourceInstance *aVoice, float *aBuffer, unsigned int aSamplesToRead, unsigned int aBufferSize, float *aScratch,
	                          unsigned int aChannels);
};

} // namespace mixing

} // namespace SoLoud

#endif
