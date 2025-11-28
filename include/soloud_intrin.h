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

#ifndef SOLOUD_INTRIN_H
#define SOLOUD_INTRIN_H

#include <array>

#include "soloud_config.h"

#ifdef SOLOUD_AVX_INTRINSICS
#include <immintrin.h>
#endif

#ifdef SOLOUD_SSE_INTRINSICS
#include <emmintrin.h>
#include <xmmintrin.h>
#endif

#include <cstddef>

namespace SoLoud
{
// SIMD configuration constants - adjust based on available instruction sets
#ifdef SOLOUD_AVX_INTRINSICS
// AVX2 configuration
static constexpr size_t SIMD_ALIGNMENT_BYTES = 32;    // 32-byte alignment for AVX2
static constexpr size_t SIMD_ALIGNMENT_MASK = 31;     // Mask for 32-byte alignment check
static constexpr size_t OPTIMAL_CHUNK_SAMPLES = 8;    // 8 floats per AVX2 register
static constexpr size_t MEMORY_ALIGNMENT_MASK = 7;    // Align to 8-float boundaries
static constexpr size_t TINY_BUFFER_FLOAT_COUNT = 32; // Buffer space for temporary operations
#elif defined(SOLOUD_SSE_INTRINSICS)
// SSE configuration
static constexpr size_t SIMD_ALIGNMENT_BYTES = 16;    // 16-byte alignment for SSE
static constexpr size_t SIMD_ALIGNMENT_MASK = 15;     // Mask for 16-byte alignment check (0xf)
static constexpr size_t OPTIMAL_CHUNK_SAMPLES = 4;    // 4 floats per SSE register
static constexpr size_t MEMORY_ALIGNMENT_MASK = 3;    // Align to 4-float boundaries
static constexpr size_t TINY_BUFFER_FLOAT_COUNT = 16; // Buffer space for temporary operations
#else
// Scalar fallback configuration
static constexpr size_t SIMD_ALIGNMENT_BYTES = 4;     // Basic float alignment
static constexpr size_t SIMD_ALIGNMENT_MASK = 3;      // Mask for 4-byte alignment check
static constexpr size_t OPTIMAL_CHUNK_SAMPLES = 4;    // Still process 4 samples at a time even for scalar
static constexpr size_t MEMORY_ALIGNMENT_MASK = 3;    // Align to 4-byte boundaries
static constexpr size_t TINY_BUFFER_FLOAT_COUNT = 16; // Buffer space for temporary operations (keeping this the same as SSE, this is how it was before)
#endif

// Resampling algorithm constants and helper functions
namespace ResamplingConstants
{
// Lookahead samples required for each resampling algorithm
constexpr unsigned int POINT_LOOKAHEAD = 1;      // Point sampling: just current sample
constexpr unsigned int LINEAR_LOOKAHEAD = 2;     // Linear: current + next sample
constexpr unsigned int CATMULLROM_LOOKAHEAD = 4; // Catmull-Rom: 4-point cubic interpolation

// Safety margin for high sample rate ratios to prevent buffer underruns
constexpr unsigned int LOOKAHEAD_SAFETY_MARGIN = 8;

// Catmull-Rom cubic interpolation coefficients
// Formula: 0.5 * ((2*p1) + (-p0+p2)*t + (2*p0-5*p1+4*p2-p3)*t^2 + (-p0+3*p1-3*p2+p3)*t^3)
constexpr float CATMULLROM_SCALE = 0.5f;
constexpr float CATMULLROM_LINEAR_COEFF = 2.0f;
constexpr float CATMULLROM_QUAD_COEFFS[4] = {2.0f, -5.0f, 4.0f, -1.0f};  // p0, p1, p2, p3 coefficients
constexpr float CATMULLROM_CUBIC_COEFFS[4] = {-1.0f, 3.0f, -3.0f, 1.0f}; // p0, p1, p2, p3 coefficients
} // namespace ResamplingConstants

// Lightweight class that handles small aligned buffer to support vectorized operations
// Used for temporary SIMD register-sized data
class TinyAlignedFloatBuffer
{
public:
	float *mData; // SIMD-aligned pointer
	std::array<unsigned char, sizeof(float) * TINY_BUFFER_FLOAT_COUNT + SIMD_ALIGNMENT_BYTES>
	    mActualData; // Space for appropriate number of floats + alignment padding

	// Constructor - automatically aligns mData to proper boundary
	TinyAlignedFloatBuffer();
};

class AudioSourceInstance;
class Soloud;

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
void resample_channels(float **srcChannels, float *outputBuffer, unsigned int outputStride, unsigned int numChannels, unsigned int outputSamples, double srcPosition,
                       double stepSize, unsigned int availableInput, unsigned int lookahead);

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
void clip_samples(const float *const aBuffer, float *aDestBuffer, unsigned int aSamples, unsigned int aChannels, float aVolume0, float aVolume1, float aScaler,
                  bool aRoundoff);

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
void interlace_samples(void *outputBuffer, const float *const rawBuffer, unsigned int aSamples, unsigned int stride, unsigned int aChannels,
                       detail::SAMPLE_FORMAT aFormat = detail::SAMPLE_FLOAT32);

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
void panAndExpand(AudioSourceInstance *aVoice, float *aBuffer, unsigned int aSamplesToRead, unsigned int aBufferSize, float *aScratch, unsigned int aChannels);
} // namespace SoLoud

#endif
