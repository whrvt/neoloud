/*
SoLoud audio engine - ffmpeg interface
Copyright (c) 2013-2020 Jari Komppa
Copyright (c) 2025 William Horvath (ffmpeg interface)

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

#ifndef SOLOUD_FFMPEG_H
#define SOLOUD_FFMPEG_H

#if defined(WITH_FFMPEG) && __has_include(<libavcodec/avcodec.h>) && ((defined(_WIN32) || defined(_WIN64)) || defined(__linux__))

#include "soloud.h"

namespace SoLoud
{
class File;

namespace FFmpeg
{

// ffmpeg decoder handle (opaque)
struct FFmpegDecoder;

// load audio data to memory using ffmpeg
result loadToMemory(File *aFile, float **aData, unsigned int *aSampleCount, unsigned int *aChannels, float *aSampleRate);

// direct decoder interface (similar to dr_* decoders)
FFmpegDecoder *open(File *aFile, bool oneshot = false);
void close(FFmpegDecoder *decoder);
unsigned int getChannels(FFmpegDecoder *decoder);
unsigned int getSampleRate(FFmpegDecoder *decoder);
unsigned long long getTotalFrameCount(FFmpegDecoder *decoder);
unsigned long long getCurrentFrame(FFmpegDecoder *decoder);
bool seekToFrame(FFmpegDecoder *decoder, unsigned long long frameIndex);
unsigned long long readFrames(FFmpegDecoder *decoder, unsigned long long framesToRead, float *buffer);

} // namespace FFmpeg

} // namespace SoLoud

#else

namespace SoLoud
{
class File;
typedef unsigned int result;

namespace FFmpeg
{
struct FFmpegDecoder
{};

inline result loadToMemory(File *, float **, unsigned int *, unsigned int *, float *)
{
	return 7;
}

inline FFmpegDecoder *open(File *, bool = false)
{
	return nullptr;
}
inline void close(FFmpegDecoder *) {}
inline unsigned int getChannels(FFmpegDecoder *)
{
	return 0;
}
inline unsigned int getSampleRate(FFmpegDecoder *)
{
	return 0;
}
inline unsigned long long getTotalFrameCount(FFmpegDecoder *)
{
	return 0;
}
inline unsigned long long getCurrentFrame(FFmpegDecoder *)
{
	return 0;
}
inline bool seekToFrame(FFmpegDecoder *, unsigned long long)
{
	return false;
}
inline unsigned long long readFrames(FFmpegDecoder *, unsigned long long, float *)
{
	return 0;
}
} // namespace FFmpeg
} // namespace SoLoud

#endif
#endif // SOLOUD_FFMPEG_H
