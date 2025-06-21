/*
SoLoud audio engine
Copyright (c) 2013-2018 Jari Komppa
Copyright (c) 2025 William Horvath (mpg123 interface)

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

#ifndef SOLOUD_MPG123_H
#define SOLOUD_MPG123_H

#include <stdio.h> // for off_t

struct mpg123_handle_struct;
typedef struct mpg123_handle_struct mpg123_handle;

namespace SoLoud
{
class File;
namespace MPG123
{
struct MPG123Decoder
{
	mpg123_handle *handle;
	File *file;
	int channels;
	long rate;
	off_t totalFrames;
	unsigned char *tempBuffer;
	size_t tempBufferSize;
};

// decoder functions
MPG123Decoder *open(File *aFile);
void close(MPG123Decoder *aDecoder);
int getChannels(MPG123Decoder *aDecoder);
int getSampleRate(MPG123Decoder *aDecoder);
off_t getTotalFrameCount(MPG123Decoder *aDecoder);
size_t readFrames(MPG123Decoder *aDecoder, size_t aFrameCount, float *aBuffer);
off_t seekToFrame(MPG123Decoder *aDecoder, off_t aFrame);
off_t getCurrentFrame(MPG123Decoder *aDecoder);

// custom I/O callbacks for mpg123
ssize_t readCallback(void *handle, void *buf, size_t count);
off_t seekCallback(void *handle, off_t offset, int whence);
} // namespace MPG123
} // namespace SoLoud
#endif
