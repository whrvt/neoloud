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

#include <cstdio> // for off_t, size_t, ssize_t
#ifdef _MSC_VER
#include <sys/types.h>
#endif

namespace SoLoud
{
class File;
}

#ifdef WITH_LIBMPG123

struct mpg123_handle_struct;
typedef struct mpg123_handle_struct mpg123_handle;

namespace SoLoud::MPG123
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
	bool ended;
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
bool isAtEnd(MPG123Decoder *aDecoder);
} // namespace SoLoud::MPG123
#else

namespace SoLoud::MPG123
{
// clang-format off
	struct MPG123Decoder{};
	inline MPG123Decoder *open(File * /*aFile*/) { return nullptr; }
	inline void close(MPG123Decoder * /*aDecoder*/) { ; }
	inline int getChannels(MPG123Decoder * /*aDecoder*/) { return -1; }
	inline int getSampleRate(MPG123Decoder * /*aDecoder*/) { return -1; }
	inline off_t getTotalFrameCount(MPG123Decoder * /*aDecoder*/) { return -1; }
	inline size_t readFrames(MPG123Decoder * /*aDecoder*/, size_t /*aFrameCount*/, float * /*aBuffer*/) { return 0; }
	inline off_t seekToFrame(MPG123Decoder * /*aDecoder*/, off_t /*aFrame*/) { return -1; }
	inline off_t getCurrentFrame(MPG123Decoder * /*aDecoder*/) { return -1; }
	inline bool isAtEnd(MPG123Decoder * /*aDecoder*/) { return true; }
	//clang-format on
} // namespace SoLoud::MPG123
#endif // WITH_LIBMPG123
#endif
