/*
SoLoud audio engine
Copyright (c) 2013-2018 Jari Komppa

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

#ifndef SOLOUD_WAVSTREAM_H
#define SOLOUD_WAVSTREAM_H

#include "soloud_audiosource.h"
#include <memory> // unique_ptr

struct stb_vorbis;
#ifndef dr_flac_h
struct drflac;
#endif
#ifndef dr_mp3_h
struct drmp3;
struct drmp3_seek_point;
typedef unsigned int drmp3_uint32;
#endif
#ifndef dr_wav_h
struct drwav;
#endif

namespace SoLoud
{
namespace FFmpeg
{
struct FFmpegDecoder;
}
namespace MPG123
{
struct MPG123Decoder;
}
class WavStream;
class File;

enum WAVSTREAM_FILETYPE
{
	WAVSTREAM_WAV = 0,
	WAVSTREAM_OGG = 1,
	WAVSTREAM_FLAC = 2,
	WAVSTREAM_MPG123 = 3,
	WAVSTREAM_DRMP3 = 4,
	WAVSTREAM_FFMPEG = 5
};

class WavStreamInstance : public AudioSourceInstance
{
	struct FileInstanceHandle : public std::unique_ptr<File, void (*)(File *)>
	{
		using unique_ptr::unique_ptr;
		operator File *() const { return this->get(); }

		static void fileDeleter(File *rawFile);
		static void noopFileDeleter(File * /* rawFile */) {}
	};

	WavStream *mParent;
	unsigned int mOffset;
	FileInstanceHandle mFile;
	union {
		stb_vorbis *mOgg;
		drflac *mFlac;
		MPG123::MPG123Decoder *mMpg123;
		drmp3 *mDrmp3;
		drwav *mWav;
		FFmpeg::FFmpegDecoder *mFfmpeg;
		// default nullptr init, not to be used
		void *dummy{nullptr};
	};
	bool haveCodec() const { return mDrmp3 || mFfmpeg || mFlac || mMpg123 || mOgg || mWav; }

	unsigned int mOggFrameSize;
	unsigned int mOggFrameOffset;
	float **mOggOutputs;

public:
	WavStreamInstance() = delete;
	WavStreamInstance(WavStream *aParent);
	~WavStreamInstance() override;

	WavStreamInstance(const WavStreamInstance &) = delete;
	WavStreamInstance &operator=(const WavStreamInstance &) = delete;
	WavStreamInstance(WavStreamInstance &&) = delete;
	WavStreamInstance &operator=(WavStreamInstance &&) = delete;

	unsigned int getAudio(float *aBuffer, unsigned int aSamplesToRead, unsigned int aBufferSize) override;
	result seek(double aSeconds, float *mScratch, unsigned int mScratchSize) override;
	result rewind() override;
	bool hasEnded() override;
};

class WavStream : public AudioSource
{
	result loadwav(File *fp);
	result loadogg(File *fp);
	result loadflac(File *fp);
	result loadmpg123(File *fp);
	result loaddrmp3(File *fp);
	result loadffmpeg(File *fp);

	// mp3 seek tables
	friend class WavStreamInstance;
	std::unique_ptr<drmp3_seek_point[]> mMp3SeekPoints;
	drmp3_uint32 mMp3SeekPointCount;

	bool mPreferFFmpeg;

public:
	WAVSTREAM_FILETYPE mFiletype;
	std::unique_ptr<char[]> mFilename;
	std::unique_ptr<File> mMemFile;
	File *mStreamFile;
	unsigned int mSampleCount;

	WavStream(bool preferFFmpeg = false);
	~WavStream() override;

	WavStream(const WavStream &) = delete;
	WavStream &operator=(const WavStream &) = delete;
	WavStream(WavStream &&) = delete;
	WavStream &operator=(WavStream &&) = delete;

	result load(const char *aFilename);
	result loadMem(const unsigned char *aData, unsigned int aDataLen, bool aCopy = false, bool aTakeOwnership = true);
	result loadToMem(const char *aFilename);
	result loadFile(File *aFile);
	result loadFileToMem(File *aFile);
	AudioSourceInstance *createInstance() override;
	time getLength();

public:
	result parse(File *aFile);
};
}; // namespace SoLoud

#endif
