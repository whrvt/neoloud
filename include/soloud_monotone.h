/*
MONOTONE module for SoLoud audio engine
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

#ifndef MONOTONE_H
#define MONOTONE_H

#include "soloud_audiosource.h"

namespace SoLoud
{
class Monotone;
class File;

struct MonotoneSong
{
	char *mTitle;
	char *mComment;
	unsigned char mVersion; // must be 1
	unsigned char mTotalPatterns;
	unsigned char mTotalTracks;
	unsigned char mCellSize; // must be 2 for version 1
	std::array<unsigned char, 256> mOrder;
	unsigned int *mPatternData; // 64 rows * mTotalPatterns * mTotalTracks
};

struct MonotoneChannel
{
	int mEnabled;
	int mActive;
	std::array<int, 3> mFreq;
	int mPortamento;
	int mArpCounter;
	int mArp;
	int mLastNote;
	int mPortamentoToNote;
	int mVibrato;
	int mVibratoIndex;
	int mVibratoDepth;
	int mVibratoSpeed;
};

struct MonotoneHardwareChannel
{
	int mEnabled;
	float mSamplePos;
	float mSamplePosInc;
};

class MonotoneInstance : public AudioSourceInstance
{
	Monotone *mParent;

public:
	std::array<MonotoneChannel, 12> mChannel;
	std::array<MonotoneHardwareChannel, 12> mOutput;
	int mNextChannel;
	int mTempo; // ticks / row. Tick = 60hz. Default 4.
	int mOrder;
	int mRow;
	int mSampleCount;
	int mRowTick;

	MonotoneInstance(Monotone *aParent);
	virtual unsigned int getAudio(float *aBuffer, unsigned int aSamples, unsigned int aBufferSize);
	virtual bool hasEnded();
};

class Monotone : public AudioSource
{
public:
	std::array<int, 800> mNotesHz;
	std::array<int, 32> mVibTable;
	int mHardwareChannels;
	int mWaveform;
	MonotoneSong mSong;
	Monotone();
	~Monotone();
	result setParams(int aHardwareChannels, int aWaveform = 0); /* WAVE_SQUARE */
	result load(const char *aFilename);
	result loadMem(const unsigned char *aMem, unsigned int aLength, bool aCopy = false, bool aTakeOwnership = true);
	result loadFile(File *aFile);
	virtual AudioSourceInstance *createInstance();

public:
	void clear();
};
}; // namespace SoLoud

#endif