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

#include "sanity.h"

#include <cstdarg>

// global definitions
int gErrorCount = 0;
int gTests = 0;
int gVerbose = 2;

float gLastKnownScratch[2048]{};
FILE *gLastKnownFile = nullptr;
int gLastKnownWrite = 0;

#if defined(WINDOWS_VERSION)
#include <windows.h>

long getmsec()
{
	LARGE_INTEGER ts, freq;

	QueryPerformanceCounter(&ts);
	QueryPerformanceFrequency(&freq);
	ts.QuadPart *= 1000;
	ts.QuadPart /= freq.QuadPart;

	return (long)ts.QuadPart;
}
#else
#include <ctime>

long getmsec()
{
	struct timespec ts;

	clock_gettime(CLOCK_MONOTONIC, &ts);

	return (ts.tv_sec * 1000) + (ts.tv_nsec / 1000000);
}
#endif

static void writeHeader()
{
	alignas(int) unsigned char buf[46] = {
	    0x52, 0x49, 0x46,
	    0x46, // RIFF
	    0xa4, 0x3e, 0x00,
	    0x00, // length of file - 8
	    0x57, 0x41, 0x56,
	    0x45, // WAVE
	    0x66, 0x6d, 0x74,
	    0x20, // fmt
	    0x12, 0x00, 0x00,
	    0x00,       // PCM block
	    0x03, 0x00, // Uncompressed PCM, float
	    0x02, 0x00, // Channels (stereo)
	    0x44, 0xac, 0x00,
	    0x00, // 44100Hz
	    0x20, 0x62, 0x05,
	    0x00,       // 44100*8=352800 bytes per sec
	    0x08, 0x00, // Number of bytes per sample for all channels
	    0x20, 0x00, // Bits per channel (32)
	    0x00, 0x00, // 0 bytes of extension
	    0x64, 0x61, 0x74,
	    0x61, // data
	    0x80, 0x3e, 0x00,
	    0x00 // bytes of data
	         // 46 bytes up to this point
	};
	unsigned char *flen = buf + 4;
	unsigned char *dlen = buf + 42;
	if (!gLastKnownFile || !gLastKnownWrite)
		return;
	int len = ftell(gLastKnownFile);

	int temp = len - 8;
	memcpy(flen, &temp, sizeof(int));

	temp = len - 46;
	memcpy(dlen, &temp, sizeof(int));

	fseek(gLastKnownFile, 0, SEEK_SET);
	fwrite(buf, 1, 46, gLastKnownFile);
}

void generateTestWave(SoLoud::Wav &aWav)
{
	alignas(int) unsigned char charbuf[16044] = {
	    0x52, 0x49, 0x46, 0x46, // RIFF
	    0xa4, 0x3e, 0x00, 0x00, // length of file - 8
	    0x57, 0x41, 0x56, 0x45, // WAVE
	    0x66, 0x6d, 0x74, 0x20, // fmt
	    0x10, 0x00, 0x00, 0x00, // PCM block
	    0x01, 0x00,             // Uncompressed PCM
	    0x01, 0x00,             // Channels (mono)
	    0x40, 0x1f, 0x00, 0x00, // 8000Hz
	    0x40, 0x1f, 0x00, 0x00, // 8000 bytes per sec
	    0x01, 0x00,             // Number of bytes per sample for all channels
	    0x08, 0x00,             // Bits per channel
	    0x64, 0x61, 0x74, 0x61, // data
	    0x80, 0x3e, 0x00, 0x00, // bytes of data
	                            // 44 bytes up to this point
	};
	float floatbuf[(sizeof(charbuf) / sizeof(float)) + sizeof(float)];
	memcpy(floatbuf, charbuf, sizeof(charbuf));

	unsigned int buflen = sizeof(charbuf);
	int i;
	for (i = 0; i < (16000 / sizeof(float)); i++)
	{
		floatbuf[i + (44 / sizeof(float))] = ((i & 1) ? 1 : -1) * ((std::sin(i * i * 0.000001f) * 0x7f) + i);
	}
	aWav.loadMem(reinterpret_cast<const unsigned char *>(floatbuf), buflen, true, false);
}

// generate a longer test wave with known precise sample count for timing tests
void generateTimingTestWave(SoLoud::Wav &aWav, unsigned int sampleCount, unsigned int sampleRate)
{
	// WAV header: 44 bytes for standard PCM
	unsigned int dataSize = sampleCount * 2; // 16-bit = 2 bytes per sample
	unsigned int fileSize = 44 + dataSize;

	auto *buf = new unsigned char[fileSize];

	// RIFF header
	buf[0] = 'R';
	buf[1] = 'I';
	buf[2] = 'F';
	buf[3] = 'F';
	unsigned int chunkSize = fileSize - 8;
	memcpy(buf + 4, &chunkSize, 4);
	buf[8] = 'W';
	buf[9] = 'A';
	buf[10] = 'V';
	buf[11] = 'E';

	// fmt subchunk
	buf[12] = 'f';
	buf[13] = 'm';
	buf[14] = 't';
	buf[15] = ' ';
	unsigned int fmtSize = 16;
	memcpy(buf + 16, &fmtSize, 4);
	unsigned short audioFormat = 1; // PCM
	memcpy(buf + 20, &audioFormat, 2);
	unsigned short numChannels = 1; // mono
	memcpy(buf + 22, &numChannels, 2);
	memcpy(buf + 24, &sampleRate, 4);
	unsigned int byteRate = sampleRate * 2; // sampleRate * numChannels * bytesPerSample
	memcpy(buf + 28, &byteRate, 4);
	unsigned short blockAlign = 2; // numChannels * bytesPerSample
	memcpy(buf + 32, &blockAlign, 2);
	unsigned short bitsPerSample = 16;
	memcpy(buf + 34, &bitsPerSample, 2);

	// data subchunk
	buf[36] = 'd';
	buf[37] = 'a';
	buf[38] = 't';
	buf[39] = 'a';
	memcpy(buf + 40, &dataSize, 4);

	// Generate samples: simple sine wave
	auto *samples = reinterpret_cast<short *>(buf + 44);
	for (unsigned int i = 0; i < sampleCount; i++)
	{
		float t = (float)i / (float)sampleRate;
		float sample = std::sin(2.0f * (float)M_PI * 440.0f * t);
		samples[i] = (short)(sample * 16000.0f);
	}

	aWav.loadMem(buf, fileSize, true, false);
	delete[] buf;
}

namespace
{
struct TestEntry
{
	const char *name;
	void (*func)();
	bool isBench;
	bool usesLastKnown;
};
} // namespace

static const TestEntry gTestTable[] = {
    {.name = "misc",           .func = testMisc,                    .isBench = false, .usesLastKnown = false},
    {.name = "getters",        .func = testGetters,                 .isBench = false, .usesLastKnown = false},
    {.name = "vis",            .func = testVis,                     .isBench = false, .usesLastKnown = false},
    {.name = "timing",         .func = testRelativePlaySpeedTiming, .isBench = false, .usesLastKnown = false},
    {.name = "play",           .func = testPlay,                    .isBench = false, .usesLastKnown = true },
    {.name = "3d",             .func = test3d,                      .isBench = false, .usesLastKnown = true },
    {.name = "filters",        .func = testFilters,                 .isBench = false, .usesLastKnown = true },
    {.name = "core",           .func = testCore,                    .isBench = false, .usesLastKnown = true },
    {.name = "speech",         .func = testSpeech,                  .isBench = false, .usesLastKnown = true },
    {.name = "loudness",       .func = testLoudness,                .isBench = false, .usesLastKnown = false},
    {.name = "mixer",          .func = testMixer,                   .isBench = false, .usesLastKnown = false},
    {.name = "bench",          .func = testSpeedThings,             .isBench = true,  .usesLastKnown = false},
    {.name = "bench-loudness", .func = benchLoudness,               .isBench = true,  .usesLastKnown = false},
};

static const int gTestCount = sizeof(gTestTable) / sizeof(gTestTable[0]);

static void printUsage(const char *argv0)
{
	SoLoud::logStdout("Usage: %s [options] [test names...]\n\n", argv0);
	SoLoud::logStdout("Options:\n");
	SoLoud::logStdout("    -t NAME    run only the named test (repeatable)\n");
	SoLoud::logStdout("    -b         include benchmarks (excluded by default)\n");
	SoLoud::logStdout("    -l         list available tests and exit\n");
	SoLoud::logStdout("    -v N       verbosity (0=quiet, 1=normal, 2=verbose; default 2)\n");
	SoLoud::logStdout("    -h         show usage\n");
}

int main(int parc, char **pars)
{
	bool includeBench = false;
	bool selective = false;
	bool selected[sizeof(gTestTable) / sizeof(gTestTable[0])]{};

	for (int a = 1; a < parc; a++)
	{
		if (strcmp(pars[a], "-h") == 0 || strcmp(pars[a], "--help") == 0)
		{
			printUsage(pars[0]);
			return 0;
		}
		else if (strcmp(pars[a], "-l") == 0)
		{
			SoLoud::logStdout("Available tests:\n");
			for (int i = 0; i < gTestCount; i++)
			{
				SoLoud::logStdout("  %-12s%s%s\n", gTestTable[i].name, gTestTable[i].isBench ? " [bench]" : "", gTestTable[i].usesLastKnown ? " [lastknown]" : "");
			}
			return 0;
		}
		else if (strcmp(pars[a], "-b") == 0)
		{
			includeBench = true;
		}
		else if (strcmp(pars[a], "-v") == 0)
		{
			if (a + 1 < parc)
			{
				gVerbose = atoi(pars[++a]);
			}
			else
			{
				SoLoud::logStdout("Error: -v requires an argument\n");
				return 1;
			}
		}
		else if (strcmp(pars[a], "-t") == 0)
		{
			if (a + 1 >= parc)
			{
				SoLoud::logStdout("Error: -t requires an argument\n");
				return 1;
			}
			a++;
			selective = true;
			bool found = false;
			for (int i = 0; i < gTestCount; i++)
			{
				if (strcmp(pars[a], gTestTable[i].name) == 0)
				{
					selected[i] = true;
					found = true;
					break;
				}
			}
			if (!found)
			{
				SoLoud::logStdout("Error: unknown test '%s'\n", pars[a]);
				return 1;
			}
		}
		else
		{
			SoLoud::logStdout("Error: unknown option '%s'\n", pars[a]);
			printUsage(pars[0]);
			return 1;
		}
	}

	if (selective)
	{
		SoLoud::logStdout("Selective execution: lastknown.wav checks disabled.\n");
	}

#ifndef NO_LASTKNOWN_CHECK
	if (!selective)
	{
		gLastKnownFile = fopen("lastknown.wav", "rb");
		if (!gLastKnownFile)
		{
			SoLoud::logStdout("lastknown.wav not found, writing one.\n");
			gLastKnownFile = fopen("lastknown.wav", "wb");
			gLastKnownWrite = 1;
		}
		else
		{
			// skip header
			int i;
			for (i = 0; i < 46; i++)
				fgetc(gLastKnownFile);
		}
		writeHeader();
	}
#endif

	for (int i = 0; i < gTestCount; i++)
	{
		if (selective && !selected[i])
			continue;
		if (gTestTable[i].isBench && !includeBench)
			continue;

		gTestTable[i].func();
	}

	SoLoud::logStdout("\n%d tests, %d error(s) ", gTests, gErrorCount);
	if (!gLastKnownWrite && gErrorCount)
		SoLoud::logStdout("(To rebuild lastknown.wav, simply delete it)\n");
	SoLoud::logStdout("\n");
#ifndef NO_LASTKNOWN_CHECK
	if (!selective)
	{
		writeHeader();
		if (gLastKnownFile)
			fclose(gLastKnownFile);
		if (gLastKnownFile && gLastKnownWrite)
			SoLoud::logStdout("lastknown.wav written.\n");
	}
#endif
	return 0;
}
