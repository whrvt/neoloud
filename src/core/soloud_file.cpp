/*
SoLoud audio engine
Copyright (c) 2013-2015 Jari Komppa
Copyright (c) 2026 William Horvath

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

#undef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS

#include "soloud_file.h"
#include "soloud_error.h"

#include <cstdio>
#include <cstring>

#ifdef WINDOWS_VERSION
#include <windows.h>
#include <stringapiset.h>
#endif

namespace SoLoud
{
unsigned int File::read8()
{
	unsigned char d = 0;
	read(&d, 1);
	return d;
}

unsigned int File::read16()
{
	unsigned short d = 0;
	read((unsigned char *)&d, 2);
	return d;
}

unsigned int File::read32()
{
	unsigned int d = 0;
	read((unsigned char *)&d, 4);
	return d;
}

unsigned int DiskFile::read(unsigned char *aDst, unsigned int aBytes)
{
	return (unsigned int)fread(aDst, 1, aBytes, mFileHandle);
}

unsigned int DiskFile::length()
{
	if (!mFileHandle)
		return 0;

	// Return the cached value if we already have it,
	// which assumes that the FILE* is not being modified after initially opening it.
	// This should be a safe assumption for practical purposes.
	if (mLength == (unsigned int)-1)
	{
		// Save our current position.
		const auto pos = ftell(mFileHandle);
		// Errored, sets errno
		if (pos < 0)
			return 0;

		// Errored, should set errno
		if (fseek(mFileHandle, 0, SEEK_END) != 0)
			return 0;

		const auto len = ftell(mFileHandle);
		// Errored.
		if (len < 0)
			return 0;

		// Set the cached length
		mLength = (unsigned int)len;

		// Seek back to where we were
		(void)fseek(mFileHandle, pos, SEEK_SET);
	}

	return mLength;
}

void DiskFile::seek(int aOffset)
{
	fseek(mFileHandle, aOffset, SEEK_SET);
}

unsigned int DiskFile::pos()
{
	return (unsigned int)ftell(mFileHandle);
}

DiskFile::~DiskFile()
{
	if (mFileHandle)
		fclose(mFileHandle);
}

DiskFile::DiskFile(const char *aFilename)
    : File()
{
	mLength = (unsigned int)-1;
	open(aFilename);
}

FILE *DiskFile::openWithConversion(const char *const aFilenameSrc, char *&aFilenameOut)
{
	SOLOUD_ASSERT(aFilenameSrc);
	SOLOUD_ASSERT(!aFilenameOut);

	FILE *out = nullptr;
	const char *aFilename = (const char *)aFilenameSrc;
	// On Windows, convert to wide char internally to avoid failing on Unicode filenames.
#ifdef WINDOWS_VERSION
	std::unique_ptr<wchar_t[]> wideFilename;
	const int wideFilenameLen = MultiByteToWideChar(CP_UTF8, 0, aFilename, -1, nullptr, 0);
	if (wideFilenameLen > 0)
	{
		wideFilename = std::make_unique_for_overwrite<wchar_t[]>(wideFilenameLen);
		if (wideFilename && MultiByteToWideChar(CP_UTF8, 0, aFilename, -1, wideFilename.get(), wideFilenameLen) != 0)
		{
			out = _wfopen(wideFilename.get(), L"rb");
		}
	}
	if (out)
	{
		aFilenameOut = new char[wideFilenameLen * sizeof(wchar_t)];
		std::memcpy(aFilenameOut, wideFilename.get(), wideFilenameLen * sizeof(wchar_t));
	}
	else // Fall back to fopen, if anything failed for whatever reason.
#endif
	{
		out = fopen(aFilename, "rb");
		if (out)
		{
			const size_t len = std::strlen(aFilename);
			aFilenameOut = new char[len + 1];
			std::memcpy(aFilenameOut, aFilename, len);
			aFilenameOut[len] = '\0';
		}
	}

	return out;
}

DiskFile::DiskFile(const File &other)
    : File(other)
{
	mLength = (unsigned int)-1;

	if (!other.getFileName())
		return;

	char *tempName = nullptr;
	mFileHandle = openWithConversion(other.getFileName(), tempName);
	if (!mFileHandle)
		return;

	mFileName.reset(tempName);
}

DiskFile::DiskFile(const DiskFile &other)
    : File(other)
{
	mLength = (unsigned int)-1;

	if (!other.mFileName)
		return;

	char *tempName = nullptr;
	mFileHandle = openWithConversion(other.mFileName.get(), tempName);
	if (!mFileHandle)
		return;

	mFileName.reset(tempName);
}

DiskFile &DiskFile::operator=(const DiskFile &other)
{
	if (this == &other)
		return *this;

	// Clean up existing resources
	if (mFileHandle)
		fclose(mFileHandle);
	mFileHandle = nullptr;
	mFileName.reset();

	if (!other.mFileName)
	{
		mLength = (unsigned int)-1;
		return *this;
	}

	char *tempName = nullptr;
	mFileHandle = openWithConversion(other.mFileName.get(), tempName);
	if (!mFileHandle)
	{
		mLength = (unsigned int)-1;
		return *this;
	}

	mFileName.reset(tempName);
	mLength = (unsigned int)-1; // Re-cache on next length() call
	return *this;
}

DiskFile::DiskFile(DiskFile &&other) noexcept
    : mFileHandle(other.mFileHandle),
      mFileName(std::move(other.mFileName))
{
	mLength = other.mLength;

	other.mFileHandle = nullptr;
	other.mLength = (unsigned int)-1;
}

DiskFile &DiskFile::operator=(DiskFile &&other) noexcept
{
	if (this == &other)
		return *this;

	// Clean up existing resources
	if (mFileHandle)
		fclose(mFileHandle);

	mFileHandle = other.mFileHandle;
	mFileName = std::move(other.mFileName);
	mLength = other.mLength;

	other.mFileHandle = nullptr;
	other.mLength = (unsigned int)-1;
	return *this;
}

result DiskFile::open(const char *aFilename)
{
	if (!aFilename)
		return INVALID_PARAMETER;

	char *tempName = nullptr;
	mFileHandle = openWithConversion(aFilename, tempName);
	if (!mFileHandle)
		return FILE_NOT_FOUND;

	mFileName.reset(tempName);

	return SO_NO_ERROR;
}

bool DiskFile::eof()
{
	return feof(mFileHandle);
}

unsigned int MemoryFile::read(unsigned char *aDst, unsigned int aBytes)
{
	if (mOffset + aBytes >= mLength)
		aBytes = mLength - mOffset;

	memcpy(aDst, mDataPtr + mOffset, aBytes);
	mOffset += aBytes;

	return aBytes;
}

void MemoryFile::seek(int aOffset)
{
	if (aOffset >= 0)
		mOffset = aOffset;
	else
		mOffset = mLength + aOffset;
	if (mOffset > mLength - 1)
		mOffset = mLength - 1;
}

MemoryFile::~MemoryFile()
{
	if (mDataOwned)
		delete[] mDataPtr;
}

MemoryFile::MemoryFile(const MemoryFile &other)
    : File(other),
      mOffset(other.mOffset),
      mDataOwned(other.mDataOwned)
{
	if (other.mDataOwned && other.mDataPtr)
	{
		// Deep copy owned data
		unsigned char *newData = new unsigned char[mLength];
		std::memcpy(newData, other.mDataPtr, mLength);
		mDataPtr = newData;
	}
	else
	{
		// Non-owned: just copy the pointer (original owner manages lifetime)
		mDataPtr = other.mDataPtr;
	}
}

MemoryFile &MemoryFile::operator=(const MemoryFile &other)
{
	if (this == &other)
		return *this;

	// Clean up existing owned data
	if (mDataOwned)
		delete[] mDataPtr;

	mLength = other.mLength;
	mOffset = other.mOffset;
	mDataOwned = other.mDataOwned;

	if (other.mDataOwned && other.mDataPtr)
	{
		unsigned char *newData = new unsigned char[mLength];
		std::memcpy(newData, other.mDataPtr, mLength);
		mDataPtr = newData;
	}
	else
	{
		mDataPtr = other.mDataPtr;
	}

	return *this;
}

MemoryFile::MemoryFile(MemoryFile &&other) noexcept
    : File(other),
      mDataPtr(other.mDataPtr),
      mOffset(other.mOffset),
      mDataOwned(other.mDataOwned)
{
	other.mDataPtr = nullptr;
	other.mLength = 0;
	other.mOffset = 0;
	other.mDataOwned = false;
}

MemoryFile &MemoryFile::operator=(MemoryFile &&other) noexcept
{
	if (this == &other)
		return *this;

	// Clean up existing owned data
	if (mDataOwned)
		delete[] mDataPtr;

	mDataPtr = other.mDataPtr;
	mLength = other.mLength;
	mOffset = other.mOffset;
	mDataOwned = other.mDataOwned;

	other.mDataPtr = nullptr;
	other.mLength = 0;
	other.mOffset = 0;
	other.mDataOwned = false;

	return *this;
}

result MemoryFile::openMem(const unsigned char *aData, unsigned int aDataLength, bool aCopy, bool aTakeOwnership)
{
	if (aData == nullptr || aDataLength == 0)
		return INVALID_PARAMETER;

	if (mDataOwned)
		delete[] mDataPtr;
	mDataPtr = nullptr;
	mOffset = 0;

	mLength = aDataLength;

	if (aCopy)
	{
		mDataOwned = true;
		mDataPtr = new unsigned char[aDataLength];
		if (mDataPtr == nullptr)
			return OUT_OF_MEMORY;
		memcpy((void *)mDataPtr, aData, aDataLength);
		return SO_NO_ERROR;
	}

	mDataPtr = aData;
	mDataOwned = aTakeOwnership;
	return SO_NO_ERROR;
}

result MemoryFile::openToMem(const char *aFile)
{
	if (!aFile)
		return INVALID_PARAMETER;
	if (mDataOwned)
		delete[] mDataPtr;
	mDataPtr = nullptr;
	mOffset = 0;

	DiskFile df;
	result res = df.open(aFile);
	if (res != SO_NO_ERROR)
		return res;

	mLength = df.length();
	mDataPtr = new unsigned char[mLength];
	if (mDataPtr == nullptr)
		return OUT_OF_MEMORY;
	df.read((unsigned char *)mDataPtr, mLength);
	mDataOwned = true;
	return SO_NO_ERROR;
}

result MemoryFile::openFileToMem(File *aFile)
{
	if (!aFile)
		return INVALID_PARAMETER;
	if (mDataOwned)
		delete[] mDataPtr;
	mDataPtr = nullptr;
	mOffset = 0;

	mLength = aFile->length();
	mDataPtr = new unsigned char[mLength];
	if (mDataPtr == nullptr)
		return OUT_OF_MEMORY;
	aFile->read((unsigned char *)mDataPtr, mLength);
	mDataOwned = true;
	return SO_NO_ERROR;
}

bool MemoryFile::eof()
{
	if (mOffset >= mLength)
		return true;
	return false;
}
} // namespace SoLoud

extern "C"
{
	int Soloud_Filehack_fgetc(Soloud_Filehack *f)
	{
		SoLoud::File *fp = (SoLoud::File *)f;
		if (fp->eof())
			return EOF;
		return fp->read8();
	}

	int Soloud_Filehack_fread(void *dst, int s, int c, Soloud_Filehack *f)
	{
		SoLoud::File *fp = (SoLoud::File *)f;
		return fp->read((unsigned char *)dst, s * c) / s;
	}

	int Soloud_Filehack_fseek(Soloud_Filehack *f, int idx, int base)
	{
		SoLoud::File *fp = (SoLoud::File *)f;
		switch (base)
		{
		case SEEK_CUR:
			fp->seek(fp->pos() + idx);
			break;
		case SEEK_END:
			fp->seek(fp->length() + idx);
			break;
		default:
			fp->seek(idx);
		}
		return 0;
	}

	int Soloud_Filehack_ftell(Soloud_Filehack *f)
	{
		SoLoud::File *fp = (SoLoud::File *)f;
		return fp->pos();
	}

	int Soloud_Filehack_fclose(Soloud_Filehack *f)
	{
		SoLoud::File *fp = (SoLoud::File *)f;
		delete fp;
		return 0;
	}

	Soloud_Filehack *Soloud_Filehack_fopen(const char *aFilename, char * /*aMode*/)
	{
		SoLoud::DiskFile *df = new SoLoud::DiskFile();
		SoLoud::result res = df->open(aFilename);
		if (res != SoLoud::SO_NO_ERROR)
		{
			delete df;
			df = nullptr;
		}
		return (Soloud_Filehack *)df;
	}

	int Soloud_Filehack_fopen_s(Soloud_Filehack **f, const char *aFilename, char * /*aMode*/)
	{
		*f = Soloud_Filehack_fopen(aFilename, nullptr);
		return 1;
	}
}
