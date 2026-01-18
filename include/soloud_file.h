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

#ifndef SOLOUD_FILE_H
#define SOLOUD_FILE_H

#include "soloud_config.h"
#include <cstdio>
#include <memory> // unique_ptr

typedef void *Soloud_Filehack;

namespace SoLoud
{
class File
{
public:
	File() = default;
	virtual ~File() = default;

	File(const File &) = default;
	File &operator=(const File &) = default;
	File(File &&) noexcept = default;
	File &operator=(File &&) noexcept = default;

	unsigned int read8();
	unsigned int read16();
	unsigned int read32();
	virtual bool eof() = 0;
	virtual unsigned int read(unsigned char *aDst, unsigned int aBytes) = 0;
	virtual unsigned int length() = 0;
	virtual void seek(int aOffset) = 0;
	virtual unsigned int pos() = 0;
	virtual FILE *getFilePtr() { return nullptr; }
	virtual const unsigned char *getMemPtr() { return nullptr; }
	[[nodiscard]] virtual const char *getFileName() const { return nullptr; }

protected:
	unsigned int mLength{0};
};

class DiskFile : public File
{
public:
	DiskFile() = default;
	DiskFile(FILE *fp)
	    : File(),
	      mFileHandle(fp)
	{
		// Sentinel value, for caching length.
		mLength = (unsigned int)-1;
	}

	// Construct with a filename directly
	DiskFile(const char *aFilename);
	~DiskFile() override;

	// Valid if getFileName is non-null for File.
	DiskFile(const File &);

	DiskFile(const DiskFile &);
	DiskFile &operator=(const DiskFile &);
	DiskFile(DiskFile &&) noexcept;
	DiskFile &operator=(DiskFile &&) noexcept;

	// Overridden methods.
	bool eof() override;
	unsigned int read(unsigned char *aDst, unsigned int aBytes) override;
	unsigned int length() override;
	void seek(int aOffset) override;
	unsigned int pos() override;
	FILE *getFilePtr() override { return mFileHandle; }
	[[nodiscard]] const char *getFileName() const override { return mFileName.get(); }

	// Methods unique to DiskFile.
	result open(const char *aFilename);

protected:
	FILE *mFileHandle{nullptr};
	std::unique_ptr<char[]> mFileName{nullptr};
	static FILE *openWithConversion(const char *const aFilenameSrc, char *&aFilenameOut);
};

class MemoryFile : public File
{
public:
	MemoryFile() = default;
	~MemoryFile() override;

	MemoryFile(const MemoryFile &);
	MemoryFile &operator=(const MemoryFile &);
	MemoryFile(MemoryFile &&) noexcept;
	MemoryFile &operator=(MemoryFile &&) noexcept;

	// Overridden methods.
	bool eof() override;
	unsigned int read(unsigned char *aDst, unsigned int aBytes) override;
	unsigned int length() override { return mLength; }
	void seek(int aOffset) override;
	unsigned int pos() override { return mOffset; }
	const unsigned char *getMemPtr() override { return mDataPtr; }

	// Methods unique to MemoryFile.
	result openMem(const unsigned char *aData, unsigned int aDataLength, bool aCopy = false, bool aTakeOwnership = true);
	result openToMem(const char *aFilename);
	result openFileToMem(File *aFile);

protected:
	const unsigned char *mDataPtr{nullptr};
	unsigned int mOffset{0};
	bool mDataOwned{false};
};
}; // namespace SoLoud

#endif
