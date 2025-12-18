/*
SoLoud audio engine - logging facilities
Copyright (c) 2025 William Horvath

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

#include "soloud_config.h"

#include <cstdarg>
#include <cstdio>
#include <vector>

namespace SoLoud
{
namespace // anonymous
{
logFunctionType customStdoutFunction{nullptr};
void *customStdoutUserdata{nullptr};

logFunctionType customStderrFunction{nullptr};
void *customStderrUserdata{nullptr};
} // namespace

void setStdoutLogFunction(logFunctionType aCustomStdout, void *aUserdata)
{
	customStdoutFunction = aCustomStdout;
	customStdoutUserdata = aUserdata;
}

void setStderrLogFunction(logFunctionType aCustomStderr, void *aUserdata)
{
	customStderrFunction = aCustomStderr;
	customStderrUserdata = aUserdata;
}

void logStdout(const char *fmt, ...)
{
	va_list varargs;
	va_start(varargs, fmt);
	if (customStdoutFunction) // DISCLAIMER: I have no idea what I'm doing here.
	{
		va_list args_copy;
		va_copy(args_copy, varargs);

		int len = vsnprintf(nullptr, 0, fmt, varargs); // measure

		std::vector<char> buffer(len + 1);
		vsnprintf(buffer.data(), buffer.size(), fmt, args_copy);

		va_end(args_copy);
		customStdoutFunction(buffer.data(), customStdoutUserdata);
	}
	else
	{
		vfprintf(stdout, fmt, varargs);
	}
	va_end(varargs);
}

void logStderr(const char *fmt, ...)
{
	va_list varargs;
	va_start(varargs, fmt);
	if (customStderrFunction)
	{
		va_list args_copy;
		va_copy(args_copy, varargs);

		int len = vsnprintf(nullptr, 0, fmt, varargs); // measure

		std::vector<char> buffer(len + 1);
		vsnprintf(buffer.data(), buffer.size(), fmt, args_copy);

		va_end(args_copy);
		customStderrFunction(buffer.data(), customStderrUserdata);
	}
	else
	{
		vfprintf(stderr, fmt, varargs);
	}
	va_end(varargs);
}

} // namespace SoLoud
