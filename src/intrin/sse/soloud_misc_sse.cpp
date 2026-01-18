/*
SoLoud audio engine - miscellaneous SSE-related things
Copyright (c) 2013-2020 Jari Komppa
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

#include "soloud_cpu.h"

#if defined(SOLOUD_SUPPORT_SSE2) || defined(SOLOUD_IS_X86_64)

#ifdef _MSC_VER
#include <intrin.h>
#endif
#include <xmmintrin.h>

namespace SoLoud
{
extern void setCTZDAZ();
void setCTZDAZ()
{
	_mm_setcsr(_mm_getcsr() | 0x8040);
}

} // namespace SoLoud

#endif