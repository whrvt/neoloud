
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

#include "soloud_intrin.h"
#include "soloud.h"
#include "soloud_audiosource.h"

#ifdef SOLOUD_SSE_INTRINSICS
#include <xmmintrin.h>

#include <cstddef>
#ifdef _M_IX86
#include <emmintrin.h>
#endif
#endif

namespace SoLoud
{

#if defined(SOLOUD_SSE_INTRINSICS)
void clip_internal(const Soloud *aSoloud, AlignedFloatBuffer &aBuffer, AlignedFloatBuffer &aDestBuffer, unsigned int aSamples, float aVolume0, float aVolume1)
{
	float vd = (aVolume1 - aVolume0) / aSamples;
	float v = aVolume0;
	unsigned int i, j, c, d;
	unsigned int samplequads = (aSamples + 3) / 4; // rounded up

	// Clip
	if (aSoloud->mFlags & Soloud::CLIP_ROUNDOFF)
	{
		float nb = -1.65f;
		__m128 negbound = _mm_load_ps1(&nb);
		float pb = 1.65f;
		__m128 posbound = _mm_load_ps1(&pb);
		float ls = 0.87f;
		__m128 linearscale = _mm_load_ps1(&ls);
		float cs = -0.1f;
		__m128 cubicscale = _mm_load_ps1(&cs);
		float nw = -0.9862875f;
		__m128 negwall = _mm_load_ps1(&nw);
		float pw = 0.9862875f;
		__m128 poswall = _mm_load_ps1(&pw);
		__m128 postscale = _mm_load_ps1(&aSoloud->mPostClipScaler);
		TinyAlignedFloatBuffer volumes;
		volumes.mData[0] = v;
		volumes.mData[1] = v + vd;
		volumes.mData[2] = v + vd + vd;
		volumes.mData[3] = v + vd + vd + vd;
		vd *= 4;
		__m128 vdelta = _mm_load_ps1(&vd);
		c = 0;
		d = 0;
		for (j = 0; j < aSoloud->mChannels; j++)
		{
			__m128 vol = _mm_load_ps(volumes.mData);

			for (i = 0; i < samplequads; i++)
			{
				// float f1 = origdata[c] * v;	c++; v += vd;
				__m128 f = _mm_load_ps(&aBuffer.mData[c]);
				c += 4;
				f = _mm_mul_ps(f, vol);
				vol = _mm_add_ps(vol, vdelta);

				// float u1 = (f1 > -1.65f);
				__m128 u = _mm_cmpgt_ps(f, negbound);

				// float o1 = (f1 < 1.65f);
				__m128 o = _mm_cmplt_ps(f, posbound);

				// f1 = (0.87f * f1 - 0.1f * f1 * f1 * f1) * u1 * o1;
				__m128 lin = _mm_mul_ps(f, linearscale);
				__m128 cubic = _mm_mul_ps(f, f);
				cubic = _mm_mul_ps(cubic, f);
				cubic = _mm_mul_ps(cubic, cubicscale);
				f = _mm_add_ps(cubic, lin);

				// f1 = f1 * u1 + !u1 * -0.9862875f;
				__m128 lowmask = _mm_andnot_ps(u, negwall);
				__m128 ilowmask = _mm_and_ps(u, f);
				f = _mm_add_ps(lowmask, ilowmask);

				// f1 = f1 * o1 + !o1 * 0.9862875f;
				__m128 himask = _mm_andnot_ps(o, poswall);
				__m128 ihimask = _mm_and_ps(o, f);
				f = _mm_add_ps(himask, ihimask);

				// outdata[d] = f1 * postclip; d++;
				f = _mm_mul_ps(f, postscale);
				_mm_store_ps(&aDestBuffer.mData[d], f);
				d += 4;
			}
		}
	}
	else
	{
		float nb = -1.0f;
		__m128 negbound = _mm_load_ps1(&nb);
		float pb = 1.0f;
		__m128 posbound = _mm_load_ps1(&pb);
		__m128 postscale = _mm_load_ps1(&aSoloud->mPostClipScaler);
		TinyAlignedFloatBuffer volumes;
		volumes.mData[0] = v;
		volumes.mData[1] = v + vd;
		volumes.mData[2] = v + vd + vd;
		volumes.mData[3] = v + vd + vd + vd;
		vd *= 4;
		__m128 vdelta = _mm_load_ps1(&vd);
		c = 0;
		d = 0;
		for (j = 0; j < aSoloud->mChannels; j++)
		{
			__m128 vol = _mm_load_ps(volumes.mData);
			for (i = 0; i < samplequads; i++)
			{
				// float f1 = aBuffer.mData[c] * v; c++; v += vd;
				__m128 f = _mm_load_ps(&aBuffer.mData[c]);
				c += 4;
				f = _mm_mul_ps(f, vol);
				vol = _mm_add_ps(vol, vdelta);

				// f1 = (f1 <= -1) ? -1 : (f1 >= 1) ? 1 : f1;
				f = _mm_max_ps(f, negbound);
				f = _mm_min_ps(f, posbound);

				// aDestBuffer.mData[d] = f1 * mPostClipScaler; d++;
				f = _mm_mul_ps(f, postscale);
				_mm_store_ps(&aDestBuffer.mData[d], f);
				d += 4;
			}
		}
	}
}
#else // fallback code
void clip_internal(const Soloud *aSoloud, AlignedFloatBuffer &aBuffer, AlignedFloatBuffer &aDestBuffer, unsigned int aSamples, float aVolume0, float aVolume1)
{
	float vd = (aVolume1 - aVolume0) / aSamples;
	float v = aVolume0;
	unsigned int i, j, c, d;
	unsigned int samplequads = (aSamples + 3) / 4; // rounded up
	// Clip
	if (aSoloud->mFlags & Soloud::CLIP_ROUNDOFF)
	{
		c = 0;
		d = 0;
		for (j = 0; j < aSoloud->mChannels; j++)
		{
			v = aVolume0;
			for (i = 0; i < samplequads; i++)
			{
				float f1 = aBuffer.mData[c] * v;
				c++;
				v += vd;
				float f2 = aBuffer.mData[c] * v;
				c++;
				v += vd;
				float f3 = aBuffer.mData[c] * v;
				c++;
				v += vd;
				float f4 = aBuffer.mData[c] * v;
				c++;
				v += vd;

				f1 = (f1 <= -1.65f) ? -0.9862875f : (f1 >= 1.65f) ? 0.9862875f : (0.87f * f1 - 0.1f * f1 * f1 * f1);
				f2 = (f2 <= -1.65f) ? -0.9862875f : (f2 >= 1.65f) ? 0.9862875f : (0.87f * f2 - 0.1f * f2 * f2 * f2);
				f3 = (f3 <= -1.65f) ? -0.9862875f : (f3 >= 1.65f) ? 0.9862875f : (0.87f * f3 - 0.1f * f3 * f3 * f3);
				f4 = (f4 <= -1.65f) ? -0.9862875f : (f4 >= 1.65f) ? 0.9862875f : (0.87f * f4 - 0.1f * f4 * f4 * f4);

				aDestBuffer.mData[d] = f1 * aSoloud->mPostClipScaler;
				d++;
				aDestBuffer.mData[d] = f2 * aSoloud->mPostClipScaler;
				d++;
				aDestBuffer.mData[d] = f3 * aSoloud->mPostClipScaler;
				d++;
				aDestBuffer.mData[d] = f4 * aSoloud->mPostClipScaler;
				d++;
			}
		}
	}
	else
	{
		c = 0;
		d = 0;
		for (j = 0; j < aSoloud->mChannels; j++)
		{
			v = aVolume0;
			for (i = 0; i < samplequads; i++)
			{
				float f1 = aBuffer.mData[c] * v;
				c++;
				v += vd;
				float f2 = aBuffer.mData[c] * v;
				c++;
				v += vd;
				float f3 = aBuffer.mData[c] * v;
				c++;
				v += vd;
				float f4 = aBuffer.mData[c] * v;
				c++;
				v += vd;

				f1 = (f1 <= -1) ? -1 : (f1 >= 1) ? 1 : f1;
				f2 = (f2 <= -1) ? -1 : (f2 >= 1) ? 1 : f2;
				f3 = (f3 <= -1) ? -1 : (f3 >= 1) ? 1 : f3;
				f4 = (f4 <= -1) ? -1 : (f4 >= 1) ? 1 : f4;

				aDestBuffer.mData[d] = f1 * aSoloud->mPostClipScaler;
				d++;
				aDestBuffer.mData[d] = f2 * aSoloud->mPostClipScaler;
				d++;
				aDestBuffer.mData[d] = f3 * aSoloud->mPostClipScaler;
				d++;
				aDestBuffer.mData[d] = f4 * aSoloud->mPostClipScaler;
				d++;
			}
		}
	}
}
#endif

void panAndExpand(AudioSourceInstance *aVoice, float *aBuffer, unsigned int aSamplesToRead, unsigned int aBufferSize, float *aScratch, unsigned int aChannels)
{
#ifdef SOLOUD_SSE_INTRINSICS
	SOLOUD_ASSERT(((size_t)aBuffer & 0xf) == 0);
	SOLOUD_ASSERT(((size_t)aScratch & 0xf) == 0);
	SOLOUD_ASSERT(((size_t)aBufferSize & 0xf) == 0);
#endif
	float pan[MAX_CHANNELS];  // current speaker volume
	float pand[MAX_CHANNELS]; // destination speaker volume
	float pani[MAX_CHANNELS]; // speaker volume increment per sample
	unsigned int j, k;
	for (k = 0; k < aChannels; k++)
	{
		pan[k] = aVoice->mCurrentChannelVolume[k];
		pand[k] = aVoice->mChannelVolume[k] * aVoice->mOverallVolume;
		pani[k] = (pand[k] - pan[k]) / aSamplesToRead; // TODO: this is a bit inconsistent.. but it's a hack to begin with
	}

	int ofs = 0;
	switch (aChannels)
	{
	case 1: // Target is mono. Sum everything. (1->1, 2->1, 4->1, 6->1, 8->1)
		for (j = 0, ofs = 0; j < aVoice->mChannels; j++, ofs += aBufferSize)
		{
			pan[0] = aVoice->mCurrentChannelVolume[0];
			for (k = 0; k < aSamplesToRead; k++)
			{
				pan[0] += pani[0];
				aBuffer[k] += aScratch[ofs + k] * pan[0];
			}
		}
		break;
	case 2:
		switch (aVoice->mChannels)
		{
		case 8: // 8->2, just sum lefties and righties, add a bit of center and sub?
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				float s1 = aScratch[j];
				float s2 = aScratch[aBufferSize + j];
				float s3 = aScratch[aBufferSize * 2 + j];
				float s4 = aScratch[aBufferSize * 3 + j];
				float s5 = aScratch[aBufferSize * 4 + j];
				float s6 = aScratch[aBufferSize * 5 + j];
				float s7 = aScratch[aBufferSize * 6 + j];
				float s8 = aScratch[aBufferSize * 7 + j];
				aBuffer[j + 0] += 0.2f * (s1 + s3 + s4 + s5 + s7) * pan[0];
				aBuffer[j + aBufferSize] += 0.2f * (s2 + s3 + s4 + s6 + s8) * pan[1];
			}
			break;
		case 6: // 6->2, just sum lefties and righties, add a bit of center and sub?
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				float s1 = aScratch[j];
				float s2 = aScratch[aBufferSize + j];
				float s3 = aScratch[aBufferSize * 2 + j];
				float s4 = aScratch[aBufferSize * 3 + j];
				float s5 = aScratch[aBufferSize * 4 + j];
				float s6 = aScratch[aBufferSize * 5 + j];
				aBuffer[j + 0] += 0.3f * (s1 + s3 + s4 + s5) * pan[0];
				aBuffer[j + aBufferSize] += 0.3f * (s2 + s3 + s4 + s6) * pan[1];
			}
			break;
		case 4: // 4->2, just sum lefties and righties
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				float s1 = aScratch[j];
				float s2 = aScratch[aBufferSize + j];
				float s3 = aScratch[aBufferSize * 2 + j];
				float s4 = aScratch[aBufferSize * 3 + j];
				aBuffer[j + 0] += 0.5f * (s1 + s3) * pan[0];
				aBuffer[j + aBufferSize] += 0.5f * (s2 + s4) * pan[1];
			}
			break;
		case 2: // 2->2
#if defined(SOLOUD_SSE_INTRINSICS)
		{
			int c = 0;
			// if ((aBufferSize & 3) == 0)
			{
				unsigned int samplequads = aSamplesToRead / 4; // rounded down
				TinyAlignedFloatBuffer pan0;
				pan0.mData[0] = pan[0] + pani[0];
				pan0.mData[1] = pan[0] + pani[0] * 2;
				pan0.mData[2] = pan[0] + pani[0] * 3;
				pan0.mData[3] = pan[0] + pani[0] * 4;
				TinyAlignedFloatBuffer pan1;
				pan1.mData[0] = pan[1] + pani[1];
				pan1.mData[1] = pan[1] + pani[1] * 2;
				pan1.mData[2] = pan[1] + pani[1] * 3;
				pan1.mData[3] = pan[1] + pani[1] * 4;
				pani[0] *= 4;
				pani[1] *= 4;
				__m128 pan0delta = _mm_load_ps1(&pani[0]);
				__m128 pan1delta = _mm_load_ps1(&pani[1]);
				__m128 p0 = _mm_load_ps(pan0.mData);
				__m128 p1 = _mm_load_ps(pan1.mData);

				for (j = 0; j < samplequads; j++)
				{
					__m128 f0 = _mm_load_ps(aScratch + c);
					__m128 c0 = _mm_mul_ps(f0, p0);
					__m128 f1 = _mm_load_ps(aScratch + c + aBufferSize);
					__m128 c1 = _mm_mul_ps(f1, p1);
					__m128 o0 = _mm_load_ps(aBuffer + c);
					__m128 o1 = _mm_load_ps(aBuffer + c + aBufferSize);
					c0 = _mm_add_ps(c0, o0);
					c1 = _mm_add_ps(c1, o1);
					_mm_store_ps(aBuffer + c, c0);
					_mm_store_ps(aBuffer + c + aBufferSize, c1);
					p0 = _mm_add_ps(p0, pan0delta);
					p1 = _mm_add_ps(p1, pan1delta);
					c += 4;
				}
			}

			// If buffer size or samples to read are not divisible by 4, handle leftovers
			for (j = c; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				float s1 = aScratch[j];
				float s2 = aScratch[aBufferSize + j];
				aBuffer[j + 0] += s1 * pan[0];
				aBuffer[j + aBufferSize] += s2 * pan[1];
			}
		}
#else // fallback
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				float s1 = aScratch[j];
				float s2 = aScratch[aBufferSize + j];
				aBuffer[j + 0] += s1 * pan[0];
				aBuffer[j + aBufferSize] += s2 * pan[1];
			}
#endif
		break;
		case 1: // 1->2
#if defined(SOLOUD_SSE_INTRINSICS)
		{
			int c = 0;
			// if ((aBufferSize & 3) == 0)
			{
				unsigned int samplequads = aSamplesToRead / 4; // rounded down
				TinyAlignedFloatBuffer pan0;
				pan0.mData[0] = pan[0] + pani[0];
				pan0.mData[1] = pan[0] + pani[0] * 2;
				pan0.mData[2] = pan[0] + pani[0] * 3;
				pan0.mData[3] = pan[0] + pani[0] * 4;
				TinyAlignedFloatBuffer pan1;
				pan1.mData[0] = pan[1] + pani[1];
				pan1.mData[1] = pan[1] + pani[1] * 2;
				pan1.mData[2] = pan[1] + pani[1] * 3;
				pan1.mData[3] = pan[1] + pani[1] * 4;
				pani[0] *= 4;
				pani[1] *= 4;
				__m128 pan0delta = _mm_load_ps1(&pani[0]);
				__m128 pan1delta = _mm_load_ps1(&pani[1]);
				__m128 p0 = _mm_load_ps(pan0.mData);
				__m128 p1 = _mm_load_ps(pan1.mData);

				for (j = 0; j < samplequads; j++)
				{
					__m128 f = _mm_load_ps(aScratch + c);
					__m128 c0 = _mm_mul_ps(f, p0);
					__m128 c1 = _mm_mul_ps(f, p1);
					__m128 o0 = _mm_load_ps(aBuffer + c);
					__m128 o1 = _mm_load_ps(aBuffer + c + aBufferSize);
					c0 = _mm_add_ps(c0, o0);
					c1 = _mm_add_ps(c1, o1);
					_mm_store_ps(aBuffer + c, c0);
					_mm_store_ps(aBuffer + c + aBufferSize, c1);
					p0 = _mm_add_ps(p0, pan0delta);
					p1 = _mm_add_ps(p1, pan1delta);
					c += 4;
				}
			}
			// If buffer size or samples to read are not divisible by 4, handle leftovers
			for (j = c; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				float s = aScratch[j];
				aBuffer[j + 0] += s * pan[0];
				aBuffer[j + aBufferSize] += s * pan[1];
			}
		}
#else // fallback
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				float s = aScratch[j];
				aBuffer[j + 0] += s * pan[0];
				aBuffer[j + aBufferSize] += s * pan[1];
			}
#endif
		break;
		}
		break;
	case 4:
		switch (aVoice->mChannels)
		{
		case 8: // 8->4, add a bit of center, sub?
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				pan[2] += pani[2];
				pan[3] += pani[3];
				float s1 = aScratch[j];
				float s2 = aScratch[aBufferSize + j];
				float s3 = aScratch[aBufferSize * 2 + j];
				float s4 = aScratch[aBufferSize * 3 + j];
				float s5 = aScratch[aBufferSize * 4 + j];
				float s6 = aScratch[aBufferSize * 5 + j];
				float s7 = aScratch[aBufferSize * 6 + j];
				float s8 = aScratch[aBufferSize * 7 + j];
				float c = (s3 + s4) * 0.7f;
				aBuffer[j + 0] += s1 * pan[0] + c;
				aBuffer[j + aBufferSize] += s2 * pan[1] + c;
				aBuffer[j + aBufferSize * 2] += 0.5f * (s5 + s7) * pan[2];
				aBuffer[j + aBufferSize * 3] += 0.5f * (s6 + s8) * pan[3];
			}
			break;
		case 6: // 6->4, add a bit of center, sub?
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				pan[2] += pani[2];
				pan[3] += pani[3];
				float s1 = aScratch[j];
				float s2 = aScratch[aBufferSize + j];
				float s3 = aScratch[aBufferSize * 2 + j];
				float s4 = aScratch[aBufferSize * 3 + j];
				float s5 = aScratch[aBufferSize * 4 + j];
				float s6 = aScratch[aBufferSize * 5 + j];
				float c = (s3 + s4) * 0.7f;
				aBuffer[j + 0] += s1 * pan[0] + c;
				aBuffer[j + aBufferSize] += s2 * pan[1] + c;
				aBuffer[j + aBufferSize * 2] += s5 * pan[2];
				aBuffer[j + aBufferSize * 3] += s6 * pan[3];
			}
			break;
		case 4: // 4->4
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				pan[2] += pani[2];
				pan[3] += pani[3];
				float s1 = aScratch[j];
				float s2 = aScratch[aBufferSize + j];
				float s3 = aScratch[aBufferSize * 2 + j];
				float s4 = aScratch[aBufferSize * 3 + j];
				aBuffer[j + 0] += s1 * pan[0];
				aBuffer[j + aBufferSize] += s2 * pan[1];
				aBuffer[j + aBufferSize * 2] += s3 * pan[2];
				aBuffer[j + aBufferSize * 3] += s4 * pan[3];
			}
			break;
		case 2: // 2->4
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				pan[2] += pani[2];
				pan[3] += pani[3];
				float s1 = aScratch[j];
				float s2 = aScratch[aBufferSize + j];
				aBuffer[j + 0] += s1 * pan[0];
				aBuffer[j + aBufferSize] += s2 * pan[1];
				aBuffer[j + aBufferSize * 2] += s1 * pan[2];
				aBuffer[j + aBufferSize * 3] += s2 * pan[3];
			}
			break;
		case 1: // 1->4
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				pan[2] += pani[2];
				pan[3] += pani[3];
				float s = aScratch[j];
				aBuffer[j + 0] += s * pan[0];
				aBuffer[j + aBufferSize] += s * pan[1];
				aBuffer[j + aBufferSize * 2] += s * pan[2];
				aBuffer[j + aBufferSize * 3] += s * pan[3];
			}
			break;
		}
		break;
	case 6:
		switch (aVoice->mChannels)
		{
		case 8: // 8->6
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				pan[2] += pani[2];
				pan[3] += pani[3];
				pan[4] += pani[4];
				pan[5] += pani[5];
				float s1 = aScratch[j];
				float s2 = aScratch[aBufferSize + j];
				float s3 = aScratch[aBufferSize * 2 + j];
				float s4 = aScratch[aBufferSize * 3 + j];
				float s5 = aScratch[aBufferSize * 4 + j];
				float s6 = aScratch[aBufferSize * 5 + j];
				float s7 = aScratch[aBufferSize * 6 + j];
				float s8 = aScratch[aBufferSize * 7 + j];
				aBuffer[j + 0] += s1 * pan[0];
				aBuffer[j + aBufferSize] += s2 * pan[1];
				aBuffer[j + aBufferSize * 2] += s3 * pan[2];
				aBuffer[j + aBufferSize * 3] += s4 * pan[3];
				aBuffer[j + aBufferSize * 4] += 0.5f * (s5 + s7) * pan[4];
				aBuffer[j + aBufferSize * 5] += 0.5f * (s6 + s8) * pan[5];
			}
			break;
		case 6: // 6->6
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				pan[2] += pani[2];
				pan[3] += pani[3];
				pan[4] += pani[4];
				pan[5] += pani[5];
				float s1 = aScratch[j];
				float s2 = aScratch[aBufferSize + j];
				float s3 = aScratch[aBufferSize * 2 + j];
				float s4 = aScratch[aBufferSize * 3 + j];
				float s5 = aScratch[aBufferSize * 4 + j];
				float s6 = aScratch[aBufferSize * 5 + j];
				aBuffer[j + 0] += s1 * pan[0];
				aBuffer[j + aBufferSize] += s2 * pan[1];
				aBuffer[j + aBufferSize * 2] += s3 * pan[2];
				aBuffer[j + aBufferSize * 3] += s4 * pan[3];
				aBuffer[j + aBufferSize * 4] += s5 * pan[4];
				aBuffer[j + aBufferSize * 5] += s6 * pan[5];
			}
			break;
		case 4: // 4->6
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				pan[2] += pani[2];
				pan[3] += pani[3];
				pan[4] += pani[4];
				pan[5] += pani[5];
				float s1 = aScratch[j];
				float s2 = aScratch[aBufferSize + j];
				float s3 = aScratch[aBufferSize * 2 + j];
				float s4 = aScratch[aBufferSize * 3 + j];
				aBuffer[j + 0] += s1 * pan[0];
				aBuffer[j + aBufferSize] += s2 * pan[1];
				aBuffer[j + aBufferSize * 2] += 0.5f * (s1 + s2) * pan[2];
				aBuffer[j + aBufferSize * 3] += 0.25f * (s1 + s2 + s3 + s4) * pan[3];
				aBuffer[j + aBufferSize * 4] += s3 * pan[4];
				aBuffer[j + aBufferSize * 5] += s4 * pan[5];
			}
			break;
		case 2: // 2->6
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				pan[2] += pani[2];
				pan[3] += pani[3];
				pan[4] += pani[4];
				pan[5] += pani[5];
				float s1 = aScratch[j];
				float s2 = aScratch[aBufferSize + j];
				aBuffer[j + 0] += s1 * pan[0];
				aBuffer[j + aBufferSize] += s2 * pan[1];
				aBuffer[j + aBufferSize * 2] += 0.5f * (s1 + s2) * pan[2];
				aBuffer[j + aBufferSize * 3] += 0.5f * (s1 + s2) * pan[3];
				aBuffer[j + aBufferSize * 4] += s1 * pan[4];
				aBuffer[j + aBufferSize * 5] += s2 * pan[5];
			}
			break;
		case 1: // 1->6
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				pan[2] += pani[2];
				pan[3] += pani[3];
				pan[4] += pani[4];
				pan[5] += pani[5];
				float s = aScratch[j];
				aBuffer[j + 0] += s * pan[0];
				aBuffer[j + aBufferSize] += s * pan[1];
				aBuffer[j + aBufferSize * 2] += s * pan[2];
				aBuffer[j + aBufferSize * 3] += s * pan[3];
				aBuffer[j + aBufferSize * 4] += s * pan[4];
				aBuffer[j + aBufferSize * 5] += s * pan[5];
			}
			break;
		}
		break;
	case 8:
		switch (aVoice->mChannels)
		{
		case 8: // 8->8
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				pan[2] += pani[2];
				pan[3] += pani[3];
				pan[4] += pani[4];
				pan[5] += pani[5];
				pan[6] += pani[6];
				pan[7] += pani[7];
				float s1 = aScratch[j];
				float s2 = aScratch[aBufferSize + j];
				float s3 = aScratch[aBufferSize * 2 + j];
				float s4 = aScratch[aBufferSize * 3 + j];
				float s5 = aScratch[aBufferSize * 4 + j];
				float s6 = aScratch[aBufferSize * 5 + j];
				float s7 = aScratch[aBufferSize * 6 + j];
				float s8 = aScratch[aBufferSize * 7 + j];
				aBuffer[j + 0] += s1 * pan[0];
				aBuffer[j + aBufferSize] += s2 * pan[1];
				aBuffer[j + aBufferSize * 2] += s3 * pan[2];
				aBuffer[j + aBufferSize * 3] += s4 * pan[3];
				aBuffer[j + aBufferSize * 4] += s5 * pan[4];
				aBuffer[j + aBufferSize * 5] += s6 * pan[5];
				aBuffer[j + aBufferSize * 6] += s7 * pan[6];
				aBuffer[j + aBufferSize * 7] += s8 * pan[7];
			}
			break;
		case 6: // 6->8
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				pan[2] += pani[2];
				pan[3] += pani[3];
				pan[4] += pani[4];
				pan[5] += pani[5];
				pan[6] += pani[6];
				pan[7] += pani[7];
				float s1 = aScratch[j];
				float s2 = aScratch[aBufferSize + j];
				float s3 = aScratch[aBufferSize * 2 + j];
				float s4 = aScratch[aBufferSize * 3 + j];
				float s5 = aScratch[aBufferSize * 4 + j];
				float s6 = aScratch[aBufferSize * 5 + j];
				aBuffer[j + 0] += s1 * pan[0];
				aBuffer[j + aBufferSize] += s2 * pan[1];
				aBuffer[j + aBufferSize * 2] += s3 * pan[2];
				aBuffer[j + aBufferSize * 3] += s4 * pan[3];
				aBuffer[j + aBufferSize * 4] += 0.5f * (s5 + s1) * pan[4];
				aBuffer[j + aBufferSize * 5] += 0.5f * (s6 + s2) * pan[5];
				aBuffer[j + aBufferSize * 6] += s5 * pan[6];
				aBuffer[j + aBufferSize * 7] += s6 * pan[7];
			}
			break;
		case 4: // 4->8
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				pan[2] += pani[2];
				pan[3] += pani[3];
				pan[4] += pani[4];
				pan[5] += pani[5];
				pan[6] += pani[6];
				pan[7] += pani[7];
				float s1 = aScratch[j];
				float s2 = aScratch[aBufferSize + j];
				float s3 = aScratch[aBufferSize * 2 + j];
				float s4 = aScratch[aBufferSize * 3 + j];
				aBuffer[j + 0] += s1 * pan[0];
				aBuffer[j + aBufferSize] += s2 * pan[1];
				aBuffer[j + aBufferSize * 2] += 0.5f * (s1 + s2) * pan[2];
				aBuffer[j + aBufferSize * 3] += 0.25f * (s1 + s2 + s3 + s4) * pan[3];
				aBuffer[j + aBufferSize * 4] += 0.5f * (s1 + s3) * pan[4];
				aBuffer[j + aBufferSize * 5] += 0.5f * (s2 + s4) * pan[5];
				aBuffer[j + aBufferSize * 6] += s3 * pan[4];
				aBuffer[j + aBufferSize * 7] += s4 * pan[5];
			}
			break;
		case 2: // 2->8
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				pan[2] += pani[2];
				pan[3] += pani[3];
				pan[4] += pani[4];
				pan[5] += pani[5];
				pan[6] += pani[6];
				pan[7] += pani[7];
				float s1 = aScratch[j];
				float s2 = aScratch[aBufferSize + j];
				aBuffer[j + 0] += s1 * pan[0];
				aBuffer[j + aBufferSize] += s2 * pan[1];
				aBuffer[j + aBufferSize * 2] += 0.5f * (s1 + s2) * pan[2];
				aBuffer[j + aBufferSize * 3] += 0.5f * (s1 + s2) * pan[3];
				aBuffer[j + aBufferSize * 4] += s1 * pan[4];
				aBuffer[j + aBufferSize * 5] += s2 * pan[5];
				aBuffer[j + aBufferSize * 6] += s1 * pan[6];
				aBuffer[j + aBufferSize * 7] += s2 * pan[7];
			}
			break;
		case 1: // 1->8
			for (j = 0; j < aSamplesToRead; j++)
			{
				pan[0] += pani[0];
				pan[1] += pani[1];
				pan[2] += pani[2];
				pan[3] += pani[3];
				pan[4] += pani[4];
				pan[5] += pani[5];
				pan[6] += pani[6];
				pan[7] += pani[7];
				float s = aScratch[j];
				aBuffer[j + 0] += s * pan[0];
				aBuffer[j + aBufferSize] += s * pan[1];
				aBuffer[j + aBufferSize * 2] += s * pan[2];
				aBuffer[j + aBufferSize * 3] += s * pan[3];
				aBuffer[j + aBufferSize * 4] += s * pan[4];
				aBuffer[j + aBufferSize * 5] += s * pan[5];
				aBuffer[j + aBufferSize * 6] += s * pan[6];
				aBuffer[j + aBufferSize * 7] += s * pan[7];
			}
			break;
		}
		break;
	}

	for (k = 0; k < aChannels; k++)
		aVoice->mCurrentChannelVolume[k] = pand[k];
}

} // namespace SoLoud
