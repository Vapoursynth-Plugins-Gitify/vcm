/* This file contains a  f1quiver function of FFTQuiver plugin for vapoursynth
// Row by row the image is transormed into frequency domain, frequency filtered and 
transformed back into row. In addition to a large number of Butterworth
filters, filter can be custom designed.

  This plugin needs any one of libfftw3f-3.dll,
  FFTW3 dll, fftw.dll to reside in path (may be windows\system32 folder)
  
Author V.C.Mohan. 
jun 2015, 14 sep 2020, 18 May 2021

  Copyright (C) <2014 - 2021>  <V.C.Mohan>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    A copy of the GNU General Public License is at
    see <http://www.gnu.org/licenses/>.
*/	
//---------------------------------------------------------------------------
/*
#include <stdlib.h>
#include "VapourSynth.h"
#include "VSHelper.h"

#include  <math.h>
#define _USE_MATH_DEFINES
#include <vector>
#include "windows.h"
#include "fftwlite.h"
#include "FQDomainHelper.h"
*/
//-----------------------------------------------------------------------------


typedef struct
{
		VSNodeRef *node;
		const VSVideoInfo *vi;
	
		bool  	test;		// is this a test?
		int  	row;		// test this start row
		int 	nrows;		//  number of rows to average to test
		float 	gamma;		// for spectrum scaling
		bool 	morph;		// is homomorphic process required	
		bool 	custom;		// filter is custom designed
		int		Fspec[64];
		int 	npoints;		// number of filters or for custo pairs specified

		float* FreqFilter;	//  freq response of filter 
		float* logLUT;
		fftwf_plan 	pf, pin;		// fftwf creates a plan of process. pointer to it

		int		wbest;	// wbest dimension for speed
		int 	nfft;	// number of points in this fft
		

#include "fftLateBindingClassParams.cpp"
		float* inBuf;
		fftwf_complex* outBuf;

	
} F1QuiverData;

	// filter designing functions
	//void f1BuildFilterCascade(float * FreqFilter, int * filterSpec, int nfft, int npoints);

	//void f1BuildCustomFilter(float * FrqFilt, int * specs, int nfft, int nval);
	
	//template <typename finc>
	//void f1DisplayHorizontalScale(int nyq, int wbest, int panelh, int wd, int pitch, finc * dp, finc max);
	template <typename finc>
	void f1DisplayPowerSpectrumAndFilter(float * powerspect, float * FreqFilter, float pscale, float pmax, float gamma, int panelh,int nfft,
						int wd, int pitch, finc * dp, finc max );
	template <typename finc>
	float f1GetSummedPowerspectrum(void** instanceData, float * in, fftwf_complex* out, float * powerspect,
		const finc * sp, const int pitch, const int wd);

	template <typename finc>
	void f1ProcessFullFrame(F1QuiverData* d, float * in, fftwf_complex* out,
		const finc * sp, finc *dp, const int pitch, const int wd , const int ht,finc min, finc max, float * LUT );


// This function is called immediately after vsapi->createFilter(). This is the only place where the video
// properties may be set. In this case we simply use the same as the input clip. You may pass an array
// of VSVideoInfo if the filter has more than one output, like rgb+alpha as two separate clips.
static void VS_CC f1quiverInit(VSMap *in, VSMap *out, void **instanceData, VSNode *node, 
								VSCore *core, const VSAPI *vsapi) 
{
    F1QuiverData *d = (F1QuiverData *) * instanceData;
    vsapi->setVideoInfo(d->vi, 1, node);

	int * facbuf = (int *) vs_aligned_malloc <int>(sizeof( int) *64, 32);	//maximum 64 factors, in this buf values filled are pairs of first is factor, second is dividend to be factored. At 
								// a value of 1 no more factors	
	//	wbest dimensions for speed. make sure starting with even number for width
	int wdEven = ((d->vi->width + 1) >> 1) << 1;
	d->wbest = getBestDim(wdEven, facbuf);
	
	vs_aligned_free(facbuf);	


#include "ConstructorCodeForLateBindingfft.cpp"

	if (!ok)
	{
		vsapi->setError(out, "vcm.f1quiver: could not load any of the dll or get required fnctions");
		if (d->hinstLib != NULL)
			FreeLibrary(d->hinstLib);

		vsapi->freeNode(d->node);
		return;
	}

	
	 // create fft plans. Requires buffers temporarily
	d->inBuf =  (float *)d->fftwf_malloc (sizeof(float) * d->wbest );

	d->outBuf = (fftwf_complex*) d->fftwf_malloc (sizeof(fftwf_complex) * (d->wbest/2+1));

	d->FreqFilter = (float*)d->fftwf_malloc(sizeof(float) * (d->wbest / 2 + 1));	// filter buffer

			// get fft sine cosine config buffers allocated by plan
	
	d->pf = d-> fftwf_plan_dft_r2c_1d( d->wbest, d->inBuf, d->outBuf, FFTW_MEASURE | FFTW_DESTROY_INPUT);

	d->pin = d->fftwf_plan_dft_c2r_1d( d->wbest, d->outBuf, d->inBuf, FFTW_MEASURE | FFTW_DESTROY_INPUT);

	
			// initialize freq response buffer with value of one
	for(int i = 0; i < d->wbest / 2 + 1; i ++)
		// all ones. will be multiplied by each filter with its constants (cascading)
		d->FreqFilter[i] = 1.0;			
	
	if (d->custom)
	{
		
		f1BuildCustomFilter( d->FreqFilter, d->Fspec, d->wbest,  d->npoints);
	}
	else
	{
		// Freq filters of all Input params are built and cascaded
		f1BuildFilterCascade(d->FreqFilter, d->Fspec, d->wbest, d->npoints);	
	}

	// normalize and scale by 1/ nfft 
	float fmax = 0;

	for( int i = 0; i <= d->wbest / 2; i ++)
	{
		fmax = fmax >= d->FreqFilter[i] ? fmax : d->FreqFilter[i];
	}

	float fscaler = 1.0f / (d->wbest * fmax);	// as fft scales up by nfft we compensate it here

	for(int i = 0; i <= d->wbest / 2; i ++)
		d->FreqFilter[i] *= fscaler;

	int nbits = d->vi->format->bitsPerSample;	

	if (d->morph && nbits >= 8 && nbits <= 12)
	{
		int nl = 1 << nbits;
		d->logLUT = (float*)d->fftwf_malloc(sizeof(float) * nl);

		for (int i = 0; i < nl; i++)
			d->logLUT[i] = (float)log((float)i);
	}
	else
		d->logLUT = NULL;
}

//---------------------------------------------------------------------------------------------------------------------------------

template <typename finc>
void f1DisplayPowerSpectrumAndFilter(float * powerspect, float * FreqFilter, float pscale, 
									float pmax, float gamma, int panelh,int nfft,
									int wd, int pitch, finc * dp, finc max )
{
	// pmax = 1.0;
	for(int i = 0; i < wd / 2; i++)
	{
					//power spectrum
		for (int h = 0; h <  (powerspect[i] * pscale) / pmax; h ++)
		{

			dp[(2 * panelh + 20 - h) * pitch +  i] = max; // normal scale
		}				
						// exp gamma scale
		for(int h = 0; h < pow(powerspect[i]  / pmax, gamma) * pscale; h ++)
		{

			dp[(panelh - h) * pitch + i] = max; // gamma scale
		}
					
					// designed filter
		dp[(panelh - (int)(pscale * FreqFilter[i] * nfft)) * pitch + i] = (max + max) / 3;
		// to make it bold 2 pixel wide
		dp[(panelh + 1 - (int)(pscale * FreqFilter[i] * nfft)) * pitch +  i] = (max + max) / 3;
		
	}
}
//--------------------------------------------------------------------------------------------------------------------
template <typename finc>
float f1GetSummedPowerspectrum(void** instanceData,  float * in, fftwf_complex* out,
							float * powerspect, const finc * sp, const int pitch, const int wd)
{
	F1QuiverData* d = (F1QuiverData*)*instanceData;
	// zero power spectrum buffer

	for(int i = 0; i < d->wbest / 2 + 1; i++)
	{				

		powerspect[i] = 0;
	}

	for( int r = d->row; r <  d->row + d->nrows; r ++)
	{
		// process each row and get average powerspectrum
	 
		 if (d->morph)

			getRowMorphInput(in, sp + r * pitch, d->wbest, wd, false, 1, d->logLUT);

		else

			getRowInput(in, sp + r * pitch, d->wbest, wd);

		// d->fftwf_execute_dft_r2c(d->pf, in, out);
		 d->fftwf_execute(d->pf);
						// sum powerspectra
		for(int i = 0; i < d->wbest / 2 + 1; i++)
		{				

			powerspect[i] += out[i][0] * out[i][0] + out[i][1] * out[i][1];
		}
	}
		// get max value
	float pmax = 0.0f;

	for(int i = 0; i < d->wbest/2 + 1; i ++)
	{
		if ( pmax < powerspect[i])

			pmax = powerspect[i];
	}


	return pmax;
}
//---------------------------------------------------------------------------------------------
template <typename finc>

void f1ProcessFullFrame(F1QuiverData* d,   float * in, fftwf_complex* out,
					const finc * sp, finc * dp, const int pitch, const int wd ,
					const int ht, finc min, finc max, float *LUT )
{
	
	bool center = false;
	int start = 1;

	for(int h = 0; h < ht ; h++)
	{
		if ( d->morph)
			getRowMorphInput(in, sp, d->wbest, wd, center, start, LUT);
		else
			getRowInput(in, sp, d->wbest, wd );

		//d->fftwf_execute_dft_r2c(d->pf, in, out);
		d->fftwf_execute(d->pf);

		F1ApplyFilter(out,d->FreqFilter, d->wbest /2 + 1);		

		//d->fftwf_execute_dft_c2r(d->pin,  out, in);
		d->fftwf_execute(d->pin);

		if( d->morph)

			getRowMorphOutput(in, dp, wd, min, max);
		else

			getRowOutput(in, dp, wd, min, max);

		sp += pitch;
		dp += pitch;

	}	
}

//---------------------------------------------------------------------------------------------------------------------
void scaleFloatInput( float * fp, float scale, int nval)
{
	for ( int i = 0; i < nval; i ++)
	{
		fp[i] *= scale;
	}
}
//---------------------------------------------------------------------------------------------------------
// This is the main function that gets called when a frame should be produced. It will, in most cases, get
// called several times to produce one frame. This state is being kept track of by the value of
// activationReason. The first call to produce a certain frame n is always arInitial. In this state
// you should request all the input frames you need. Always do it in ascending order to play nice with the
// upstream filters.
// Once all frames are ready, the filter will be called with arAllFramesReady. It is now time to
// do the actual processing.
static const VSFrameRef *VS_CC f1quiverGetFrame(int n, int activationReason, void **instanceData, void **frameData,
						VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi)
{
    F1QuiverData *d = (F1QuiverData *) * instanceData;

    if (activationReason == arInitial)
	{
        // Request the source frame on the first call
        vsapi->requestFrameFilter(n, d->node, frameCtx);
    } 
	else if (activationReason == arAllFramesReady)
	{
		//we are creating buffers here so that each thread has separate areas
		// so need to free them at end of GetFrame

		//-------------------------------------------------------------------

		const VSFrameRef* src = vsapi->getFrameFilter(n, d->node, frameCtx);

		const VSFormat* fi = d->vi->format;
		// process Green or Y component
		int plane = fi->colorFamily == cmRGB ? 1 : 0;
		int height = vsapi->getFrameHeight(src, plane);
		int width = vsapi->getFrameWidth(src, plane);
		VSFrameRef* dst = vsapi->copyFrame(src, core); //newVideoFrame(fi, width, height, src, core); 
		const uint8_t* srcp = vsapi->getReadPtr(src, plane);
		int src_stride = vsapi->getStride(src, plane);
		uint8_t* dstp = vsapi->getWritePtr(dst, plane);
		int nbytes = fi->bytesPerSample;
		int nbits = fi->bitsPerSample;
		int pitch = src_stride / nbytes;
		int ht = height;
		int wd = width; 		//
		int iwidth = d->wbest;
		int owidth = 2 + iwidth;
		// we can use those pointers in struct as not using threads
		//float* inBuf = (float*)d->fftwf_malloc (sizeof(float) * iwidth);

		//fftwf_complex* outBuf = (fftwf_complex*)d->fftwf_malloc(sizeof(fftwf_complex) * owidth);

		
		if (fi->sampleType == stInteger && nbits == 8)
		{
			uint8_t max = (1 << nbits) - 1, min = 0;

			f1ProcessFullFrame(d, d->inBuf, d->outBuf, srcp, dstp, pitch, wd, ht, min, max, d->logLUT);
		}

		else if (fi->sampleType == stInteger && nbits > 8)
		{
			uint16_t* dp = (uint16_t*)dstp;

			const uint16_t* sp = (const uint16_t*)srcp;

			uint16_t max = (1 << nbits) - 1, min = 0;

			f1ProcessFullFrame(d, d->inBuf, d->outBuf, sp, dp, pitch, wd, ht, min, max, d->logLUT);
		}

		else // float
		{
			float* dp = (float*)dstp;

			const float* sp = (const float*)srcp;

			float max = 1.0f, min = 0.0f;

			f1ProcessFullFrame(d, d->inBuf, d->outBuf, sp, dp, pitch, wd, ht, min, max, d->logLUT);
		}

			// Release the source frame
		vsapi->freeFrame(src);				
		return dst;
    }

    return 0;
}
//--------------------------------------------------------------------------------------------------
// test process
static const VSFrameRef* VS_CC f1qtestGetFrame(int n, int activationReason, void** instanceData, void** frameData,
	VSFrameContext* frameCtx, VSCore* core, const VSAPI* vsapi)
{
	F1QuiverData* d = (F1QuiverData*)*instanceData;

	if (activationReason == arInitial)
	{
		// Request the source frame on the first call
		vsapi->requestFrameFilter(n, d->node, frameCtx);
	}
	else if (activationReason == arAllFramesReady)
	{
		//we are creating buffers here so that each thread has separate areas
		// so need to free them at end of GetFrame

		//-------------------------------------------------------------------

		const VSFrameRef* src = vsapi->getFrameFilter(n, d->node, frameCtx);

		const VSFormat* fi = d->vi->format;
		// process Green or Y
		int plane = fi->colorFamily == cmRGB ? 1 : 0;
		int height = vsapi->getFrameHeight(src, plane);
		int width = vsapi->getFrameWidth(src, plane);
		VSFrameRef* dst = vsapi->copyFrame(src, core); //newVideoFrame(fi, width, height, src, core); 
		const uint8_t* srcp = vsapi->getReadPtr(src, plane);
		int src_stride = vsapi->getStride(src, plane);
		int nbits = fi->bitsPerSample;
		int nbytes = fi->bytesPerSample;
		uint8_t* dstp = vsapi->getWritePtr(dst, plane);
		int pitch = src_stride / nbytes;
		int ht = height;
		int wd = width; 		//
		int iwidth = d->wbest;
		int owidth = 2 + iwidth;

		//float* inBuf = (float*)d->fftwf_malloc(sizeof(float) * iwidth);

		//fftwf_complex* outBuf = (fftwf_complex*)d->fftwf_malloc(sizeof(fftwf_complex) * owidth);
		// in test processing we are not particular about time optimization
		float* powerspect = (float*)d->fftwf_malloc(sizeof(float) * iwidth);

		if (fi->sampleType == stInteger)
		{
			if (nbytes == 1)
			{
				// 8 bit samples
				const uint8_t* sp = (uint8_t*)srcp;
				uint8_t* dp = (uint8_t*)dstp;
				uint8_t gray = 1 << (nbits - 1);
				uint8_t max = (1 << nbits) - 1;
				uint8_t zero = 0;


				// get row values in to in, do  forward fft, sum square of complex numbers, find maximum
				float pmax = f1GetSummedPowerspectrum(instanceData, d->inBuf, d->outBuf, powerspect, sp, pitch, wd);

				// forward fft, apply filter, inverse fft plane, use right half of it in display

				f1ProcessFullFrame(d, d->inBuf, d->outBuf, sp, dp, pitch, wd, ht / 2, zero, max, d->logLUT);

				// zero out left half luma

				fillPlaneWithVal(dp, pitch, wd / 2, ht, zero);
				if (fi->colorFamily != cmRGB)
				{
					for (int p = 1; p < fi->numPlanes; p++)
					{
						int ht = vsapi->getFrameHeight(src, p);
						int wd = vsapi->getFrameWidth(src, p);
						const uint8_t* sp = vsapi->getReadPtr(src, p);
						int src_stride = vsapi->getStride(src, p);
						int pitch = src_stride / nbytes;
						uint8_t* dp = vsapi->getWritePtr(dst, p);
						// blacken left half background
						fillPlaneWithVal(dp, pitch, wd / 2, ht, gray);
					}
				}

				int panelh = (ht - 40) / 2;
				int pscale = panelh;

				if (pmax > 0.1f)	// pmax normally should be a large value dc value * nrows. zero only for a black clip
				{
					f1DisplayPowerSpectrumAndFilter(powerspect, d->FreqFilter, pscale,
						pmax, d->gamma, panelh, d->wbest,
						wd, pitch, dp, max);
				}

				f1DisplayHorizontalScale(NYQUIST, d->wbest, panelh, wd, pitch, dp, max);
			}

			else	// nb 9 to 16 bits per sample
			{
				const uint16_t* sp = (uint16_t*)srcp;
				uint16_t* dp = (uint16_t*)dstp;

				uint16_t gray = 1 << (nbits - 1);
				uint16_t max = (1 << nbits) - 1;
				uint16_t zero = 0;

				float pmax;
				// get row values in to in, do  forward fft, sum square of complex numbers, find maximum
				pmax = f1GetSummedPowerspectrum(instanceData, d->inBuf, d->outBuf, powerspect, sp, pitch, wd);

				// forward fft, apply filter, inverse fft top half frame, use right half of it in display

				f1ProcessFullFrame(d, d->inBuf, d->outBuf, sp, dp, pitch, wd, ht / 2, zero, max, d->logLUT);

				// zero out left half luma

				fillPlaneWithVal(dp, pitch, wd / 2, ht, zero);
				if (fi->colorFamily != cmRGB)
				{
					for (int p = 1; p < fi->numPlanes; p++)
					{
						int ht = vsapi->getFrameHeight(src, p);
						int wd = vsapi->getFrameWidth(src, p);
						const uint16_t* sp = (const uint16_t*)vsapi->getReadPtr(src, p);
						int src_stride = vsapi->getStride(src, p);
						int dst_stride = vsapi->getStride(src, p);
						int pitch = src_stride / fi->bytesPerSample;
						uint16_t* dp = (uint16_t*)vsapi->getWritePtr(dst, p);

						// blacken left half background
						fillPlaneWithVal(dp, pitch, wd / 2, ht, gray);
					}
				}

				int panelh = (ht - 40) / 2;
				int pscale = panelh;

				if (pmax > 0.1f)	// pmax normally should be a large value dc value * nrows. zero only for a black clip
				{
					f1DisplayPowerSpectrumAndFilter(powerspect, d->FreqFilter, pscale, pmax, d->gamma, panelh, d->wbest,
						wd, pitch, dp, max);
				}

				f1DisplayHorizontalScale(NYQUIST, d->wbest, panelh, wd, pitch, dp, max);
			}	// 16 bit

		}	// integer


		else	// if(fi->sampleType == stFloat)
		{
			const float* sp = (const float*)srcp;
			float* dp = (float*)dstp;

			float  gray = fi->colorFamily == cmRGB ? 0.5f : 0.0; // plane 1 & 2
			float  max = 1.0f;	// for plane 0
			float  zero = 0.0f; // plane 0

			float pmax;
			// get row values in to in, do  forward fft, sum square of complex numbers, find maximum
			pmax = f1GetSummedPowerspectrum(instanceData, d->inBuf, d->outBuf, powerspect, sp, pitch, wd);

			// forward fft, apply filter, inverse fft top half frame, use right half of it in display

			f1ProcessFullFrame(d, d->inBuf, d->outBuf, sp, dp, pitch, wd, ht / 2, zero, max, d->logLUT);

			// zero out left half luma

			fillPlaneWithVal(dp, pitch, wd / 2, ht, zero);
			if (fi->colorFamily != cmRGB)
			{
				for (int p = 1; p < fi->numPlanes; p++)
				{
					int ht = vsapi->getFrameHeight(src, p);
					int wd = vsapi->getFrameWidth(src, p);
					const float* sp = (const float*)vsapi->getReadPtr(src, p);
					int src_stride = vsapi->getStride(src, p);
					int dst_stride = vsapi->getStride(src, p);
					int pitch = src_stride / fi->bytesPerSample;
					float* dp = (float*)vsapi->getWritePtr(dst, p);

					// blacken left half background
					fillPlaneWithVal(dp, pitch, wd / 2, ht, gray);
				}
			}

			int panelh = (ht - 40) / 2;
			int pscale = panelh;

			if (pmax > 0.1f)	// pmax normally should be a large value dc value * nrows. zero only for a black clip
			{
				f1DisplayPowerSpectrumAndFilter(powerspect, d->FreqFilter, pscale, pmax, d->gamma, panelh, d->wbest,
					wd, pitch, dp, max);
			}

			f1DisplayHorizontalScale(NYQUIST, d->wbest, panelh, wd, pitch, dp, max);
		}

		if (fi->colorFamily == cmRGB)
		{
			// copy Green on to Blu and Red planes
			vs_bitblt(vsapi->getWritePtr(dst, 0), vsapi->getStride(dst, 0),
				vsapi->getWritePtr(dst, 1), vsapi->getStride(dst, 1),
				wd * nbytes, ht);
			vs_bitblt(vsapi->getWritePtr(dst, 2), vsapi->getStride(dst, 2),
				vsapi->getWritePtr(dst, 1), vsapi->getStride(dst, 1),
				wd * nbytes, ht);
		}

		d->fftwf_free(powerspect);
		vsapi->freeFrame(src);
		return dst;
	}
	return 0;
}


// Free all allocated data on filter destruction
static void VS_CC f1quiverFree(void *instanceData, VSCore *core, const VSAPI *vsapi)
{
    F1QuiverData *d = (F1QuiverData *)instanceData;
    vsapi->freeNode(d->node);
	d->fftwf_free (d->FreqFilter);
	d->fftwf_free(d->inBuf);
	d->fftwf_free(d->outBuf);
	if (d->logLUT != NULL)
		d->fftwf_free(d->logLUT);
	d->fftwf_destroy_plan(d->pf);
	d->fftwf_destroy_plan(d->pin);
	if (d->hinstLib != NULL)
		FreeLibrary(d->hinstLib);
    free(d);
}

// This function is responsible for validating arguments and creating a new filter
static void VS_CC f1quiverCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi)
{
    F1QuiverData d;
    F1QuiverData *data;
    int err;
	int temp;
    // Get a clip reference from the input arguments. This must be freed later.
    d.node = vsapi->propGetNode(in, "clip", 0, 0);
    d.vi = vsapi->getVideoInfo(d.node);

    // In this first version we only want to handle 8bit integer formats. Note that
    // vi->format can be 0 if the input clip can change format midstream.
    if (!isConstantFormat(d.vi) || d.vi->width == 0 || d.vi->height == 0 
		|| (d.vi->format->colorFamily != cmYUV 	&& d.vi->format->colorFamily != cmGray
			 &&  d.vi->format->colorFamily != cmRGB) )
	{
        vsapi->setError(out, "F1Quiver: only RGB, Yuv or Gray color constant formats and const frame dimensions input supported");
        vsapi->freeNode(d.node);
        return;
    }
	
	if (d.vi->format->sampleType == stFloat && d.vi->format->bitsPerSample == 16)
	{
		vsapi->setError(out, "F1Quiver: Half float formats not allowed ");
		vsapi->freeNode(d.node);
		return;
	}

    // If a property read fails for some reason (index out of bounds/wrong type)
    // then err will have flags set to indicate why and 0 will be returned. This
    // can be very useful to know when having optional arguments. Since we have
    // strict checking because of what we wrote in the argument string, the only
    // reason this could fail is when the value wasn't set by the user.
    // And when it's not set we want it to default to enabled.
    temp =  vsapi->propGetInt(in, "test", 0, &err);
    if (err)
        d.test = false;
	else
    // Let's pretend the only allowed values are 1 or 0...
    if (temp < 0 || temp > 1) 
	{
        vsapi->setError(out, "F1Quiver: test must be 0 or 1");
        vsapi->freeNode(d.node);
        return;
    } 
	else
		d.test = temp == 0? false : true;

	temp =   vsapi->propGetInt(in, "custom", 0, &err);
    if (err)
        d.custom = false;
	else
    // Let's pretend the only allowed values are 1 or 0...
    if (temp < 0 || temp > 1) 
	{
        vsapi->setError(out, "F1Quiver: custom must be 0 or 1");
        vsapi->freeNode(d.node);
        return;
    } 
	else
		d.custom = temp == 0? false : true;

	temp =   vsapi->propGetInt(in, "morph", 0, &err);
    if (err)
        d.morph = false;
	else
    // Let's pretend the only allowed values are 1 or 0...
    if (temp < 0 || temp > 1) 
	{
        vsapi->setError(out, "F1Quiver: morph must be 0 or 1");
        vsapi->freeNode(d.node);
        return;
    } 

	else
		d.morph = temp == 0? false : true;

	d.row = vsapi->propGetInt(in, "strow", 0, &err);
	if(err)
		d.row = 0;
	else
	{
		if(d.row < 0 || d.row >= d.vi->height)
		{
			
			vsapi->setError(out, "F1Quiver: strow must be in frame");
			vsapi->freeNode(d.node);
			return;
		}
    }

	d.nrows = vsapi->propGetInt(in, "nrows", 0, &err);

	if(err)
		d.nrows = d.vi->height / 2;
	else
	{
		if(d.nrows < 0 || d.row + d.nrows >= d.vi->height)
		{
			
			vsapi->setError(out, "F1Quiver: nrows must be one or more and strow +  nrows must be within frame height");
			vsapi->freeNode(d.node);
			return;
		}
    }

	if (d.test)
	{
		if (d.vi->height < 80)
		{
			vsapi->setError(out, "F1Quiver: frame height must be atleast 80 for test display. may add border  to make up");
			vsapi->freeNode(d.node);
			return;
		}

		d.gamma = vsapi->propGetFloat(in, "gamma", 0, &err);

		if(err)
			d.gamma = 0.05;
		else
		{
			if(d.gamma < 0.00001 || d.gamma > 1.0f)
			{
			
				vsapi->setError(out, "F1Quiver: gamma must be +ve and less than 1.0");
				vsapi->freeNode(d.node);
				return;
			}
		}
    }

	d.npoints = vsapi->propNumElements( in,"filter");
	if (d.npoints > 64 || d.npoints < 2 || (d.custom && ((d.npoints & 1) != 0 )) || ( ! d.custom && (d.npoints & 3) != 0 ) )
	{
		vsapi->setError(out, "F1Quiver: filter entries should not be more than 64, even number for custom and otherwise multiple of 4  ");
		vsapi->freeNode(d.node);
		return;
	}

	if ( d.custom)
	{
		temp = -1;
		for ( int i = 0; i < d.npoints; i += 2)
		{
			d.Fspec[i] = vsapi->propGetInt(in, "filter", i, 0);

			if(d.Fspec[i] <= temp || d.Fspec[i] > NYQUIST)
			{
				vsapi->setError(out, "F1Quiver: first value of filter pair must be in ascending order and less than NYQUIST ");
				vsapi->freeNode(d.node);
				return;
			}

			temp = d.Fspec[i];

			d.Fspec[i + 1] = vsapi->propGetInt(in, "filter", i + 1, 0);

			if(d.Fspec[i+ 1] <= 0 || d.Fspec[i+ 1] > 100)
			{
				vsapi->setError(out, "F1Quiver: second value of custom filter pair should be zero to 100 only ");
				vsapi->freeNode(d.node);
				return;
			}
		}
	}

	else		// not custom
	{
		for ( int i = 0; i < d.npoints; i += 4)
		{
			d.Fspec[i ] = vsapi->propGetInt(in, "filter", i , 0);

			if(d.Fspec[i] < 0 || d.Fspec[i] > 4)
			{
				vsapi->setError(out, "F1Quiver: first value of filter quartet should be 0 to 4 only ");
				vsapi->freeNode(d.node);
				return;
			}
			d.Fspec[i + 1] = vsapi->propGetInt(in, "filter", i + 1, 0);

			if(d.Fspec[i+ 1] <= 0 || d.Fspec[i+ 1] > NYQUIST)
			{
				vsapi->setError(out, "F1Quiver: Frequency the second value of filter pair should be zero to 100 only ");
				vsapi->freeNode(d.node);
				return;
			}

			d.Fspec[i + 2] = vsapi->propGetInt(in, "filter", i + 2, 0);

			if ( d.Fspec[0] == 3 && (d.Fspec[i+ 2] < d.Fspec[i + 1] || d.Fspec[i+ 2] > NYQUIST) )
			{
				vsapi->setError(out, "F1Quiver:  freq2 the third value of filter pair should not be less than freq or more than NYQUIST ");
				vsapi->freeNode(d.node);
				return;
			}
			else if (d.Fspec[0] == 3 && (d.Fspec[i + 2] <= 0 || d.Fspec[i + 2] > 100))
			{
				vsapi->setError(out, "F1Quiver: bandwidth being %age of freq the third value of filter pair should be 1 to 100 only ");
				vsapi->freeNode(d.node);
				return;
			}

			d.Fspec[i + 3] = vsapi->propGetInt(in, "filter", i + 3, 0);

			if(d.Fspec[i+ 3] <= 0 || d.Fspec[i+ 3] > 12)
			{
				vsapi->setError(out, "F1Quiver: degree the sharpness  value of filter pair should be 1 to 12 only ");
				vsapi->freeNode(d.node);
				return;
			}
		}
	}

	
    // I usually keep the filter data struct on the stack and don't allocate it
    // until all the input validation is done.
    data = (F1QuiverData *) malloc(sizeof(d));
    *data = d;

	if(d.test)
		vsapi->createFilter(in, out, "F1Quiver", f1quiverInit, f1qtestGetFrame, f1quiverFree, fmParallelRequests, 0, data, core);
	else
		vsapi->createFilter(in, out, "F1Quiver", f1quiverInit, f1quiverGetFrame, f1quiverFree, fmParallelRequests, 0, data, core);

}


