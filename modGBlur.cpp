/******************************************************************************
gBlur filter plugin for vapoursynth by V.C.Mohan
 Calculates a Gaussian seperable kernel 3X3 or 5X5  up to 11x11 and applies to image
 The code for kernel computation is courtesy Mike Heath's published material
 GBlur is recommended to be applied on input image (Grey scale or color)
 to remove noise, so that cleaner segmentation, or masks or edges can be obtained.
 If input clip is constant format and integer type checks whether given ksize and sd combination
 is effective and gives a message if not. In case format is variable, on integer input if ksize
 is more than required reduces it to optimum size. For float input no checks or changes are made.
 Author V.C.Mohan	
 Oct, 2014, 20 Aug 2020

  Copyright (C) <2014- 2020>  <V.C.Mohan>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    A copy of the GNU General Public License is at
    see <http://www.gnu.org/licenses/>.
	 
*************************************************************************************************/  

/*
#include <stdlib.h>
#include "VapourSynth.h"
#include "VSHelper.h"
#include "math.h"
*/
//----------------------------
typedef struct {
    VSNodeRef *node;
    const VSVideoInfo *vi;
    int ksize;
	float sd;		// standard deviation
	
	float *kern;
	
} GBlurData;
//------------------------------------------------------------------------------------------------------

template <typename finc>
void blurPlane ( finc * dptr, float * kern, int ksize, int pitch, int ht, int wd, int * koffset);
//------------------------------------------------------------------------------------------------------

// This function is called immediately after vsapi->createFilter(). This is the only place where the video
// properties may be set. In this case we simply use the same as the input clip. You may pass an array
// of VSVideoInfo if the filter has more than one output, like rgb+alpha as two separate clips.
static void VS_CC gblurInit(VSMap *in, VSMap *out, void **instanceData, VSNode *node, VSCore *core, const VSAPI *vsapi) {
    GBlurData *d = (GBlurData *) * instanceData;
    vsapi->setVideoInfo(d->vi, 1, node);

	d->kern = (float*) vs_aligned_malloc <float>( sizeof(float) * d->ksize, 32);
		// this kern will be used by float formats

		float sum = 0;		

		for( int i = 0; i < d->ksize; i++)
		{
			d->kern[i] = pow(2.71828, -(0.5*(i - d-> ksize / 2) * (i - d->ksize / 2) ) / ( d->sd * d->sd) ) / (d->sd * sqrt(6.2831853));
			
			sum += d->kern[i];
		}

		for(int i = 0; i < d->ksize; i ++)

			d->kern[i] /=  sum ;	// so that in blur we only do integer multiplications and a final div by bit shift

		
}
//---------------------------------------------------------------------------------------------------------

// convolve with kernel first alonng width and then along height of frame
template <typename finc>
void blurPlane ( finc * dptr, float * kern, int ksize, int pitch, int ht, int wd, int * koffset)
{
	finc * dp = dptr;

	for ( int h = 0; h < ht; h ++)
	{
		for( int w = 0; w < wd - ksize ; w ++)
		{
			float sum = 0.0;

			for ( int k = 0 ; k < ksize; k ++)
			{
				sum += dp[w + k] * kern[k];
			}

			dp[w + ksize / 2] = sum; 
		}

		dp += pitch;
	}

			// blur along height direction
	int kcenterOffset = ksize / 2 * pitch;
						
	for ( int w = 0; w < wd; w ++)
	{
		finc * dp = dptr + w;

		for( int h = 0; h < ht - ksize ; h ++)
		{
			float sum = 0.0;							

			for ( int k = 0 ; k < ksize; k ++)
			{
				sum += dp[koffset[k] ] * kern[k];
			}

			dp[kcenterOffset] = sum ; /// d->ksize;

			dp += pitch;
		}
							
	}
}

//----------------------------------------------------------------------------------------------------------

// This is the main function that gets called when a frame should be produced. It will, in most cases, get
// called several times to produce one frame. This state is being kept track of by the value of
// activationReason. The first call to produce a certain frame n is always arInitial. In this state
// you should request all the input frames you need. Always do it in ascending order to play nice with the
// upstream filters.
// Once all frames are ready, the filter will be called with arAllFramesReady. It is now time to
// do the actual processing.
static const VSFrameRef *VS_CC gblurGetFrame(int n, int activationReason, void **instanceData, 
					void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) 
{
    GBlurData *d = (GBlurData *) * instanceData;

    if (activationReason == arInitial)
	{
        // Request the source frame on the first call
        vsapi->requestFrameFilter(n, d->node, frameCtx);	
    } 
	else if (activationReason == arAllFramesReady) 
	{
        const VSFrameRef *src = vsapi->getFrameFilter(n, d->node, frameCtx);
        // The reason we query this on a per frame basis is because we want our filter
        // to accept clips with varying dimensions. If we reject such content using d->vi
        // would be better.
        const VSFormat *fi = d->vi->format;
        int height = vsapi->getFrameHeight(src, 0);
        int width = vsapi->getFrameWidth(src, 0);
		int kb = fi->bytesPerSample;
		int ksize = d->ksize;

		if (fi->sampleType == stInteger)
		{
			int nb = fi->bitsPerSample;

			int maxval =  (1 << nb) - 1;
			
			while (ksize > 1)
			{
				// check if atleast it makes a difference of one to the addition. If not reduce ksize
				if( maxval * (pow(2.71828, -(0.5*( ksize / 2) * ( ksize / 2) ) / ( d->sd * d->sd) ) / (d->sd * sqrt(6.2831853)) ) >= 1)

					break;
				ksize -= 2;
			}

			if (ksize < 3)
			{
				// nothing to be done
				return src;
			}
		}

        // When creating a new frame for output it is VERY EXTREMELY SUPER IMPORTANT to
        // supply the "dominant" source frame to copy properties from. Frame props
        // are an essential part of the filter chain and you should NEVER break it.
        VSFrameRef *dst = vsapi->copyFrame(src, core);

        // It's processing loop time!
        // Loop over all the planes
        
        for (int plane = 0; plane < fi->numPlanes; plane++)
		{
			if( plane == 0 || ( fi->colorFamily == cmYUV && fi->subSamplingH == 0 && fi->subSamplingW == 0)  || fi->colorFamily == cmRGB)			
			{
				const uint8_t *srcp = vsapi->getReadPtr(src, plane);
				int src_stride = vsapi->getStride(src, plane);
				uint8_t *dstp = vsapi->getWritePtr(dst, plane);
				int dst_stride = vsapi->getStride(dst, plane);
				int ht = vsapi->getFrameHeight(src, plane);            
				int wd = vsapi->getFrameWidth(src, plane);
				int pitch = dst_stride / kb;
				int koffset[11];

				for ( int i = 0; i < d->ksize; i ++)
				{

					koffset[i] = i * pitch;
				}

			
				if (fi->sampleType == stInteger)
				{
					int nb = fi->bitsPerSample;	

					float kern[11];
					float sum = 0;

					for ( int i = 0; i < ksize; i ++)
					{				

						kern[i] = pow(2.71828, -(0.5*(i - ksize / 2) * (i - ksize / 2) ) / ( d->sd * d->sd) ) / (d->sd * sqrt(6.2831853));
			
						sum += kern[i];
					}

					for(int i = 0; i < ksize; i ++)
					{

						kern[i] /=  sum ;	// so that in blur we only do integer multiplications and a final div by bit shift
					}

					if ( nb == 8)
					{
						uint8_t * dp = dstp;
						
						blurPlane (dp, kern, ksize, pitch,  ht,  wd,  koffset);
					
					}

					else  // if (nb > 8)
					{
						
						uint16_t * dp = (uint16_t *)dstp;
						
						blurPlane (dp, kern, ksize, pitch,  ht,  wd,  koffset);
					
					}
				}	// sample type integer

				else // float
				{
					float * dp = (float *)dstp;
						
					blurPlane (dp, d->kern, d->ksize, pitch,  ht,  wd,  koffset);
					
				}

			}	// if rgb, plane,-----
			
		}	// plane
		vsapi->freeFrame(src);
		return dst;
			
	}	// all frames ready
	return 0;
}



// Free all allocated data on filter destruction
static void VS_CC gblurFree(void *instanceData, VSCore *core, const VSAPI *vsapi) 
{
    GBlurData *d = (GBlurData *)instanceData;
    vsapi->freeNode(d->node);
	vs_aligned_free(d->kern);
    free(d);
}

// This function is responsible for validating arguments and creating a new filter
static void VS_CC gblurCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi) {
    GBlurData d;
    GBlurData *data;
    int err;
	
    // Get a clip reference from the input arguments. This must be freed later.
    d.node = vsapi->propGetNode(in, "clip", 0, 0);
    d.vi = vsapi->getVideoInfo(d.node);
	if (d.vi->format->colorFamily == cmCompat)
	{
		vsapi->setError(out, "gBlur: compat format input is not supported");
		vsapi->freeNode(d.node);
		return;
	}
	d.ksize = vsapi->propGetInt(in, "ksize", 0, &err);
	if (err)
	{
		d.ksize = 5;
	}
	else
	{
		if (d.ksize < 3 || d.ksize > 11 || (d.ksize % 2) == 0)
		{
			vsapi->setError(out, "gBlur: ksize need to be an odd number between 3 and 11");
			vsapi->freeNode(d.node);
			return;
		}
	}

	d.sd = vsapi->propGetFloat(in, "sd", 0, &err);

	if (err)
	{
		d.sd = 1.5;
	}
	else
	{
		if (d.sd < 0.01)
		{
			vsapi->setError(out, "gBlur: sd must have a value above 0.01");
			vsapi->freeNode(d.node);
			return;
		}
	}

	if (isConstantFormat(d.vi) && d.vi->format->sampleType == stInteger)
	{
		int maxval =  (1 << d.vi->format->bitsPerSample) - 1;
		int temp = d.ksize;
		while (temp > 1)
		{
			if( maxval * (pow(2.71828, -(0.5*( temp / 2) * ( temp / 2) ) / ( d.sd * d.sd) ) / (d.sd * sqrt(6.2831853)) ) >= 1)

				break;
			temp -= 2;
		}

		if (temp < d.ksize)
		{
			vsapi->setError(out, "gBlur: either decrease ksize or increase sd to be effective");
			vsapi->freeNode(d.node);
			return;
		}
	}


    // I usually keep the filter data struct on the stack and don't allocate it
    // until all the input validation is done.
    data = ( GBlurData *) malloc(sizeof(d));
    *data = d;

   
    // If your filter is really fast (such as a filter that only resorts frames) you should set the
    // nfNoCache flag to make the caching work smoother.
    vsapi->createFilter(in, out, "gBlur", gblurInit, gblurGetFrame, gblurFree, fmParallel, 0, data, core);
}

//////////////////////////////////////////

/*
VS_EXTERNAL_API(void) VapourSynthPluginInit(VSConfigPlugin configFunc, VSRegisterFunction registerFunc, VSPlugin *plugin)
{
    configFunc("com.example.gblur", "vmod", "VapourSynth GBlur Example", VAPOURSYNTH_API_VERSION, 1, plugin);
    registerFunc("GBlur", "clip:clip;ksize:int:opt;sd:float:opt;", gblurCreate, 0, plugin);
}
*/