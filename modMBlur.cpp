/******************************************************************************
MBlur filter plugin for vapoursynth by V.C.Mohan
 blurs in either circular (out of focus) or rectangular or linear (motion blur)
 borders are  unblurred. No input format restrictions
 Author V.C.Mohan	
 6 Nov 2020

  Copyright (C) <2020>  <V.C.Mohan>

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
typedef struct 
{
    VSNodeRef *node;
    const VSVideoInfo *vi;
	int type;	// 1 linear, 2 rectangular, 3.circular 
    int x, y;	
	int npoints;
} MBlurData;
//------------------------------------------------------------------------------------------------------
/*int makeLinearLUT(int* offsets, int pitch, int xcoord, int ycoord);
template <typename finc>
finc blurPoint ( int * offsets,  const finc * sp, int npts);
int makeRectLUT(int* offsets, int pitch,  int xcoord, int ycoord);
int makeCircularLUTforY(int* offsets, int pitch,  int radius);*/
int makeCircularLUTUV(int* offsets, int pitchUV,  int radius, int subW, int subH);


//--------------------------------------------------------------------------------
/*
template <typename finc>
finc blurPoint(int* offsets,  const finc* sp, int npts )
{
	float sum = 0;

	for (int i = 0; i < npts; i++)
		sum += (float)sp[offsets[i]];
	return  (finc)(sum / npts);
}
//------------------------------------------------------------------
*/
/*
int makeRectLUT(int* offsets, int pitch,  int xcoord, int ycoord)

{
	int count = 0;

	for (int h = -ycoord; h <= ycoord; h++)
	{

		for (int w = -xcoord; w <= xcoord; w++)
		{

			offsets[count] = h * pitch + w;

			count++;
		}
	}
	return count;
}
//--------------------------------------------------------------------------------------
*/
/*
int makeCircularLUTforY(int* offsets, int pitch,  int radius)
{
	int count = 0;

	for (int h = -radius; h <= radius; h++)
	{
		for (int w = -radius; w <= radius; w++)
		{
			if (h * h + w * w <= radius * radius)
			{				
				offsets[count] = h * pitch +  w;

				count++;
			}
		}
	}

	return count;		// we created buffer for a square shape
}
*/
//---------------------------------------------------------------------------------
int makeCircularLUTUV(int* offsets, int pitchUV,  int radius,  int subW, int subH)
{
	int count = 0;
	int andH = (1 << subH) - 1;
	int andW = (1 << subW) - 1;

	for (int h = -radius; h <= radius; h++)
	{
		if ((h & subH) == 0)
		{

			for (int w = -radius; w <= radius; w++)
			{
				if ((w & subW) == 0)
				{
					if (h * h + w * w <= radius * radius)
					{
						offsets[count] = (h >> subH) * pitchUV + (w >> subW);

						count++;
					}
				}
			}
		}
	}

	return count;		
}
//-------------------------------------------------------------------------------
/*
int makeLinearLUT(int* offsets, int pitch, int xcoord, int ycoord)
{
	int count = 0, npoints = 0;

	int absx = abs(xcoord), absy = abs(ycoord);

	if (absx < absy)
	{
		npoints = 2 * absy + 1;

		if (xcoord == 0)
		{
			// xcoord = 0. So avoid divide by zero
			for (int i = 0; i < npoints; i++)
			{
				offsets[i] = (i - npoints / 2) * pitch;
				count++;
			}
		}

		else
		{
			for (int i = 0; i < npoints; i++)
			{
				offsets[i] = (absy / ycoord) * (i - npoints / 2) * pitch
					+ absx / xcoord * (((i - npoints / 2) * absx
						+ (absy / 2)) / absy);
				count++;
			}
		}

	}
	else
	{
		npoints = 2 * absx + 1;
		if (ycoord == 0)	// avoid div by zero
		{

			for (int i = 0; i < npoints; i++)
			{
				offsets[i] = (i - npoints / 2);
				count++;
			}
		}

		else
		{
			for (int i = 0; i < npoints; i++)
			{
				offsets[i] = (absx / xcoord) * (i - npoints / 2)
					+ (absy / ycoord) * pitch * (((i - npoints / 2) * absy
						+ (absx / 2)) / absx);
				count++;
			}
		}
	}

	return count;
}

*/
//------------------------------------------------------------------------------------------------------

// This function is called immediately after vsapi->createFilter(). This is the only place where the video
// properties may be set. In this case we simply use the same as the input clip. You may pass an array
// of VSVideoInfo if the filter has more than one output, like rgb+alpha as two separate clips.
static void VS_CC mblurInit(VSMap *in, VSMap *out, void **instanceData, VSNode *node, VSCore *core, const VSAPI *vsapi)
{
    MBlurData *d = (MBlurData *) * instanceData;
    vsapi->setVideoInfo(d->vi, 1, node);
	int xcoord = d->x, ycoord = d->y, npoints;
	
	if (d->type == 1)	// linear
	{
		if (abs(xcoord) > abs(ycoord))		

			npoints = 2 * abs(xcoord) + 1;

		else

			npoints = 2 * abs(ycoord) + 1;
	}

	else if (d->type == 2)	// rectangular
	{

		npoints = (2 * xcoord + 1) * (2 * ycoord + 1);
		
	}

	else		// circular focus
	{
		xcoord = abs(xcoord);

		npoints = (2 * xcoord + 1) * (2 * xcoord + 1);
	}
	d->npoints = npoints;
	
}
//---------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------

// This is the main function that gets called when a frame should be produced. It will, in most cases, get
// called several times to produce one frame. This state is being kept track of by the value of
// activationReason. The first call to produce a certain frame n is always arInitial. In this state
// you should request all the input frames you need. Always do it in ascending order to play nice with the
// upstream filters.
// Once all frames are ready, the filter will be called with arAllFramesReady. It is now time to
// do the actual processing.
static const VSFrameRef *VS_CC mblurGetFrame(int n, int activationReason, void **instanceData, 
					void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) 
{
    MBlurData *d = (MBlurData *) * instanceData;

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
		int nbytes = fi->bytesPerSample;
        VSFrameRef *dst = vsapi->copyFrame(src, core);

		int * offset = (int *)vs_aligned_malloc<int>(sizeof(int) * d->npoints, 32);
		
		int subH = fi->subSamplingH;
		int subW = fi->subSamplingW;
		int count = d->npoints;
		int nplanes = fi->numPlanes > 3 ? 3 : fi->numPlanes;

		for (int plane = 0; plane < nplanes; plane++)
		{

			const uint8_t* srcp = vsapi->getReadPtr(src, plane);
			int src_stride = vsapi->getStride(src, plane);
			uint8_t* dstp = vsapi->getWritePtr(dst, plane);
			int dst_stride = vsapi->getStride(dst, plane);
			int ht = vsapi->getFrameHeight(src, plane);
			int wd = vsapi->getFrameWidth(src, plane);
			int pitch = dst_stride / nbytes;

			// create LUT of offsets to ease blurring
			if (d->type == 1)
			{
				if (plane == 0)
					count = makeLinearLUT(offset, pitch, d->x, d->y);
				else if (plane == 1)
				{
					if (fi->colorFamily == cmYUV && (fi->subSamplingH != 0 || fi->subSamplingW != 0))
						count = makeLinearLUT(offset, pitch, d->x >> subW, d->y >> subH);
				}
			}

			else if (d->type == 2)
			{
				if (plane == 0)					// rectangle
					count = makeRectGridLUT(offset, pitch, d->x, d->y);
				else if (plane == 1)
				{
					if (fi->colorFamily == cmYUV && (fi->subSamplingH != 0 || fi->subSamplingW != 0))
						count = makeRectGridLUT(offset, pitch, d->x >> subW, d->y >> subH);
				}
			}

			else if (d->type == 3)
			{								// create blur offsets due to bad focus
				if (plane == 0)
					count = makeCircularLUT(offset, pitch, d->x);
				else if (plane == 1)
				{
					if (fi->colorFamily == cmYUV && (fi->subSamplingH != 0 || fi->subSamplingW != 0))
						count = makeCircularLUTUV(offset, pitch, d->x, subW, subH);
				}
			}

			// now process planes

			int x = plane == 0 ? abs(d->x) : abs(d->x) >> subW;
			int y = plane == 0 ? abs(d->y) : abs(d->y) >> subH;

			dstp += y * src_stride;
			srcp += y * src_stride;

			for (int h = y; h < ht - y; h++)
			{
				for (int w = x; w < wd - x; w++)
				{

					if (fi->sampleType == stInteger)
					{
						int nb = fi->bitsPerSample;

						if (nb == 8)
						{
							uint8_t* dp = dstp;
							uint8_t* sp = (uint8_t*)srcp;
							dp[w] = (unsigned char) getMeanValue( sp + w, offset, count);

							//dp[w] = blurPoint(offset, sp + w, count);

						}
						else  // if (nb  9 to 16)
						{

							uint16_t* dp = ((uint16_t*)dstp);
							uint16_t* sp = ((uint16_t*)srcp);

							dp[w] = (uint16_t)getMeanValue(sp + w, offset, count);
						//	dp[w] = blurPoint(offset, sp + w, count);
						}
					} // sample type integer


					else if (nbytes == 4)
					{

						float* dp = ((float*)dstp);
						float* sp = ((float*)srcp);
						dp[w] = (float)getMeanValue(sp + w, offset, count);
						//dp[w] = blurPoint(offset, sp + w, count);

					}

				}

				srcp += src_stride;
				dstp += dst_stride;
			}			
		}
		vsapi->freeFrame(src);
		vs_aligned_free(offset);
		return dst;
			
	}	// all frames ready
	return 0;
}



// Free all allocated data on filter destruction
static void VS_CC mblurFree(void *instanceData, VSCore *core, const VSAPI *vsapi) 
{
    MBlurData *d = (MBlurData *)instanceData;
    vsapi->freeNode(d->node);
	
    free(d);
}

// This function is responsible for validating arguments and creating a new filter
static void VS_CC mblurCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi) {
    MBlurData d;
    MBlurData *data;
    int err;
	
    // Get a clip reference from the input arguments. This must be freed later.
    d.node = vsapi->propGetNode(in, "clip", 0, 0);
    d.vi = vsapi->getVideoInfo(d.node);
	if (d.vi->format->colorFamily != cmRGB && d.vi->format->colorFamily != cmYUV && d.vi->format->colorFamily != cmGray)
	{
		vsapi->setError(out, "mBlur: RGB, YUV and Gray color formats only for input allowed ");
		vsapi->freeNode(d.node);
		return;
	}
	if (d.vi->format->sampleType == stFloat && d.vi->format->bitsPerSample == 16)
	{
		vsapi->setError(out, "mBlur: Half float formats not allowed ");
		vsapi->freeNode(d.node);
		return;
	}
	d.type = vsapi->propGetInt(in, "type", 0, &err);
	if (err)
	{
		d.type = 1;
	}
	else
	{
		if (d.type < 1 || d.type > 3 )
		{
			vsapi->setError(out, "mBlur: Type value may be either 1 for linear, or  2 for rectangular or 3 for circular focal blur only");
			vsapi->freeNode(d.node);
			return;
		}
	}

	d.x = vsapi->propGetInt(in, "x", 0, &err);

	if (err)
	{
		d.x = 5;
	}
	else
	{
		if (abs (d.x) > 100)
		{
			vsapi->setError(out, "mBlur: x must be -100 to 100");
			vsapi->freeNode(d.node);
			return;
		}
	}
	d.y = vsapi->propGetInt(in, "y", 0, &err);

	if (err)
	{
		d.y = d.x;
	}
	else
	{
		if (abs (d.y) > 100 || (d.x == 0 && d.y == 0) )
		{
			vsapi->setError(out, "mBlur: y must be -100 to 100. Both x and y should not be zero");
			vsapi->freeNode(d.node);
			return;
		}
	}

    // I usually keep the filter data struct on the stack and don't allocate it
    // until all the input validation is done.
    data = ( MBlurData *) malloc(sizeof(d));
    *data = d;

   
    // If your filter is really fast (such as a filter that only resorts frames) you should set the
    // nfNoCache flag to make the caching work smoother.
    vsapi->createFilter(in, out, "mBlur", mblurInit, mblurGetFrame, mblurFree, fmParallel, 0, data, core);
}

//////////////////////////////////////////
/*

VS_EXTERNAL_API(void) VapourSynthPluginInit(VSConfigPlugin configFunc, VSRegisterFunction registerFunc, VSPlugin *plugin)
{
    configFunc("com.example.mblur", "vmd", "VapourSynth Motion Blur ", VAPOURSYNTH_API_VERSION, 1, plugin);
    registerFunc("MBlur", "clip:clip;type:int:opt;x:int:opt;y:int:opt;", mblurCreate, 0, plugin);
}
*/