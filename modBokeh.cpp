/******************************************************************************
bokeh filter plugin for vapoursynth by V.C.Mohan
 Introduces Bokeh( Japanese term for blur) effect on an image
 which has slight blur on background.
 
 Author V.C.Mohan
 20 Dec 2020 mod on 31 Dec 2020

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
#include "statsAndOffsets.h"

//----------------------------
typedef struct
{
	VSNodeRef* nodeA;	// input image
	VSNodeRef* nodeB;	// heavily blurred  image
	const VSVideoInfo* vi;
	int grid;	// 2 to 64 
	float thresh;	// threshold for descrimination of sharpness in image
	int rgb[3];		// color components to be considered
	int yuv[3];		// color components to be considered
	
	int* circleLUT;
	int  count;
	float tsq;
	
} BokehData;
//------------------------------------------------------------------------------------------------------
template <typename finc>
float getAverage(const finc* sp, int* offsets, int noff);

template <typename finc>
bool isVarSharp(finc* sp, int* Offsets, int noff, float tsq);
//------------------------------------------------------------------------------------------------------

// This function is called immediately after vsapi->createFilter(). This is the only place where the video
// properties may be set. In this case we simply use the same as the input clip. You may pass an array
// of VSVideoInfo if the filter has more than one output, like rgb+alpha as two separate clips.
static void VS_CC bokehInit(VSMap* in, VSMap* out, void** instanceData, VSNode* node, VSCore* core, const VSAPI* vsapi)
{
	BokehData* d = (BokehData*)*instanceData;
	vsapi->setVideoInfo(d->vi, 1, node);
	const VSFrameRef* srcA = vsapi->getFrame(0, d->nodeA, NULL, 0);
	const VSFormat* fi = d->vi->format;
	if (fi->sampleType == stInteger)
		d->tsq = d->thresh * (1 << fi->bitsPerSample) * d->thresh * (1 << fi->bitsPerSample);
	else
		d->tsq = d->thresh * d->thresh;


	int nBytes = fi->bytesPerSample;
	int pitch = vsapi->getStride(srcA, 0) / nBytes;
	int gPoints = d->grid * d->grid * 4;
	// size of square encompassing circle of radius = grid
	d->circleLUT = (int*)vs_aligned_malloc <int>(sizeof(int) * gPoints, 32);

	d->count = makeCircularLUT(d->circleLUT, pitch, d->grid, 0);
	vsapi->freeFrame(srcA);

	if ( d->count >= gPoints)
	{
		vs_aligned_free(d->circleLUT);
		vsapi->setError(out, "bokeh:  count are in error");
		vsapi->freeNode(d->nodeA);
		vsapi->freeNode(d->nodeB);
		return;		
	}

}


template <typename finc>
bool isVarSharp(finc* sp, int* offsets, int noff, float tsq)
{
	float avg = getMeanValue(sp, offsets, noff);
	float var = getVariance(sp, offsets, noff, avg);
	if (var > tsq)
		return true;
	return false;
}


//----------------------------------------------------------------------------------------------------------

// This is the main function that gets called when a frame should be produced. It will, in most cases, get
// called several times to produce one frame. This state is being kept track of by the value of
// activationReason. The first call to produce a certain frame n is always arInitial. In this state
// you should request all the input frames you need. Always do it in ascending order to play nice with the
// upstream filters.
// Once all frames are ready, the filter will be called with arAllFramesReady. It is now time to
// do the actual processing.
static const VSFrameRef* VS_CC bokehGetFrame(int n, int activationReason, void** instanceData,
	void** frameData, VSFrameContext* frameCtx, VSCore* core, const VSAPI* vsapi)
{
	BokehData* d = (BokehData*)*instanceData;

	if (activationReason == arInitial)
	{
		// Request the source frame on the first call
		vsapi->requestFrameFilter(n, d->nodeA, frameCtx);
		vsapi->requestFrameFilter(n, d->nodeB, frameCtx);
	}
	else if (activationReason == arAllFramesReady)
	{
		const VSFrameRef* srcA = vsapi->getFrameFilter(n, d->nodeA, frameCtx);
		const VSFrameRef* srcB = vsapi->getFrameFilter(n, d->nodeB, frameCtx);
		VSFrameRef* dst = vsapi->copyFrame(srcB, core);
		// The reason we query this on a per frame basis is because we want our filter
		// to accept clips with varying dimensions. If we reject such content using d->vi
		// would be better.
		const VSFormat* fi = d->vi->format;

		int nBytes = fi->bytesPerSample;
		int pitch = vsapi->getStride(srcA, 0) / nBytes;
		int gPoints = d->grid * d->grid;
		// Originally created for each frame. Since vapoursynth guarentees stride 
		// to be same for frames having equal width taken to struct and init sections
		
		int* circleLUT = d->circleLUT; 
		float tsq = d->tsq;
		
		int count = d->count; // LUTCircleOffs(circleLUT, d->grid, pitch);
		// plane A skipped
		int np = fi->numPlanes > 3 ? 3 : fi->numPlanes;
		const unsigned char* ap[] = { NULL,NULL,NULL };
		unsigned char* dp[] = { NULL,NULL,NULL };
		// initialized for all planes to be processed.
		bool proc[] = { true, true, true }; 
		// mark which planes are to be skipped 
		for (int p = 0; p < 3; p++)
		{
			if (fi->colorFamily == cmRGB)
			{
				// need to convert  rgb to bgr
				if (d->rgb[p] == 0) proc[2 - p] = false;

			}
			else if (fi->colorFamily == cmYUV)
			{
				if (d->yuv[p] == 0) proc[p] = false;
			}
		}

		int ht = vsapi->getFrameHeight(srcA, 0);
		int wd = vsapi->getFrameWidth(srcA, 0);

		for (int p = 0; p < np; p++)
		{
			ap[p] = vsapi->getReadPtr(srcA, p);
			dp[p] = vsapi->getWritePtr(dst, p);
			// added these two for on 31-12-2020 Offset of h to starting value in for loop below
			ap[p] += d->grid * pitch * nBytes;
			dp[p] += d->grid * pitch * nBytes;
		}

		//  increment after each step along height in bytes
		int hInc = (d->grid / 2)* pitch * nBytes;
		

			// now process planes. in case of odd and even grid div by 2 giving problem. so - 1
		
		for (int h = d->grid; h < ht - d->grid - 1; h += d->grid / 2)
		{
			for (int w = d->grid ; w < wd - d->grid - 1; w += d->grid / 2)
			{
				bool sharp = false;

				if (fi->sampleType == stInteger)
				{

					if (nBytes == 1)
					{
						for (int p = 0; p < np; p++)
						{
							if (proc[p])
								sharp = isVarSharp(ap[p] + w, circleLUT, count, tsq);

							if (sharp)
								break;
						}

						if (sharp)
						{
							// copy surrounding circular area. if rectangular results look blocky 
							for (int p = 0; p < np; p++)
							{
								for (int i = 0; i < count; i++)
								{
									*(dp[p] + w + circleLUT[i]) 
										= *(ap[p] + w + circleLUT[i]);
								}
							}
						}

					}

					else if ( nBytes == 2) //   (nbits  9 to 16)
					{

						for (int p = 0; p < np; p++)
						{
							if( proc[p])

								sharp = isVarSharp(((uint16_t*)(ap[p])) + w, circleLUT, count, tsq);

							if (sharp)
								break;
						}

						if (sharp)
						{
							for (int p = 0; p < np; p++)
							{
								for (int i = 0; i < count; i++)
								{
									*(((uint16_t*)(dp[p])) + w + circleLUT[i])
										= *(((uint16_t*)(ap[p])) + w + circleLUT[i]);
								}
							}
						}
					}

				} // sample type integer


				else  if (nBytes == 4) // float
				{

					for (int p = 0; p < np; p++)
					{
						if (proc[p])

							sharp = isVarSharp((float*)(ap[p]) + w, circleLUT, count, tsq);

						if (sharp)
							break;
					}

					if (sharp)
					{
						for (int p = 0; p < np; p++)
						{
							for (int i = 0; i < count; i++)
							{
								*((float*)(dp[p]) + w + circleLUT[i])
									= *((float*)(ap[p]) + w + circleLUT[i]);
							}
						}
					}
				}
			}

			for (int p = 0; p < np; p++)
			{
				ap[p] += hInc;
				dp[p] += hInc;
			}
		}	// for int h=

		// vs_aligned_free(gridLUT);
		vsapi->freeFrame(srcA);
		vsapi->freeFrame(srcB);
		
		return dst;

	}	// all frames ready
	return 0;
}



// Free all allocated data on filter destruction
static void VS_CC bokehFree(void* instanceData, VSCore* core, const VSAPI* vsapi)
{
	BokehData* d = (BokehData*)instanceData;
	vsapi->freeNode(d->nodeA);
	vsapi->freeNode(d->nodeB);
	vs_aligned_free(d->circleLUT);
	free(d);
}

// This function is responsible for validating arguments and creating a new filter
static void VS_CC bokehCreate(const VSMap* in, VSMap* out, void* userData, VSCore* core, const VSAPI* vsapi) {
	BokehData d;
	BokehData* data;
	int err;

	// Get a clip reference from the input arguments. This must be freed later.
	d.nodeA = vsapi->propGetNode(in, "clip", 0, 0);
	d.vi = vsapi->getVideoInfo(d.nodeA);
	if (d.vi->format->colorFamily != cmRGB 
		&& d.vi->format->colorFamily != cmYUV
		&& d.vi->format->colorFamily != cmGray 		
		)
	{
		vsapi->setError(out, "bokeh: input must be rgb, yuv or y only.");
		vsapi->freeNode(d.nodeA);
		return;
	}
	if ( !isConstantFormat(d.vi)  || d.vi->format->subSamplingH != 0
		|| d.vi->format->subSamplingW != 0)
	{
		vsapi->setError(out, "bokeh: input must be constant and not subsampled");
		vsapi->freeNode(d.nodeA);
		return;
	}
	d.nodeB = vsapi->propGetNode(in, "clipb", 0, 0);
	const VSVideoInfo* bvi = vsapi->getVideoInfo(d.nodeB);

	if ( !isSameFormat ( d.vi, bvi) )
	{
		vsapi->setError(out, "bokeh: both clips must have identical formats");
		vsapi->freeNode(d.nodeA);
		vsapi->freeNode(d.nodeB);
		return;
	}

	d.grid = vsapi->propGetInt(in, "grid", 0, &err);

	if (err)
	{
		d.grid = 25;
	}
	else
	{
		if (d.grid < 2 || d.grid > 64)
		{
			vsapi->setError(out, "bokeh: value of grid must be between 2 and 64");
			vsapi->freeNode(d.nodeA);
			vsapi->freeNode(d.nodeB);
			return;
		}
	}
	d.thresh = vsapi->propGetFloat(in, "thresh", 0, &err);

	if (err)
	{
		d.thresh = 30.0f;
	}
	else
	{
		if (d.thresh < 0.0  || d.thresh > 1.0)
		{
			vsapi->setError(out, "bokeh: thresh must be 0 to 1 only");
			vsapi->freeNode(d.nodeA);
			vsapi->freeNode(d.nodeB);
			return;
		}
	}

	if (d.vi->format->colorFamily == cmRGB)
	{
		int count = vsapi->propNumElements(in, "rgb");
		if (count > 3)
		{
			vsapi->setError(out, "f2qBokeh: rgb array cannot have more than 3 entries.");
			vsapi->freeNode(d.nodeA);
			vsapi->freeNode(d.nodeB);
			return;
		}
		else if (count < 1)
		{
			d.rgb[0] = 1;
			d.rgb[1] = 1;
			d.rgb[2] = 1;
		}

		for (int p = 0; p < 3; p++)
		{
			d.rgb[p] = vsapi->propGetInt(in, "rgb", p, &err);

			if (err)
			{
				d.rgb[p] = d.rgb[p - 1];
			}

			else if (d.rgb[p] < 0 || d.rgb[p] > 1)
			{
				vsapi->setError(out, "f2qBokeh: rgb array can have values of 0 or 1 only.");
				vsapi->freeNode(d.nodeA);
				vsapi->freeNode(d.nodeB);
				return;
			}

		}
		if (d.rgb[0] == 0 && d.rgb[1] == 0 && d.rgb[2] == 0)
		{
			vsapi->setError(out, "f2qBokeh: rgb array all values should not be 0");
			vsapi->freeNode(d.nodeA);
			vsapi->freeNode(d.nodeB);
			return;
		}
	}
	else if (d.vi->format->colorFamily == cmYUV)
	{
		int count = vsapi->propNumElements(in, "yuv");
		if (count > 3)
		{
			vsapi->setError(out, "f2qBokeh: yuv array cannot have more than 3 entries.");
			vsapi->freeNode(d.nodeA);
			vsapi->freeNode(d.nodeB);
			return;
		}
		else if (count < 1)
		{
			d.yuv[0] = 1;
			d.yuv[1] = 1;
			d.yuv[2] = 1;
		}
		for (int p = 0; p < 3; p++)
		{
			d.yuv[p] = vsapi->propGetInt(in, "yuv", p, &err);

			if (err)
			{
				d.yuv[p] = d.yuv[p - 1];
			}

			else if (d.yuv[p] < 0 || d.yuv[p] > 1)
			{
				vsapi->setError(out, "f2qBokeh: yuv array can have values of 0 or 1 only.");
				vsapi->freeNode(d.nodeA);
				vsapi->freeNode(d.nodeB);
				return;
			}

		}
		if (d.yuv[0] == 0 && d.yuv[1] == 0 && d.yuv[2] == 0)
		{
			vsapi->setError(out, "f2qBokeh: yuv array all values should not be 0");
			vsapi->freeNode(d.nodeA);
			vsapi->freeNode(d.nodeB);
			return;
		}
	}

	// I usually keep the filter data struct on the stack and don't allocate it
	// until all the input validation is done.
	data = (BokehData*)malloc(sizeof(d));
	*data = d;


	// If your filter is really fast (such as a filter that only resorts frames) you should set the
	// nfNoCache flag to make the caching work smoother.
	vsapi->createFilter(in, out, "bokeh", bokehInit, bokehGetFrame, bokehFree, fmParallelRequests, 0, data, core);
}

//////////////////////////////////////////
/*

VS_EXTERNAL_API(void) VapourSynthPluginInit(VSConfigPlugin configFunc, VSRegisterFunction registerFunc, VSPlugin *plugin)
{
	configFunc("com.example.bokeh", "vmd", "VapourSynth bokeh ", VAPOURSYNTH_API_VERSION, 1, plugin);
	registerFunc("bokeh", "clip:clip;clipb:clip;grid:int:opt;thresh:float:opt;", bokehCreate, 0, plugin);
}
*/