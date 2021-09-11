/* --------------------------------------------------------------------------------------------------------------
variance vapoursynth plugin
Computes the global variance  in the window of frame or clip,  in each grid of specified size
  local mean and variance are used for this adaptive noise filtering .If variance is to be 
  computed for each frame, then the window can be linearly moved in the proceessing  range of frames.

  Instead of using  the computed global variance, its value can be specified and also 
  linearly varied in the proceessing frame range
  
 input clip must be of constant format
 last modified on 20 Aug 2020
  Copyright (C) <2006, 2020>  <V.C.Mohan>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    A copy of the GNU General Public License is at
    http://www.gnu.org/licenses/.
	
  Author V.C.Mohan

  Aug 2020

//---------------------------------------------------------------------------------------------------------------
*/

/*
#include <stdlib.h>
#include "VapourSynth.h"
#include "VSHelper.h"
*/
typedef struct {
		VSNodeRef *node;
		const VSVideoInfo *vi;	
		int xgrid, ygrid;		// size of grid x and to be used. Must be submultiples of frame width and height
		bool UV;
		float gvarsq[3];
		int lx, wd, ty, ht; // window left x, width, top y and height coordinates
		int fn;						//frame number from which global values are to be computed.For frame( clip =false) these can be coordinates in start frame 		
		} VarianceData;

//------------------------------------------------------------------------
template <typename finc >
float	getVarSq( const finc * sp, int pitch, int lx, int wd, 
					int ty, int ht, float lmean );
template <typename finc>
float getMeanVal( const finc * sp, int pitch, int lx, int wd, int ty, int ht );
template <typename finc>
void varianceGrid(finc * dp, const finc * sp, int pitch, float gvarsq,
					int ht, int wd, int xgrid, int ygrid);

//----------------------------------------------------------------------------

template <typename finc>
float getMeanVal( const finc * sp, int pitch, int lx, int wd, int ty, int ht )
{
	sp += ty * pitch;

	float sum = 0;
	for(int h = 0; h <  ht; h ++)
	{
		for ( int w = lx; w < lx + wd; w ++)
		{
			sum += sp[w];
		}

		sp += pitch;
	}
	return sum / ( wd * ht);
}
//--------------------------------------------------------------------------------------
	template <typename finc >
float	getVarSq( const finc * sp, int pitch, int lx, int wd, 
					int ty, int ht, float lmean )
{
		float sum = 0;
		sp += ty * pitch;

		for(int h = 0; h <  ht; h ++)
		{

			for(int w = lx; w < lx + wd; w ++)
			{

				sum += ( sp[w] - lmean ) * ( sp[w] - lmean );
			}

			sp += pitch;
		}

		return sum/( wd * ht);
}



// This function is called immediately after vsapi->createFilter(). This is the only place where the video
// properties may be set. In this case we simply use the same as the input clip. You may pass an array
// of VSVideoInfo if the filter has more than one output, like rgb+alpha as two separate clips.
static void VS_CC varianceInit(VSMap *in, VSMap *out, void **instanceData, VSNode *node, VSCore *core, const VSAPI *vsapi) 
{
    VarianceData *d = (VarianceData *) * instanceData;
    vsapi->setVideoInfo(d->vi, 1, node);

	
		// get variance for use in entire clip
		const VSFrameRef *src = vsapi->getFrame(d->fn, d->node,0,0);

		const VSFormat *fi = d->vi->format;
        int height = vsapi->getFrameHeight(src, 0);
        int width = vsapi->getFrameWidth(src, 0);
		int nb = fi->bitsPerSample;
		int kb = fi->bytesPerSample;
		int subW[] = {0, fi->subSamplingW, fi->subSamplingW};
		int subH[] = {0, fi->subSamplingH, fi->subSamplingH};
		if (fi->colorFamily == cmRGB)
		{
			
			d->ty = height  - d->ty - d->ht;
			
		}

		 for (int plane = 0; plane < fi->numPlanes; plane++)
		 {
            const uint8_t *srcp = vsapi->getReadPtr(src, plane);
            int src_stride = vsapi->getStride(src, plane);
            int ht = vsapi->getFrameHeight(src, plane);            
            int wd = vsapi->getFrameWidth(src, plane);
			int pitch = src_stride / kb;

			if(plane == 0 || (fi->colorFamily == cmYUV && d->UV) || fi->colorFamily == cmRGB )
			{

				// getMeanVal( const finc * sp, int pitch, int lx, int rx, int ty, int by );

				if (fi->sampleType == stInteger)
				{
					if( nb == 8)
					{
						float lmean = getMeanVal(srcp, pitch, d->lx >> subW[plane], d->wd >> subW[plane],
												d->ty >> subH[plane], d->ht >> subH[plane]);

						d->gvarsq[plane] = getVarSq(srcp, pitch, d->lx >> subW[plane], d->wd >> subW[plane],
												d->ty >> subH[plane], d->ht >> subH[plane], lmean);
					}

					else	// 10 or 12  16bit samples
					{
						uint16_t * sp = (uint16_t *) srcp;

						float lmean = getMeanVal(sp, pitch, d->lx >> subW[plane], d->wd >> subW[plane],
												d->ty >> subH[plane], d->ht >> subH[plane]);

						d->gvarsq[plane] = getVarSq(sp, pitch, d->lx >> subW[plane], d->wd >> subW[plane],
												d->ty >> subH[plane], d->ht >> subH[plane], lmean);
					}
				}
				else	// float
				{
					float * sp = (float *) srcp;

					float lmean = getMeanVal(sp, pitch, d->lx >> subW[plane], d->wd >> subW[plane],
												d->ty >> subH[plane], d->ht >> subH[plane]);

					d->gvarsq[plane] = getVarSq(sp, pitch, d->lx >> subW[plane], d->wd >> subW[plane],
												d->ty >> subH[plane], d->ht >> subH[plane], lmean);
				}
			}	// if plane

		}	// for plane

	
	vsapi->freeFrame(src);
}
//-------------------------------------------------------------------------------
template <typename finc>
void varianceGrid(finc * dp, const finc * sp, int pitch, float gvarsq,
					int ht, int wd, int xgrid, int ygrid)
{
	// grid center starting position
	int offset = (ygrid / 2) * pitch + xgrid / 2;
	

	for (int h = 0; h < ht - ygrid  ; h ++)
	{
		for ( int w = 0; w < wd - xgrid ; w ++)
		{
			float gridmean = getMeanVal(sp, pitch, w,  xgrid , 0,  ygrid);
			float gridvar =  getVarSq(sp, pitch, w,  xgrid, 0,  ygrid, gridmean);

			if(gridvar > gvarsq)
						// local is high. So may be an edge. Ensure output is very near to this value
				dp [offset + w] = sp [offset + w] - (gvarsq * (sp[ offset + w] - gridmean))/gridvar;
			else
				dp [offset + w] = gridmean;

		}

	
	
		sp += pitch;
		dp += pitch;
	}
}
//----------------------------------------------------------------------------------------------
// This is the main function that gets called when a frame should be produced. It will, in most cases, get
// called several times to produce one frame. This state is being kept track of by the value of
// activationReason. The first call to produce a certain frame n is always arInitial. In this state
// you should request all the input frames you need. Always do it in ascending order to play nice with the
// upstream filters.
// Once all frames are ready, the filter will be called with arAllFramesReady. It is now time to
// do the actual processing.
static const VSFrameRef *VS_CC varianceGetFrame(int n, int activationReason, void **instanceData, void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi)
{
    VarianceData *d = (VarianceData *) * instanceData;

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
		VSFrameRef *dst = vsapi->copyFrame( src, core);		
		int nb = fi->bitsPerSample;
		int kb = fi->bytesPerSample;
        // When creating a new frame for output it is VERY EXTREMELY SUPER IMPORTANT to
        // supply the "dominant" source frame to copy properties from. Frame props
        // are an essential part of the filter chain and you should NEVER break it.
			 

        // It's processing loop time!
        // Loop over all the planes
        
        for (int p = 0; p < fi->numPlanes; p++)
		{
            const uint8_t *srcp = vsapi->getReadPtr(src, p);
            int src_stride = vsapi->getStride(src, p);
            uint8_t *dstp = vsapi->getWritePtr(dst, p);
            int dst_stride = vsapi->getStride(dst, p); // note that if a frame has the same dimensions and format, the stride is guaranteed to be the same. int dst_stride = src_stride would be fine too in this filter.
            // Since planes may be subsampled you have to query the height of them individually
            int ht = vsapi->getFrameHeight(src, p);            
            int wd = vsapi->getFrameWidth(src, p);
			int pitch = src_stride / kb;	// pitch in pixels

			if( fi->sampleType == stInteger)
			{

				 if( d->gvarsq[p] < 1 || (fi->colorFamily == cmYUV && ! d->UV) )
				// gvarsq = 0 means no noise. UV process not opted
					continue;
				if ( nb == 8)
				{
					uint8_t * dp = (uint8_t *) dstp;
					const uint8_t * sp = (uint8_t *) srcp;

					varianceGrid(dp, sp, pitch, d->gvarsq[p], ht, wd, d->xgrid, d->ygrid);
				}
				else	// 10 or 12 bit
				{
					uint16_t * dp = (uint16_t *) dstp;
					const uint16_t * sp = (uint16_t *) srcp;

					varianceGrid(dp, sp, pitch, d->gvarsq[p], ht, wd, d->xgrid, d->ygrid);
				}
			}

			else	// float
			{
				float * dp = (float *) dstp;
				const float * sp = (float *) srcp;

				varianceGrid(dp, sp, pitch, d->gvarsq[p], ht, wd, d->xgrid, d->ygrid);
			}
		}

		vsapi->freeFrame(src);

		return dst;
	}

	return 0;
}


// Free all allocated data on filter destruction
static void VS_CC varianceFree(void *instanceData, VSCore *core, const VSAPI *vsapi) {
    VarianceData *d = (VarianceData *)instanceData;
    vsapi->freeNode(d->node);
    free(d);
}

// This function is responsible for validating arguments and creating a new filter
static void VS_CC varianceCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi)
{
    VarianceData d;
    VarianceData *data;
    int err;
	int temp;

    // Get a clip reference from the input arguments. This must be freed later.
    d.node = vsapi->propGetNode(in, "clip", 0, 0);
    d.vi = vsapi->getVideoInfo(d.node);

    
    // vi->format can be 0 if the input clip can change format midstream.
    if (!isConstantFormat(d.vi) || d.vi->format->colorFamily == cmCompat ) 
	{
        vsapi->setError(out, "variance: input clip of only constant and other than Compat format allowed");
        vsapi->freeNode(d.node);
        return;
    }
	d.lx = vsapi->propGetInt(in, "lx", 0, 0);
	d.ty = vsapi->propGetInt(in, "ty", 0, 0);
	d.wd = vsapi->propGetInt(in, "wd", 0, 0);
	d.ht = vsapi->propGetInt(in, "ht", 0, 0);

	if(d.lx < 0 || d.lx >= d.vi->width)
	{
        vsapi->setError(out, "variance: lx is out of frame");
        vsapi->freeNode(d.node);
        return;
    }
	if(d.wd <= 0 || d.lx + d.wd >= d.vi->width)
	{
        vsapi->setError(out, "variance: invalid wd. wd must be +ve number and lx + wd must be in frame");
        vsapi->freeNode(d.node);
        return;
    }
	if(d.ty < 0 || d.ty >= d.vi->height)
	{
        vsapi->setError(out, "variance: ty is out of frame");
        vsapi->freeNode(d.node);
        return;
    }
	if(d.ht <= 0 || d.ty + d.ht >= d.vi->height)
	{
        vsapi->setError(out, "variance: invalid ht. ht must be +ve number and ty + ht must be in frame");
        vsapi->freeNode(d.node);
        return;
    }
	d.fn = vsapi->propGetInt(in, "fn", 0, &err);
	if (err)
		d.fn = 0;
	else
	{
		if (d.fn < 0 || d.fn >= d.vi->numFrames)
		{
			vsapi->setError(out, "variance: invalid fn. Not within clip");
			vsapi->freeNode(d.node);
			return;
		}
    }
	d.xgrid = vsapi->propGetInt(in, "xgrid", 0, &err);
	if (err)
		d.xgrid = 5;
	else
	{
		if (d.xgrid < 3 || d.xgrid >=  d.vi->width)
		{
			vsapi->setError(out, "variance: invalid xgrid. value must be 3 to width of frame");
			vsapi->freeNode(d.node);
			return;
		}
    }
	d.ygrid = vsapi->propGetInt(in, "ygrid", 0, &err);
	if (err)
		d.ygrid = 5;
	else
	{
		if (d.ygrid < 3 || d.ygrid >=  d.vi->height)
		{
			vsapi->setError(out, "variance: invalid ygrid. value must be 3 to height of frame");
			vsapi->freeNode(d.node);
			return;
		}
    }
	temp = !!vsapi->propGetInt(in, "uv", 0, &err);
    if (err)
        d.UV = true;
	else
	{
		if (temp < 0 || temp >  1)
		{
			vsapi->setError(out, "variance: invalid uv. value must be 0 or 1");
			vsapi->freeNode(d.node);
			return;
		}
	}
    
	
    // I usually keep the filter data struct on the stack and don't allocate it
    // until all the input validation is done.
    data = (VarianceData *)malloc(sizeof(d));
    *data = d;

    // Creates a new filter and returns a reference to it. Always pass on the in and out
    // arguments or unexpected things may happen. The name should be something that's
    // easy to connect to the filter, like its function name.
    // The three function pointers handle initialization, frame processing and filter destruction.
    // The filtermode is very important to get right as it controls how threading of the filter
    // is handled. In general you should only use fmParallel whenever possible. This is if you
    // need to modify no shared data at all when the filter is running.
    // For more complicated filters, fmParallelRequests is usually easier to achieve as it can
    // be prefetched in parallel but the actual processing is serialized.
    // The others can be considered special cases where fmSerial is useful to source filters and
    // fmUnordered is useful when a filter's state may change even when deciding which frames to
    // prefetch (such as a cache filter).
    // If your filter is really fast (such as a filter that only resorts frames) you should set the
    // nfNoCache flag to make the caching work smoother.
    vsapi->createFilter(in, out, "variance", varianceInit, varianceGetFrame, varianceFree, fmParallel, 0, data, core);
}

//////////////////////////////////////////
// Init

// This is the entry point that is called when a plugin is loaded. You are only supposed
// to call the two provided functions here.
// configFunc sets the id, namespace, and long name of the plugin (the last 3 arguments
// never need to be changed for a normal plugin).
//
// id: Needs to be a "reverse" url and unique among all plugins.
// It is inspired by how android packages identify themselves.
// If you don't own a domain then make one up that's related
// to the plugin name.
//
// namespace: Should only use [a-z_] and not be too long.
//
// full name: Any name that describes the plugin nicely.
//
// registerFunc is called once for each function you want to register. Function names
// should be PascalCase. The argument string has this format:
// name:type; or name:type:flag1:flag2....;
// All argument name should be lowercase and only use [a-z_].
// The valid types are int,float,data,clip,frame,func. [] can be appended to allow arrays
// of type to be passed (numbers:int[])
// The available flags are opt, to make an argument optional, empty, which controls whether
// or not empty arrays are accepted
/*
VS_EXTERNAL_API(void) VapourSynthPluginInit(VSConfigPlugin configFunc, VSRegisterFunction registerFunc, VSPlugin *plugin)
{
    configFunc("in.vcmohan.clear", "clear", "VapourSynth Variance Example", VAPOURSYNTH_API_VERSION, 1, plugin);
    registerFunc("Variance", "clip:clip;lx:int;wd:int;ty:int;ht:int;fn:int:opt;uv:int:opt;xgrid:int:opt;ygrid:int:opt;", varianceCreate, 0, plugin);
}

*/