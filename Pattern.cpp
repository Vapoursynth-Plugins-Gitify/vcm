/******************************************************************************
The Pattern Filter creates several Patterns useful for testing filters
These patterns are overlaid on the  input clip.
Author V.C.Mohan
Created on 31 july 2020, 31 Jan 2021

Copyright (C) < 2020>  <V.C.Mohan>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is at
http://www.gnu.org/licenses/.


Author : V.C.Mohan


*************************************************************************************************/
//#include <stdlib.h>
//#include "VapourSynth.h"
//#include "VSHelper.h"


typedef struct {
	VSNodeRef *node; 
	const VSVideoInfo *vi;

	int type;				// type of Pattern 1.dirac, 2.disc, 3.zcos, 4.sine, 5.step
	int x;					//  x coord of center	
	int y;					//  y coord of center
	int rad;				// radius of disc pattern
	int wl;					// wave length of step, zcos ans sine patterns
	int orient;				// orientation of step or zcos or sine patterns  (1.circ, 2.vert, 3.hor )
	bool spk;				// For circle is center value to be spiked?
	float spike;			// value of spike
	float overlay;			// %age of overlaying of pattern
	bool stat;				// is the pattern static frame to frame?
	unsigned char color[3];	// color of pattern bgr or yuv computed value
	bool vert, hor, circ, slant;
	float * overlay_table;
}PatternData;
	
/***************************************************************/


static void VS_CC patternInit(VSMap *in, VSMap *out, void **instanceData, VSNode *node, VSCore *core, const VSAPI *vsapi)
{
	PatternData *d = (PatternData *)* instanceData;
	vsapi->setVideoInfo(d->vi, 1, node);

	d->overlay_table = NULL;
	d->vert = false;
	d->hor = false;
	d->circ = false;
	d->slant = false;
	// for zcos sine, step
	if (d->type > 2) 
	{
		if (d->orient == 2)
		{
			d->vert = true;
		}

		else if (d->orient == 3)
	 	{
			d->hor = true;
		}
		else if (d->orient == 4)
		{
			d->vert = true;
			d->slant = true;
		}

		else // if (d->orient == 1)
		{
			d->vert = true;
			d->hor = true;
			d->circ = true;
		}
	
		d->overlay_table = (float *)vs_aligned_malloc<float>(sizeof(float) * d->wl, 32);

		for (int i = 0; i < d->wl; i++)
		{
			if (d->type  == 4)
				d->overlay_table[i] = 1.0f - (((1.0f + sin(2 * M_PI * i / d->wl)) / 2.0f) *  d->overlay);		

			else if (d->type == 3)

				d->overlay_table[i] = 1.0f - (((1.0f +  cos(2 * M_PI * i / d->wl)) /2.0f) *  d->overlay);
		
			else  if (d->type == 5)  //step		

				d->overlay_table[i] = i < d->wl / 2 ? 1.0f - d->overlay : 1.0f;
		}

	} // dtype > 2  for zcos, sine, step
	
}

/***************************************************************/
// This is the main function that gets called when a frame should be produced. It will, in most cases, get
// called several times to produce one frame. This state is being kept track of by the value of
// activationReason. The first call to produce a certain frame n is always arInitial. In this state
// you should request all the input frames you need. Always do it in ascending order to play nice with the
// upstream filters.
// Once all frames are ready, the filter will be called with arAllFramesReady. It is now time to
// do the actual processing.
static const VSFrameRef *VS_CC patternGetFrame(int n, int activationReason, void **instanceData, void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi)
{
	PatternData *d = (PatternData *)* instanceData;

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
		VSFrameRef *dst = vsapi->copyFrame(src, core);
		int nplanes = fi->numPlanes;
		int nbytes = fi->bytesPerSample;
		int nbits = fi->bitsPerSample;
		// free src
		vsapi->freeFrame(src);
	
		unsigned char * dp[] = { NULL, NULL, NULL, NULL };
		int pitch[] = { 0,0,0,0 };
		int subH[] = { 0, fi->subSamplingH, fi->subSamplingH, 0 };
		int subW[] = { 0, fi->subSamplingW, fi->subSamplingW, 0 };			

		for (int p = 0; p < nplanes; p++)
		{			
			pitch[p] = vsapi->getStride(dst, p) / nbytes;
			dp[p] = vsapi->getWritePtr(dst, p);
		}
		
		if (d->type == 1)
		{
			
			for (int p = 0; p < nplanes; p++)
			{
				if (nbytes == 1)
				{
					*(dp[p] + (d->y >> subH[p]) * pitch[p] + (d->x >> subW[p])) = d->color[p];

				}
				else if (nbytes == 2)
				{
					*((uint16_t*)(dp[p]) + (d->y >> subH[p]) * pitch[p] + (d->x >> subW[p])) = d->color[p] << (nbits - 8);

				}

				else if (nbytes == 4)
				{
					if (p == 0 || fi->colorFamily == cmRGB)

						*((float*)(dp[p]) + (d->y >> subH[p]) * pitch[p] + (d->x >> subW[p])) = d->color[p] / 256.0f;
					else
						*((float*)(dp[p]) + (d->y >> subH[p]) * pitch[p] + (d->x >> subW[p])) = (d->color[p] - 127) / 256.0f;

				}
			}
			
		}	// dirac

		else if (d->type == 2)
		{
			// disc
			unsigned char  col2[3];

			if (d->spk)
			{
				col2[0] = (1.0f - d->spike) * d->color[0];

				if (fi->colorFamily == cmRGB)
				{
					col2[1] = (1.0f - d->spike) * d->color[1];
					col2[2] = (1.0f - d->spike) * d->color[2];
				}
				else //	 col
				{
					col2[1] = d->color[1];
					col2[2] = d->color[2];
				}
			}

			else
			{
				col2[0] = d->color[0];
				col2[1] = d->color[1];
				col2[2] = d->color[2];
			}

			int rsq = d->rad * d->rad;
			int startH = d->y - d->rad < 0 ? 0 : d->y - d->rad;
			int endH = d->y + d->rad > height ? height : d->y + d->rad;
			int startW = d->x - d->rad < 0 ? 0 : d->x - d->rad;
			int endW = d->x + d->rad > width ? width : d->x + d->rad;

			for (int h = startH; h < endH; h++)
			{
				int hsq = (h - d->y) * (h - d->y);

				for (int w = startW; w < endW; w++)
				{
					int wsq = (w - d->x) * (w - d->x);

					if (hsq + wsq <= rsq)
					{
						for (int p = 0; p < nplanes; p++)
						{
							if (nbytes == 1)

								*(dp[p] + (h >> subH[p]) * pitch[p]
								+ (w >> subW[p])) = col2[p];

							else if (nbytes == 2)

								*((uint16_t*)(dp[p]) + (h >> subH[p]) * pitch[p]
								+ (w >> subW[p])) = col2[p] << (nbits - 8);

							else if (nbytes == 4)
							{
								if (p == 0 || fi->colorFamily == cmRGB)

									*((float*)(dp[p]) + (h >> subH[p]) * pitch[p]
									+ (w >> subW[p])) = col2[p] / 256.0f;

								else

									*((float*)(dp[p]) + (h >> subH[p]) * pitch[p]
									+ (w >> subW[p])) = (col2[p] - 127) / 256.0f;

							}

						}
					}
				}  // for w

			}  // for h
				// central spiked value

			if (d->spk)
			{
				for (int p = 0; p < nplanes; p++)
				{

					if (nbytes == 1)
						*(dp[p] + (d->y >> subH[p]) * pitch[p]
						+ (d->x >> subW[p])) = d->color[p];

					else if (nbytes == 2)
						*((uint16_t*)(dp[p]) + (d->y >> subH[p]) * pitch[p]
						+ (d->x >> subW[p])) = d->color[p] << (nbits - 8);

					else if (nbytes == 4)
					{
						if (p == 0)

							*((float*)(dp[p]) + (d->y  * pitch[p])
							+ d->x) = d->color[p] / 256.0f;

						else
							*((float*)(dp[p]) + (d->y >> subH[p]) * pitch[p]
							+ (d->x >> subW[p])) = (d->color[p] - 127) / 256.0f;
					}
				}
			}			

		}	// else disc

		else  //if step or sine or zcos
		{
			// add this to create motion illusion in case of step, sine and zcos
			int stsq = d->stat ? 0 : n;
			bool step = d->type == 5 ? true : false;
			bool zcos = d->type == 3 ? true : false;
			bool sine = d->type == 4 ? true : false;		

			int wlsq = (d->wl * d->wl) / 4;	// wavelength / 2 squared

			float val;

			for (int h = 0; h < height; h++)
			{
				int hsq = d->hor ? (d->y - h) * (d->y - h) : 0;	// 

				for (int w = 0; w < width; w++)
				{

					int wsq = d->vert ? (d->x - w) * (d->x - w) : 0;

					int dist = ((hsq == 0 ? abs(d->x - w) : wsq == 0 ? abs(d->y - h) : sqrt(float(wsq + hsq))) + stsq);
					dist = zcos ? (wsq + hsq)/wlsq  + stsq : dist;
					if (d->slant)
						dist += h;

					val = d->overlay_table[dist % d->wl];

					for (int p = 0; p < nplanes; p++)
					{
						if (fi->colorFamily == cmYUV && p > 0)  continue;

						if (nbytes == 1)
						{
							*(dp[p]  + (w)) *= val;
						}
						else if (nbytes == 2)
						{
							*((uint16_t *)dp[p]  + (w)) *= val;
						}

						else  // float							
						{
							*((float *)dp[p] +  (w)) *= val;
						}
					}
	
				}	// for w

				dp[0] += pitch[0] * nbytes;
				dp[1] += pitch[1] * nbytes;
				dp[2] += pitch[2] * nbytes;

			}	// for h

		}	// type step sine zcos	

		return dst;
	}
	return 0;
}


/***************************************************************/
// Free all allocated data on filter destruction
static void VS_CC patternFree(void *instanceData, VSCore *core, const VSAPI *vsapi)
{
	PatternData *d = (PatternData *)instanceData;
	vsapi->freeNode(d->node);
	if (d->overlay_table != NULL)
		vs_aligned_free(d->overlay_table);
	free(d);
}

// This function is responsible for validating arguments and creating a new filter
static void VS_CC patternCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi)
{
	PatternData d;
	//   PatternData *data;
	int err;
	int temp;

	// Get a clip reference from the input arguments. This must be freed later.
	d.node = vsapi->propGetNode(in, "clip", 0, 0);
	d.vi = vsapi->getVideoInfo(d.node);

	// In this integer and float. Note that
	// vi->format can be 0 if the input clip can change format midstream.
	if (!isConstantFormat(d.vi) )
	{
		vsapi->setError(out, "Pattern: format of clip must be constant ");
		vsapi->freeNode(d.node);
		return;
	}
	if (d.vi->format->colorFamily != cmRGB && d.vi->format->colorFamily != cmYUV && d.vi->format->colorFamily != cmGray)
	{
		vsapi->setError(out, "Pattern: RGB, YUV and Gray color formats only for input allowed ");
		vsapi->freeNode(d.node);
		return;
	}
	if (d.vi->format->sampleType == stFloat && d.vi->format->bitsPerSample == 16)
	{
		vsapi->setError(out, "Pattern: Half float formats not allowed ");
		vsapi->freeNode(d.node);
		return;
	}

	d.type = vsapi->propGetInt(in, "type", 0, &err);

	if (err)
		d.type = 4;

	else if ( d.type < 1 || d.type > 5)
	{
		vsapi->setError(out, "Pattern: type can be 1 to 5 only. 1.dirac, 2 disc, 3.zcos, 4.sine or 5.step only ");
		vsapi->freeNode(d.node);
		return;
	}
	d.x = vsapi->propGetInt(in, "x", 0, &err);
	if (err)
		d.x = d.vi->width / 2;

	d.y = vsapi->propGetInt(in, "y", 0, &err);
	if (err)
		d.y = d.vi->height / 2;

	if (d.type == 1 || d.type == 2)
	{
		if (d.x < 0 || d.x >= d.vi->width || d.y < 0 || d.y >= d.vi->height)
		{
			vsapi->setError(out, "Pattern: for the type opted x and y must be within frame; ");
			vsapi->freeNode(d.node);
			return;
		}
		unsigned char bgr[3], yuv[3];

		for (int i = 0; i < 3; i++)
		{
			temp = vsapi->propGetInt(in, "bgr", i, &err);
			if (err)
			{
				if (i == 0)
					bgr[i] = 255;
				else
					bgr[i] = bgr[i - 1];
			}

			else if (temp < 0 || temp > 255)
			{
				vsapi->setError(out, "Pattern: bgr values must be between 0 and 255; ");
				vsapi->freeNode(d.node);
				return;
			}
			else
				bgr[i] = temp;
		}

		BGRtoYUV(bgr, yuv);

		for (int i = 0; i < 3; i++)
		{
			if (d.vi->format->colorFamily == cmYUV)

				d.color[i] = yuv[i];
			else
				d.color[i] = bgr[i];
		}	


		if (d.type == 2)
		{
			temp = vsapi->propGetInt(in, "rad", 0, &err);

			if (err)
				temp = 80;
			if (temp < 4 || temp  > d.vi->width || temp > d.vi->height)
			{
				vsapi->setError(out, "Pattern: value of rad must be positive and with x, y values should ensure the circle do not go out of frame");
				vsapi->freeNode(d.node);
				return;
			}
			d.rad = temp;

			temp = !!vsapi->propGetInt(in, "spk", 0, &err);
			if (err)
				d.spk = false;
			else
				d.spk = temp == 0 ? false : true;

			if (d.spk)
			{
				d.spike = vsapi->propGetFloat(in, "spike", 0, &err);
				if (err)
					d.spike = 0.1f;
				else if (d.spike <= 0.01f || d.spike > 0.99f)
				{
					vsapi->setError(out, "Pattern: spike must be between 0.01 and 0.99 ");
					vsapi->freeNode(d.node);
					return;
				}
			}
		}	// if disc
		
	}	// if disc or dirac
	else if (d.type > 2)
	{
		//3.zcos, 4.sine, 5.step
		d.orient = vsapi->propGetInt(in, "orient", 0, &err);
		if (err)
			d.orient = 1;
		else if (d.orient < 1 || d.orient > 4)
		{
			vsapi->setError(out, "Pattern: orient can be 1 to 4 for 1 circ, 2 vert, 3 hor, 4 slant ");
			vsapi->freeNode(d.node);
			return;
		}

		temp = !!vsapi->propGetInt(in, "stat", 0, &err);
		if (err)
			d.stat = true;
		else
			d.stat = temp == 0 ? false : true;

		temp = vsapi->propGetInt(in, "wl", 0, &err);

		if (err)
			temp = 16;
		else if (temp < 4 || temp > d.vi->height / 2 || temp > d.vi->width / 2)
		{
			vsapi->setError(out, "Pattern: wavelength wl should be between 4 and half of smaller dimension of frame ");
			vsapi->freeNode(d.node);
			return;
		}
		else
			d.wl = temp;

		d.overlay = vsapi->propGetFloat(in, "overlay", 0, &err);

		if (err)
			d.overlay = 0.1f;
		else if (d.overlay < 0.004 || d.overlay > 1.0)
		{
			vsapi->setError(out, "Pattern: overlay should be between 0.004 and 1.0 ");
			vsapi->freeNode(d.node);
			return;
		}
	}	// if step or sine or zcos

	else
	{
		// shoud not vome here
		vsapi->setError(out, "Pattern:unexpected error");
		vsapi->freeNode(d.node);
		return;
	}
	
	PatternData * data = (PatternData *)malloc(sizeof(d));

	*data = d;


	/***************************************************************/
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
	vsapi->createFilter(in, out, "Pattern", patternInit, patternGetFrame, patternFree, fmParallel, 0, data, core);
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
configFunc("in.vcmohan.misc", "vcMisc", "VapourSynth Pattern plugin", VAPOURSYNTH_API_VERSION, 1, plugin);
registerFunc("Pattern", "clip:clip;type:data:opt;orient:data:opt;spk:int:opt;spike:int:opt;wl:int:opt;x:int:opt,y:int:opt;rad:int:opt;stat:int:opt;overlay:int:opt;color:int[]:opt;", PatternCreate, 0, plugin);
}
*/

