/*

// Rotate vapoursynth plugin

// This file contains a  rotate
// filter. It rotates the given window of the frame around given axial coordinates
// through specified angle
// Lanczos, bicubic, bilinear interpolation or near point approximation is an option..
// 20 Aug 2020
// Copyright (C) <2006, 2020>  <V.C.Mohan>

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

  Oct 2014, july 2020
*/
/*
#include <stdlib.h>
#include "math.h"
#include "VapourSynth.h"
#include "VSHelper.h"
#include "interpolationMethods.h"

*/

typedef struct {
		 VSNodeRef *node[2];
		 const VSVideoInfo *vi[2];
			
		float dd;	// initial  angle of rotation
		float dinc;
		int lx, wd, ty, ht;		// coordinates of rectangle  to be rotated
		int axx, axy;			// axis coordinates around which frame window rotates
		int pquant;			// quantile precision		
	
		int intq;	// 0 nearest point, 1 bilinear 2 Cubic 4x4 pts 3, lanczos 6x6 
		int span;		// 6 x6 Lanczos, or 4x4 cubic or 2x2 linear interpolation or 0 for nearest neighbour
		float * lbuf;
	
	
} RotateData;
//--------------------------------------------------------------------------------------------------------
// This function is called immediately after vsapi->createFilter(). This is the only place where the video
// properties may be set. In this case we simply use the same as the input clip. You may pass an array
// of VSVideoInfo if the filter has more than one output, like rgb+alpha as two separate clips.
static void VS_CC rotateInit(VSMap *in, VSMap *out, void **instanceData,
	VSNode *node, VSCore *core, const VSAPI *vsapi) 
{
    RotateData *d = (RotateData *) * instanceData;
    vsapi->setVideoInfo(d->vi[0], 1, node);
	
		// set up all flags
	d->lbuf = NULL;	

	d->pquant = 32;	
	
	if ( d->intq == 3 )
	{
		d->span = 6;
			//lanczos  opted 6 X 6 points per interpolation
		d->lbuf = (float*) vs_aligned_malloc<float> ((d->pquant + 1) * d->span * sizeof(float), 32);		

		LanczosCoeff( d->lbuf, d->span, d->pquant);
	}

	else if ( d->intq == 2  )
	{
		d->span = 4;
			//cubic 4 X 4 opted or needed ht lanczos when we have only
			//4 points available for interpolating a value
		d->lbuf = (float*)vs_aligned_malloc<float>((d->pquant + 1) * d->span * sizeof(float), 32);

		CubicIntCoeff( d->lbuf, d->pquant);
	}

	else if( d->intq == 1)
	{
		d->span = 2;

		d->lbuf = (float*)vs_aligned_malloc<float>((d->pquant + 1) * d->span * sizeof(float), 32);

		LinearIntCoeff( d->lbuf, d->pquant);
	}
	else
	{
		d->span = 0;
		d->lbuf = NULL;
	}
}

// This is the main function that gets called when a frame should be produced. It will, in most cases, get
// called several times to produce one frame. This state is being kept track of ht the value of
// activationReason. The first call to produce a certain frame n is always arInitial. In this state
// you should request all the input frames you need. Always do it in ascending order to play nice with the
// upstream filters.
// Once all frames are ready, the filter will be called with arAllFramesReady. It is now time to
// do the actual processing.
static const VSFrameRef *VS_CC rotateGetFrame(int n, int activationReason, void **instanceData, void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) {
    RotateData *d = (RotateData *) * instanceData;

    if (activationReason == arInitial)
	{
        // Request the source frame on the first call
        vsapi->requestFrameFilter(n, d->node[0], frameCtx);
		vsapi->requestFrameFilter(n, d->node[1], frameCtx);

    }
	else if (activationReason == arAllFramesReady)
	{
        const VSFrameRef *src = vsapi->getFrameFilter(n, d->node[0], frameCtx);
        // The reason we query this on a per frame basis is because we want our filter
        // to accept clips with varying dimensions. If we reject such content using d->vi
        // would be better.
        const VSFormat *fi = d->vi[0]->format;
        int height = vsapi->getFrameHeight(src, 0);
        int width = vsapi->getFrameWidth(src, 0);

		const VSFrameRef *bkg = vsapi->getFrameFilter(n, d->node[1], frameCtx);
        // When creating a new frame for output it is VERY EXTREMELY SUPER IMPORTANT to
        // supply the "dominant" source frame to copy properties from. Frame props
        // are an essential part of the filter chain and you should NEVER break it.
        VSFrameRef *dst = vsapi->copyFrame(bkg, core);	
		
		int subW[] = { 0, fi->subSamplingW, fi->subSamplingW, 0 };
		int subH[] = { 0, fi->subSamplingH, fi->subSamplingH, 0 };
		uint8_t * dstp[] = { NULL, NULL, NULL, NULL };
		const uint8_t * srcp[] = { NULL, NULL, NULL, NULL };
		const uint8_t * bkgp[] = { NULL, NULL, NULL, NULL };
		
		int nbits = fi->bitsPerSample;
		int nbytes = fi->bytesPerSample;

		for (int plane = 0; plane < fi->numPlanes; plane++)
		{
			srcp[plane] = vsapi->getReadPtr(src, plane);
			int bpitch = vsapi->getStride(bkg, plane);
			dstp[plane] = vsapi->getWritePtr(dst, plane);
			int ht = vsapi->getFrameHeight(bkg, plane);
			int wd = vsapi->getFrameWidth(bkg, plane);
			bkgp[plane] = vsapi->getReadPtr(bkg, plane);
			
		}

		// compute rotational coordinates for the current frame
		float dd = d->dd + n * d->dinc;		

		if (dd > 0)
		{
			while ( dd >= 360 )
			{
				dd = dd - 360;
			}
		}

		else
		{
			while ( dd < 0 )
			{
				dd = dd + 360;
			}
		}
		float sinalfa = sin( (dd * M_PI ) / 180);
		float cosalfa = cos( (dd * M_PI ) / 180);

		for ( int h = d->ty; h < d->ty + d->ht; h ++)
		{
			float hsinalfa = ( h - d->axy) * sinalfa;
			float hcosalfa = ( h - d->axy) * cosalfa;

			for( int w = d->lx; w < d->lx + d->wd; w ++)
			{
				// h & w are coordinates after rotation. get their coordinates prior to rotation. 
				float newx = d->axx + ( w - d->axx) * cosalfa - hsinalfa;

				int intx =  (int) newx;			

				float newy = d->axy + (w - d->axx) * sinalfa + hcosalfa ;

				int inty = (int)newy;
				// check within frame
				if (intx >= 0 && intx < width && inty >= 0 && inty < height)
				{
									// we are digitizing the fraction
					// which quantile the fraction is
					int qx = (newx - intx) * d->pquant;
					int qy = (newy - inty) * d->pquant;

					for (int plane = 0; plane < fi->numPlanes; plane++)
					{
						int pitch = vsapi->getStride(dst, plane) / nbytes;
						int pht = vsapi->getFrameHeight(src, plane);
						int pwd = vsapi->getFrameWidth(src, plane);
						
						bool useNearPoint = true;	// flag for border pixels

						if (  subH[plane] == 0 && subW[plane] == 0) 
						{
							if (d->intq > 0)
							{
								useNearPoint = false;
								// check if there are enough pixels around for  interpolation.
								if (intx >= d->span / 2 && intx < pwd - d->span / 2
									&& inty >= d->span / 2 && inty < pht - d->span / 2)
								{
									useNearPoint = false;
									//  interpolation
									if (fi->sampleType == stInteger)
									{
										if (nbits == 8)
										{
											uint8_t min = 0, max = (1 << nbits) - 1;
											uint8_t * dp = (uint8_t *)dstp[plane];
											const uint8_t * sp = (uint8_t *)srcp[plane];
											if (needNotInterpolate(sp + inty * pitch + intx, pitch, 1))
											{
												dp[h * pitch + w] = sp[inty * pitch + intx];
											}
											else
											{
												// get interpolated value
												dp[h * pitch + w] = clamp(LaQuantile(sp + inty * pitch + intx,
													pitch, d->span, qx, qy, d->lbuf), min, max);
											}

										}
										else	// 10 or 12  16 bit samples
										{
											uint16_t min = 0, max = (1 << nbits) - 1;
											uint16_t * dp = (uint16_t *)dstp[plane];
											const uint16_t * sp = (uint16_t *)srcp[plane];
											if (needNotInterpolate(sp + inty * pitch + intx, pitch, 1))
											{
												dp[h * pitch + w] = sp[inty * pitch + intx];
											}
											else
											{
												// get interpolated value
												dp[h * pitch + w] = clamp(LaQuantile(sp + inty * pitch + intx,
													pitch, d->span, qx, qy, d->lbuf), min, max);
											}
										}
									}
									else		// floating pt samples
									{
										float min = plane == 0 ? 0.0 : fi->colorFamily == cmRGB ? 0.0 : -0.5f;
										float max = plane == 0 ? 1.0 : fi->colorFamily == cmRGB ? 1.0 : 0.5f;
										float * dp = (float *)dstp[plane];
										const float * sp = (float *)srcp[plane];
										if (needNotInterpolate(sp + inty * pitch + intx, pitch, 1))
										{
											dp[h * pitch + w] = sp[inty * pitch + intx];
										}
										else
										{
											// get interpolated value
											dp[h * pitch + w] = clamp(LaQuantile(sp + inty * pitch + intx,
												pitch, d->span, qx, qy, d->lbuf), min, max);
										}
									}
								}	// sufficient pixels for interpolation
								else
									useNearPoint = true;
							} // if intq > 0
							
						}	// no subsampling present for this plane

						else
						{
							useNearPoint = true;
						}

						if (d->intq == 0 || useNearPoint)
						{
							// in case of subsampled data or q = 0, or border pixels use nearest point
							int nearx = int(newx + 0.5f);
							int neary = int(newy + 0.5f);

							if (nearx >= 0 && nearx < width && neary >= 0 && neary < height)
							{

								if (fi->sampleType == stInteger)
								{

									if (nbits == 8)
									{
										uint8_t * dp = dstp[plane];
										const uint8_t * sp = srcp[plane];
										dp[(h >> subH[plane]) * pitch + (w >> subW[plane])]
											= sp[(neary >> subH[plane]) * pitch + (nearx >> subW[plane])];
									}
									else	// 9 to 16 bit samples
									{
										uint16_t * dp = (uint16_t *)dstp[plane];
										const uint16_t * sp = (uint16_t *)srcp[plane];
										dp[(h >> subH[plane]) * pitch + (w >> subW[plane])]
											= sp[(neary >> subH[plane]) * pitch + (nearx >> subW[plane])];
									}
								}
								else		// floating pt samples
								{
									float * dp = (float *)dstp[plane];
									const float * sp = (float *)srcp[plane];
									dp[(h >> subH[plane]) * pitch + (w >> subW[plane])]
										= sp[(neary >> subH[plane]) * pitch + (nearx >> subW[plane])];
								}

							}// if nearx and neary in frame
						}	// intq  0 or useNearPoint						
					}	// for plane
				}	// if intx inty in frame
			} // for w
		}	// for h        

        // Release the source frame
        vsapi->freeFrame(src);
		vsapi->freeFrame(bkg);
        // A reference is consumed when it is returned, so saving the dst reference somewhere
        // and reusing it is not allowed.
        return dst;
    }

    return 0;
}

// Free all allocated data on filter destruction
static void VS_CC rotateFree(void *instanceData, VSCore *core, const VSAPI *vsapi) 
{
    RotateData *d = (RotateData *)instanceData;

    vsapi->freeNode(d->node[0]);
	vsapi->freeNode(d->node[1]);

	if (d->lbuf != NULL)
		vs_aligned_free(d->lbuf);
    free(d);
}

// This function is responsible for validating arguments and creating a new filter
static void VS_CC rotateCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi) {
    RotateData d;
    RotateData *data;
    int err;

    // Get a clip reference from the input arguments. This must be freed later.
    d.node[0] = vsapi->propGetNode(in, "clip", 0, 0);
    d.vi[0] = vsapi->getVideoInfo(d.node[0]);

    // In this first version we only want to handle 8bit integer formats. Note that
    // vi->format can be 0 if the input clip can change format midstream.
    if (!isConstantFormat(d.vi[0] ) ) 
	{
        vsapi->setError(out, "rotate: only constant format input supported");
        vsapi->freeNode(d.node[0]);
        return;
    }

	d.node[1] = vsapi->propGetNode(in, "bkg", 0, 0);
    d.vi[1] = vsapi->getVideoInfo(d.node[1]);

    // In this first version we only want to handle 8bit integer formats. Note that
    // vi->format can be 0 if the input clip can change format midstream.
    if (!isSameFormat(d.vi[0] , d.vi[1]))
	{
        vsapi->setError(out, "rotate: background clip bkg must have same format as main clip");
        vsapi->freeNode(d.node[0]);
		vsapi->freeNode(d.node[1]);
        return;
    }
    // If a property read fails for some reason (index out of bounds/wrong type)
    // then err will have flags set to indicate why and 0 will be returned. This
    // can be very useful to know when having optional arguments. Since we have
    // strict checking because of what we wrote in the argument string, the only
    // reason this could fail is when the value wasn't set ht the user.
    // And when it's not set we want it to default to enabled.
	d.dd = vsapi->propGetFloat(in, "angle", 0, 0);

	d.dinc = vsapi->propGetFloat(in, "dinc", 0, &err);
	if(err)
	{
		d.dinc = 0.0;
	}

	d.lx = vsapi->propGetInt(in, "lx", 0, &err);
	if(err)
		d.lx = 0;
	if( d.lx < 0 || d.lx > d.vi[0]->width - 2)
	{
        vsapi->setError(out, "rotate: lx must be within clip and not more than frame width - 2");
        vsapi->freeNode(d.node[0]);
		vsapi->freeNode(d.node[1]);
        return;
    }

	d.wd = vsapi->propGetInt(in, "wd", 0, &err);
	if(err)
		d.wd = d.vi[0]->width - d.lx;
	if( d.wd < 2 || d.wd > d.vi[0]->width -  d.lx)
	{
        vsapi->setError(out, "rotate: wd must be atleast 2 and lx + wd within clip width");
        vsapi->freeNode(d.node[0]);
		vsapi->freeNode(d.node[1]);
        return;
    }
	d.ty = vsapi->propGetInt(in, "ty", 0, &err);
	if(err)
		d.ty = 0;
	if( d.ty < 0 || d.ty > d.vi[0]->height - 2)
	{
        vsapi->setError(out, "rotate: ty must be within clip and not more than frame height - 2");
        vsapi->freeNode(d.node[0]);
		vsapi->freeNode(d.node[1]);
        return;
    }

	d.ht = vsapi->propGetInt(in, "ht", 0, &err);
	if(err)
		d.ht = d.vi[0]->height - d.ty;
	if( d.ht < 2 || d.ht > d.vi[0]->height -  d.ty)
	{
        vsapi->setError(out, "rotate: ht must be atleast 2 and also ensure ty + ht not more than frame height");
        vsapi->freeNode(d.node[0]);
		vsapi->freeNode(d.node[1]);
        return;
    }
	d.axx = vsapi->propGetInt(in, "axx", 0, &err);
	if(err)
		d.axx = d.lx + d.wd/2;

	d.axy = vsapi->propGetInt(in, "axy", 0, &err);
	if(err)
		d.axy = d.ty + d.ht/2;
	
    d.intq = !!vsapi->propGetInt(in, "intq", 0, &err);
    if (err)
        d.intq = 2;

    // Let's pretend the only allowed values are 1 or 0...
    if (d.intq < 0 || d.intq > 3) {
        vsapi->setError(out, "rotate: intq must be 0 for near point, or 1 for bilinear or 2 bicubic or 3 for Lanczos interpolation");
        vsapi->freeNode(d.node[0]);
		vsapi->freeNode(d.node[1]);
        return;
    }

    // I usually keep the filter data struct on the stack and don't allocate it
    // until all the input validation is done.
    data = (RotateData *) malloc(sizeof(d));
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
    vsapi->createFilter(in, out, "Rotate", rotateInit, rotateGetFrame, rotateFree, fmParallel, 0, data, core);
}

//////////////////////////////////////////
// Init

// This is the entry point that is called when a plugin is loaded. You are only supposed
// to call the two provided functions here.
// configFunc sets the id, namespace, and long name of the plugin (the last 3 arguments
// never need to be changed for a normal plugin).
//
// id: Needs to be a "reverse" url and unique among all plugins.
// It is inspired ht how android packages identify themselves.
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
VS_EXTERNAL_API(void) VapourSynthPluginInit(VSConfigPlugin configFunc, VSRegisterFunction registerFunc, VSPlugin *plugin) {
    configFunc("in.vcmohan.move", "move", "VapourSynth Rotate Plugin", VAPOURSYNTH_API_VERSION, 1, plugin);
    registerFunc("Rotate", "clip:clip;bkg:clip;angle:float;dinc:float:opt;lx:int:opt;wd:int:opt;ty:int:opt;ht:int:opt;axx:int:opt;axy:int:opt;intq:int:opt;", rotateCreate, 0, plugin);
}
*/