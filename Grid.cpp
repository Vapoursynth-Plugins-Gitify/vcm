/******************************************************************************
Draws a grid over frames for getting coordinates
grid line interval,bold line and very bold line positions and colors can be specified 
Created sep 2014, 9 july 2020

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
#include "GridHelper.cpp"

typedef struct {
    VSNodeRef *node;
    const VSVideoInfo *vi;

	int gridinterval;
	int bold;	//every bold*gridinterval value line will be bold
	int vbold;	// every vbold* bold * gridinterval value will be vbold
	int color,bcolor, vbcolor;	// color of grid, bold, vbold colors
	int style;	// 0 for left top origin grid, 1 for frame centered grid and 2 for axis rulers
	
} GridData;		
/***************************************************************/


static void VS_CC gridInit(VSMap *in, VSMap *out, void **instanceData, VSNode *node, VSCore *core, const VSAPI *vsapi)
{
    GridData *d = (GridData *) * instanceData;
    vsapi->setVideoInfo(d->vi, 1, node);
	
	
}
	// This is the main function that gets called when a frame should be produced. It will, in most cases, get
// called several times to produce one frame. This state is being kept track of by the value of
// activationReason. The first call to produce a certain frame n is always arInitial. In this state
// you should request all the input frames you need. Always do it in ascending order to play nice with the
// upstream filters.
// Once all frames are ready, the filter will be called with arAllFramesReady. It is now time to
// do the actual processing.
static const VSFrameRef *VS_CC gridGetFrame(int n, int activationReason, void **instanceData, void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) 
{
    GridData *d = (GridData *) * instanceData;

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

		uint8_t color[] = { 0,0,0 }, bcolor[] = { 0,0,0 }, vbcolor[] = { 0, 0, 0 };
		//	uint8_t yuvvalue [3], byuvvalue[3], vbyuvvalue[3];		

		if (fi->colorFamily == cmRGB)
		{
			getComponentsRGB( color, d->color);
			getComponentsRGB( bcolor, d->bcolor);
			getComponentsRGB( vbcolor, d->vbcolor);
		
		}

		else if( fi->colorFamily == cmYUV)
		{
			getYUVfromRGB(color, d->color);
			getYUVfromRGB(bcolor, d->bcolor);
			getYUVfromRGB(vbcolor, d->vbcolor);
		}

		 // When creating a new frame for output it is VERY EXTREMELY SUPER IMPORTANT to
        // supply the "dominant" source frame to copy properties from. Frame props
        // are an essential part of the filter chain and you should NEVER break it.
        VSFrameRef *dst = vsapi->newVideoFrame(fi, width, height, src, core);

        // It's processing loop time!
        // Loop over all the planes
		
		for (int plane = 0; plane < fi->numPlanes; plane++)
		{
			const uint8_t *srcp = vsapi->getReadPtr(src, plane);
			int spitch = vsapi->getStride(src, plane);
			uint8_t *dstp = vsapi->getWritePtr(dst, plane);
			int dpitch = spitch / fi->bytesPerSample;	// vsapi->getStride(dst, plane); // note that if a frame has the same dimensions and format, the stride is guaranteed to be the same. int dst_stride = src_stride would be fine too in this filter.
			// Since planes may be subsampled you have to query the height of them individually
			int bht = vsapi->getFrameHeight(src, plane);
			int bwd = vsapi->getFrameWidth(src, plane);

			// copy input on to output frame (RGB and YUY2. Planar U and V will do later)
			vs_bitblt(dstp, spitch, srcp, spitch, bwd * fi->bytesPerSample, bht);

			int subW = plane == 0 ? 0 : fi->colorFamily == cmYUV ? fi->subSamplingW : 0;

			int subH = plane == 0 ? 0 : fi->colorFamily == cmYUV ? fi->subSamplingH : 0;

			// for drawing bold and very bold lines

			int modbold = d->bold * d->gridinterval;

			int modvbold = modbold * d->vbold;

			if (d->style < 2)
			{
				if ( fi->sampleType == stInteger)
				{
					const int nb = fi->bitsPerSample ;

					if ( nb == 8)
					{
						uint8_t * wp = dstp;
						// horizontal grid lines. so spacing is h 

						if (d->style == 1)
						{
							// frame centered grid
							DrawCenteredGridLines(wp + dpitch * bht / 2, dpitch, bwd, bht/2, 1,
								modvbold >> subH, modbold >> subH, d->gridinterval >> subH,
								vbcolor[plane], bcolor[plane], color[plane], 4 >> subW, 2 >> subW);
							// vertical lines
							DrawCenteredGridLines(wp + bwd / 2, 1, bht, bwd/2, dpitch,
								modvbold >> subW, modbold >> subW, d->gridinterval >> subW,
								vbcolor[plane], bcolor[plane], color[plane], 4 >> subH, 2 >> subH);
						}
						else
						{
								// left top origin grid
							DrawGridLines(wp, dpitch, bwd, bht, 1,
								modvbold >> subH, modbold >> subH, d->gridinterval >> subH,
								vbcolor[plane], bcolor[plane], color[plane], 4 >> subW, 2 >> subW);
							// vertical lines
							DrawGridLines(wp, 1, bht, bwd, dpitch,
								modvbold >> subW, modbold >> subW, d->gridinterval >> subW,
								vbcolor[plane], bcolor[plane], color[plane], 4 >> subH, 2 >> subH);
						}
					}
					else // if nb > 8
					{
						uint16_t *wp = (uint16_t *)dstp;
						uint16_t vbc = vbcolor[plane] * ( 1 << (nb - 8) );
						uint16_t  bc =  bcolor[plane] * ( 1 << (nb - 8) );
						uint16_t   c =   color[plane] * ( 1 << (nb - 8) );

						if (d->style == 1)
						{
							DrawCenteredGridLines(wp + dpitch * bht / 2, dpitch, bwd, bht / 2, 1,
								modvbold >> subH, modbold >> subH, d->gridinterval >> subH,
								vbc, bc, c, 4 >> subW, 2 >> subW);
							// vertical lines
							DrawCenteredGridLines(wp + bwd / 2, 1, bht, bwd / 2, dpitch,
								modvbold >> subW, modbold >> subW, d->gridinterval >> subW,
								vbc, bc, c, 4 >> subW, 2 >> subW);
						}
						else
						{

							DrawGridLines(wp, dpitch, bwd, bht, 1,
								modvbold >> subH, modbold >> subH, d->gridinterval >> subH,
								vbc, bc, c, 4 >> subW, 2 >> subW);

							DrawGridLines(wp, 1, bht, bwd, dpitch,
								modvbold >> subW, modbold >> subW, d->gridinterval >> subW,
								vbc, bc, c, 4 >> subH, 2 >> subH);
						}
					}
				}	// integer

				else	// float
				{
					float * wp = (float *) dstp;
					float uv = plane == 0 ? 0 : 0.5;	// luma values 0 to 1.0, chroma -0.5 to 0.5
					float vbc = vbcolor[plane] / 256.0 - uv;
					float  bc =  bcolor[plane] / 256.0 - uv;
					float   c =   color[plane] / 256.0 - uv;

					if (d->style == 1)
					{
						DrawCenteredGridLines(wp + dpitch * bht / 2, dpitch, bwd, bht / 2, 1,
							modvbold >> subH, modbold >> subH, d->gridinterval >> subH,
							vbc, bc, c, 4 >> subW, 2 >> subW);
						// vertical lines
						DrawCenteredGridLines(wp + bwd / 2, 1, bht, bwd / 2, dpitch,
							modvbold >> subW, modbold >> subW, d->gridinterval >> subW,
							vbc, bc, c, 4 >> subW, 2 >> subW);
					}
					else
					{

						DrawGridLines(wp, dpitch, bwd, bht, 1,
							modvbold >> subH, modbold >> subH, d->gridinterval >> subH,
							vbc, bc, c, 4 >> subW, 2 >> subW);

						DrawGridLines(wp, 1, bht, bwd, dpitch,
							modvbold >> subW, modbold >> subW, d->gridinterval >> subW,
							vbc, bc, c, 4 >> subH, 2 >> subH);
					}
				}

			}

			else
			{
				//  if (d->style == 2)
			
				// axial rulers only

				if ( fi->sampleType == stInteger)
				{
					const int nb = fi->bitsPerSample ;

					if ( nb == 8)
					{
						uint8_t * wp = dstp;
				
				// draw vertical line
						drawVBoldLine(wp + bwd / 2 , bht, dpitch, vbcolor[plane]);	// luma
						
						// draw short horizontal lines across vertical line
						DrawAxisScale( wp + bwd / 2 , dpitch, bwd, bht, 1,			  
									modvbold >> subH, modbold >> subH, d->gridinterval >> subH,
									vbcolor[plane], bcolor[plane], color[plane],
									8 >> subW, 6 >> subW, 4 >> subW);		// lengths of long, medium and short ticks
							
						// horizontal axis line

						drawVBoldLine(wp + bht / 2 * dpitch, bwd, 1, vbcolor[plane]);	
						
						// draw short vertical lines across horizontalal line

						DrawAxisScale( wp + bht / 2 * dpitch, 1, bht, bwd, dpitch,			  
										modvbold >> subW, modbold >> subW, d->gridinterval >> subW,
										vbcolor[plane], bcolor[plane], color[plane],
										8 >> subH, 6 >> subH, 4 >> subH);		// lengths of long, medium and short ticks

						
					}
					else // nb > 8
					{
						uint16_t *wp = (uint16_t*) dstp;
						uint16_t vbc = vbcolor[plane] * ( 1 << (nb - 8) );
						uint16_t  bc =  bcolor[plane] * ( 1 << (nb - 8) );
						uint16_t   c =   color[plane] * ( 1 << (nb - 8) );

						// draw vertical line
						drawVBoldLine(wp + bwd / 2 , bht, dpitch, vbc);	// luma
				
						// draw short horizontal lines across vertical line
						DrawAxisScale( wp + bwd / 2 , dpitch, bwd, bht, 1,			  
							modvbold >> subH, modbold >> subH, d->gridinterval >> subH,
							vbc, bc, c,
							8 >> subW, 6 >> subW, 4 >> subW);		// lengths of long, medium and short ticks
					
						// horizontal axis line

						drawVBoldLine(wp + bht / 2 * dpitch, bwd, 1, vbc);	
				
						// draw short vertical lines across horizontalal line

						DrawAxisScale( wp + bht / 2 * dpitch, 1, bht, bwd, dpitch,			  
								modvbold >> subW, modbold >> subW, d->gridinterval >> subW,
								vbc, bc, c,
								8 >> subH, 6 >> subH, 4 >> subH);		// lengths of long, medium and short ticks

					}
				}

				else	// float
				{
					float * wp = (float *) dstp;
					float uv = plane == 0 ? 0 : 0.5;	// luma values 0 to 1.0, chroma -0.5 to 0.5
					float vbc = vbcolor[plane] / 256.0 - uv;
					float  bc =  bcolor[plane] / 256.0 - uv;
					float   c =   color[plane] / 256.0 - uv;

					// draw vertical line
					drawVBoldLine(wp + bwd / 2 , bht, dpitch, vbc);	// luma
				
						// draw short horizontal lines across vertical line
					DrawAxisScale( wp + bwd / 2 , dpitch, bwd, bht, 1,			  
							modvbold >> subH, modbold >> subH, d->gridinterval >> subH,
							vbc, bc, c,
							8 >> subW, 6 >> subW, 4 >> subW);		// lengths of long, medium and short ticks
					
						// horizontal axis line

					drawVBoldLine(wp + bht / 2 * dpitch, bwd, 1, vbc);	
				
					// draw short vertical lines across horizontalal line

					DrawAxisScale( wp + bht / 2 * dpitch, 1, bht, bwd, dpitch,			  
							modvbold >> subW, modbold >> subW, d->gridinterval >> subW,
							vbc, bc, c,
							8 >> subH, 6 >> subH, 4 >> subH);		// lengths of long, medium and short ticks
				}

						
			}	// if axis

		} // int plane
		 
	
			 vsapi->freeFrame(src);

			// A reference is consumed when it is returned, so saving the dst reference somewhere
			// and reusing it is not allowed.
			return dst;
	}

		return 0;
			
	

}

	// Free all allocated data on filter destruction
static void VS_CC gridFree(void *instanceData, VSCore *core, const VSAPI *vsapi) 
{
    GridData *d = (GridData *)instanceData;
    vsapi->freeNode(d->node);
    free(d);
}

// This function is responsible for validating arguments and creating a new filter
static void VS_CC gridCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi) 
{
    GridData d;
 //   GridData *data;
    int err;
	int temp[3];
	int m;
    // Get a clip reference from the input arguments. This must be freed later.
    d.node = vsapi->propGetNode(in, "clip", 0, 0);
    d.vi = vsapi->getVideoInfo(d.node);

    // In this integer and float. Note that
    // vi->format can be 0 if the input clip can change format midstream.
   // if (!isConstantFormat(d.vi) ) // || d.vi->format->sampleType != stInteger || d.vi->format->bitsPerSample != 8)
	

	d.gridinterval = vsapi->propGetInt(in, "lineint", 0, &err);

	if (err)
		d.gridinterval = 10;

	if (d.gridinterval < 4 || d.gridinterval > 100)
	{
		vsapi->setError(out, "Grid: lineint can have values between 4 and 100 only ");
        vsapi->freeNode(d.node);
        return;
	}
	d.bold = vsapi->propGetInt(in, "bold", 0, &err);

	if (err)
		d.bold = 5;

	if (d.bold < 1 )
	{
		vsapi->setError(out, "Grid: bold can he a number 1 or more ");
        vsapi->freeNode(d.node);
        return;
	}
	d.vbold = vsapi->propGetInt(in, "vbold", 0, &err);

	if (err)
		d.vbold = 2;

	if (d.vbold < 1 )
	{
		vsapi->setError(out, "Grid: vbold must be anuber 1 or more ");
        vsapi->freeNode(d.node);
        return;
	}

	d.color = vsapi->propGetInt(in, "color", 0, &err);

	if (err)
		d.color = 0;

	else
	{
		m = vsapi->propNumElements(in, "color");
		temp[0] = 0;
		temp[1] = 0;
		temp[2] = 0;

		for ( int i = 0; i < m; i ++)
		{
			temp[i] = vsapi->propGetInt(in, "color", i, &err);

			if ( temp[i] < 0 || temp[i] > 255)
			{
				vsapi->setError(out, "Grid: color parameter values must be between 0 and 255 ");
				vsapi->freeNode(d.node);
				return;
			}
		}
	
		d.color = (((temp[0] << 16 ) | temp[1] << 8 ) | temp[2]);
	}


	d.bcolor = vsapi->propGetInt(in, "bcolor", 0, &err);

	if (err)
		d.bcolor = d.color;

	else
	{
		m = vsapi->propNumElements(in, "bcolor");
		temp[0] = 0;
		temp[1] = 0;
		temp[2] = 0;

		for ( int i = 0; i < m; i ++)
		{
			temp[i] = vsapi->propGetInt(in, "bcolor", i, &err);

			if ( temp[i] < 0 || temp[i] > 255)
			{
				vsapi->setError(out, "Grid: bcolor parameter values must be between 0 and 255 ");
				vsapi->freeNode(d.node);
				return;
			}
		}
	
		d.bcolor = (((temp[0] << 16 ) | temp[1] << 8 ) | temp[2]);
	} 

	d.vbcolor = vsapi->propGetInt(in, "vbcolor", 0, &err);

	if (err)
		d.vbcolor = d.bcolor;

	else
	{
		m = vsapi->propNumElements(in, "vbcolor");
		temp[0] = 0;
		temp[1] = 0;
		temp[2] = 0;

		for ( int i = 0; i < m; i ++)
		{
			temp[i] = vsapi->propGetInt(in, "vbcolor", i, &err);

			if ( temp[i] < 0 || temp[i] > 255)
			{
				vsapi->setError(out, "Grid: vbcolor parameter values must be between 0 and 255 ");
				vsapi->freeNode(d.node);
				return;
			}
		}
	
		d.vbcolor = (((temp[0] << 16 ) | temp[1] << 8 ) | temp[2]);
	}	

	d.style = vsapi->propGetInt(in, "style", 0, &err);

	if (err)
		d.style = 0;

	if (d.style < 0 || d.style  > 2)
	{
		vsapi->setError(out, "Grid: style can have values between 0 for Left Top origin grid, 1 for frame centered grid and 2 for centred rulers only ");
        vsapi->freeNode(d.node);
        return;
	}

	

	GridData * data = (GridData *) malloc(sizeof(d));

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
    vsapi->createFilter(in, out, "Grid", gridInit, gridGetFrame, gridFree, fmParallel, 0, data, core);
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
    configFunc("in.vcmohan.grid", "vcm", "VapourSynth grid plugin", VAPOURSYNTH_API_VERSION, 1, plugin);
	registerFunc("Grid", "clip:clip;lineint:int:opt;bold:int:opt;vbold:int:opt;color:int[]:opt;bcolor:int[]:opt;vbcolor:int[]:opt;style:int:opt;", gridCreate, 0, plugin);
}
*/    