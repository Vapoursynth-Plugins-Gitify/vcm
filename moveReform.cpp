/*
	Reform filter plugin for vapoursynth by V.C.Mohan
	deskews a quadrilateral into a rectangle or skews
	a rectangle into a quadrilateral  

	Copyright (C) <2012 - 2020>  <V.C.Mohan>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    A copy of the GNU General Public License is at
    http://www.gnu.org/licenses/.
	
	sep, 2014, 30 Aug 2020
*/
/*
#include <stdlib.h>
#include "VapourSynth.h"
#include "VSHelper.h"
#include "UnitSq2Quad_matrix.cpp"
#include "interpolationMethods.h"
*/
typedef struct {
    VSNodeRef *node[2];
    const VSVideoInfo *vi[2];
	float rect[2][2];
	int intq;	// 0 near point, 1 bilinear, 2 bicubic, 3 Lanczos interpolation
    float *lbuf;// coefficients for interpolation
	float sq[3][3];		// square matrix
	float quad[4][2];
	float inv[3][3];
	bool q2r;

	int pquant;
	int span, quantiles;
	// convertion from unit square to frame coordinates
	int transx, transy;
	float scalex, scaley;
	// limits of search
	int sxmin, sxmax, symin, symax;
	int wmin, wmax, hmin, hmax;
} ReformData;


// This function is called immediately after vsapi->createFilter(). This is the only place where the video
// properties may be set. In this case we simply use the same as the input clip. You may pass an array
// of VSVideoInfo if the filter has more than one output, like rgb+alpha as two separate clips.
static void VS_CC reformInit(VSMap *in, VSMap *out, void **instanceData, VSNode *node, VSCore *core, const VSAPI *vsapi) {
	ReformData *d = (ReformData *)* instanceData;
	vsapi->setVideoInfo(d->vi[0], 1, node);

	// 

	int width = d->vi[0]->width;
	int height = d->vi[0]->height;

	d->lbuf = NULL;

	d->pquant = 16;	// fractions binned into quantiles
	// create LUT coefficients buffers for opted interpolation
	
	if (d->intq == 3)
	{
		d->span = 6;
		d->lbuf = vs_aligned_malloc <float>((d->pquant + 1) * d->span * sizeof(float), 32);
		LanczosCoeff(d->lbuf, d->span, d->pquant);
	}
	else if (d->intq == 2)
	{
		d->span = 4;
		d->lbuf = vs_aligned_malloc<float>((d->pquant + 1) * d->span * sizeof(float), 32);
		CubicIntCoeff(d->lbuf, d->pquant);
	}

	else if (d->intq == 1)
	{
		d->span = 2;
		d->lbuf = vs_aligned_malloc <float>((d->pquant + 1) * d->span * sizeof(float), 32);
		LinearIntCoeff(d->lbuf, d->pquant);
	}
	else // d->intq == 0 nearest point
	{
		d->span = 1;
		d->lbuf = NULL;
	}

	// generate forward and inverse matrices 
	if (d->q2r)
	{
		d->transx = -d->rect[0][0];	// lx

		d->transy = -d->rect[0][1];	// topy

		d->scalex = 1.0 / (d->rect[1][0] - d->rect[0][0]); // 1 / width of rect;

		d->scaley = 1.0 / (d->rect[1][1] - d->rect[0][1]); // 1 / ht;
	}

	else
	{
		d->transx = d->rect[0][0];	// lx

		d->transy = d->rect[0][1];	// topy

		d->scalex = d->rect[1][0] - d->rect[0][0]; // wd;

		d->scaley = d->rect[1][1] - d->rect[0][1]; // ht;
	}

	if (UnitSq2Quad(d->sq, d->inv, d->quad) != 0)
	{

		vsapi->setError(out, " reform: Un invertible Matrix. Check your params");
		vsapi->freeNode(d->node[0]);
		vsapi->freeNode(d->node[1]);
		if (d->lbuf != NULL)
			free(d->lbuf);
		free(d);
		return;
	}


	int minRectX = VSMAX(VSMIN(d->rect[0][0], d->rect[1][0]), 0);
	int maxRectX = VSMIN(VSMAX(d->rect[0][0], d->rect[1][0]), width - 1);
	int minRectY = VSMAX(VSMIN(d->rect[0][1], d->rect[1][1]), 0);
	int maxRectY = VSMIN(VSMAX(d->rect[0][1], d->rect[1][1]), height - 1);


	int minQuadX = d->quad[0][0], maxQuadX = d->quad[0][0], minQuadY = d->quad[0][1], maxQuadY = d->quad[0][1];

	for (int i = 1; i < 4; i++)
	{
		minQuadX = VSMIN(minQuadX, d->quad[i][0]);
		minQuadY = VSMIN(minQuadY, d->quad[i][1]);

		maxQuadX = VSMAX(maxQuadX, d->quad[i][0]);
		maxQuadY = VSMAX(maxQuadY, d->quad[i][1]);

	}

	minQuadX = VSMAX(minQuadX, 0);
	minQuadY = VSMAX(minQuadY, 0);
	maxQuadY = VSMIN(maxQuadY, height - 1);
	maxQuadX = VSMIN(maxQuadX, width - 1);

	// the following max and min is to remove restrictions on coordinate parameters
	if (!d->q2r)
	{
		d->hmax = maxQuadY;
		d->hmin = minQuadY;
		d->wmax = maxQuadX;
		d->wmin = minQuadX;
		d->sxmin = minRectX;
		d->sxmax = maxRectX;
		d->symin = minRectY;
		d->symax = maxRectY;
	}
	else
	{
		d->hmax = maxRectY;
		d->hmin = minRectY;
		d->wmax = maxRectX;
		d->wmin = minRectX;
		d->sxmin = minQuadX;
		d->sxmax = maxQuadX;
		d->symin = minQuadY;
		d->symax = maxQuadY;
	}
}	

//-----------------------------------------------------------------

//---------------------------------------------------------------------
// This is the main function that gets called when a frame should be produced. It will, in most cases, get
// called several times to produce one frame. This state is being kept track of by the value of
// activationReason. The first call to produce a certain frame n is always arInitial. In this state
// you should request all the input frames you need. Always do it in ascending order to play nice with the
// upstream filters.
// Once all frames are ready, the filter will be called with arAllFramesReady. It is now time to
// do the actual processing.
static const VSFrameRef *VS_CC reformGetFrame(int n, int activationReason, void **instanceData, void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) 
{
    ReformData *d = (ReformData *) * instanceData;

    if (activationReason == arInitial)
	{
        // Request the source frame on the first call
        vsapi->requestFrameFilter(n, d->node[0], frameCtx);
		vsapi->requestFrameFilter(n, d->node[1], frameCtx);

    } 
	else if (activationReason == arAllFramesReady)
	{
        const VSFrameRef *src = vsapi->getFrameFilter(n, d->node[0], frameCtx);
		
        const VSFormat *fi = d->vi[0]->format;
        int height = vsapi->getFrameHeight(src, 0);
        int width = vsapi->getFrameWidth(src, 0);

        const VSFrameRef *bkg = vsapi->getFrameFilter(n, d->node[1], frameCtx);
		
		VSFrameRef *dst = vsapi->copyFrame(bkg, core);
		// VSFrameRef *dst = vsapi->newVideoFrame(fi, width, height, src, core);
		
		int subW[] = { 0, fi->subSamplingW, fi->subSamplingW, 0 };
		int subH[] = { 0, fi->subSamplingH, fi->subSamplingH, 0 };
		uint8_t * dstp[] = { NULL, NULL, NULL, NULL };
		const uint8_t * srcp[] = { NULL, NULL, NULL, NULL };
		
		int nbits = fi->bitsPerSample;
		int nbytes = fi->bytesPerSample;

		for (int plane = 0; plane < fi->numPlanes; plane++)
		{
			srcp[plane] = vsapi->getReadPtr(src, plane);
			dstp[plane] = vsapi->getWritePtr(dst, plane);
		}
		
		float sourcex, sourcey;

		for (int h = d->hmin; h < d->hmax; h++)
		{
			for (int w = d->wmin; w < d->wmax; w++)
			{
				if (d->q2r)
				{
					float hh = (h + d->transy) * d->scaley;

					float ww = (w + d->transx) * d->scalex;

					float denom = d->sq[0][2] * ww + d->sq[1][2] * hh + 1.0;

					if (denom > -1e-10 && denom < 1e-10) continue;

					sourcex = (d->sq[0][0] * ww + d->sq[1][0] * hh + d->sq[2][0]) / denom;

					sourcey = (d->sq[0][1] * ww + d->sq[1][1] * hh + d->sq[2][1]) / denom;
				}

				else
				{
					// r2q

					// matrix inversion 
					float denom = d->inv[0][2] * w + d->inv[1][2] * h + d->inv[2][2];

					if (denom > -1e-10 && denom < 1e-10) continue;

					sourcex = (d->inv[0][0] * w + d->inv[1][0] * h + d->inv[2][0]) / denom;

					sourcex = sourcex * d->scalex + d->transx;

					sourcey = (d->inv[0][1] * w + d->inv[1][1] * h + d->inv[2][1]) / denom;

					sourcey = sourcey * d->scaley + d->transy;
				}

				// get integer and fraction values

				int sx =  (int)sourcex;

				int sy =  (int)sourcey;

				if (sx >= d->sxmin && sx <= d->sxmax && sy >= d->symin && sy <= d->symax)
				{
					// inside the frame
					// get quantile corresponding to fraction value

					int qx = (sourcex - sx) * d->pquant;
					int qy = (sourcey - sy) * d->pquant;

					for (int plane = 0; plane < fi->numPlanes; plane++)
					{
						int pitch = vsapi->getStride(dst, plane) / nbytes;
						int pht = vsapi->getFrameHeight(src, plane);
						int pwd = vsapi->getFrameWidth(src, plane);

						bool useNearPoint = true;	// flag for border pixels

						if (subH[plane] == 0 && subW[plane] == 0)
						{
							if (d->intq > 0)
							{
								
								// check if there are enough pixels around for  interpolation.
								if (sx >= d->span / 2 && sx < pwd - d->span / 2
									&& sy >= d->span / 2 && sy < pht - d->span / 2)
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

											if (needNotInterpolate(sp + sy * pitch + sx, pitch, 1))
											{
												dp[h * pitch + w] = sp[sy * pitch + sx];
											}
											else
											{
												// get interpolated value
												dp[h * pitch + w] = clamp(LaQuantile(sp + sy * pitch + sx,
													pitch, d->span, qx, qy, d->lbuf), min, max);
											}

										}
										else	// 10 or 12  16 bit samples
										{
											uint16_t min = 0, max = (1 << nbits) - 1;
											uint16_t * dp = (uint16_t *)dstp[plane];
											const uint16_t * sp = (uint16_t *)srcp[plane];

											if (needNotInterpolate(sp + sy * pitch + sx, pitch, 1))
											{
												dp[h * pitch + w] = sp[sy * pitch + sx];
											}
											else
											{
												// get interpolated value
												dp[h * pitch + w] = clamp(LaQuantile(sp + sy * pitch + sx,
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

										if (needNotInterpolate(sp + sy * pitch + sx, pitch, 1))
										{
											dp[h * pitch + w] = sp[sy * pitch + sx];
										}
										else
										{
											// get interpolated value
											dp[h * pitch + w] = clamp(LaQuantile(sp + sy * pitch + sx,
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
							// subsampled plane
							useNearPoint = true;
						}

						if (d->intq == 0 || useNearPoint)
						{
							// in case of subsampled data or q = 0, or border pixels use nearest point
							int nearx = int(sourcex + 0.5f);
							int neary = int(sourcey + 0.5f);

							if (nearx >= 0 && nearx < pwd && neary >= 0 && neary < pht)
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

				}	// if within max and min values
						
			}	// for w

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
static void VS_CC reformFree(void *instanceData, VSCore *core, const VSAPI *vsapi) {
    ReformData *d = (ReformData *)instanceData;
    vsapi->freeNode(d->node[0]);
	vsapi->freeNode(d->node[1]);
	if( d->lbuf != NULL)
		vs_aligned_free (d->lbuf);
	
    free(d);
}
//-----------------------------------------------------------------------------------------------------

// This function is responsible for validating arguments and creating a new filter
static void VS_CC reformCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi)
{
    ReformData d;
    ReformData *data;
    int err;
	int temp;
	bool norm = true;	

    // Get a clip reference from the input arguments. This must be freed later.
    d.node[0] = vsapi->propGetNode(in, "clip", 0, 0);
    d.vi[0] = vsapi->getVideoInfo(d.node[0]);

	d.node[1] = vsapi->propGetNode(in, "bkg", 0, 0);
    d.vi[1] = vsapi->getVideoInfo(d.node[1]);

    // 
    // vi->format can be 0 if the input clip can change format midstream.
    if (!isConstantFormat(d.vi[0]) || !isSameFormat(d.vi[0], d.vi[1] ) )
	{
        vsapi->setError(out, "reform: only constant format input supported. Both clips must have same format");
        vsapi->freeNode(d.node[0]);
		vsapi->freeNode(d.node[1]);
        return;
    }
	int height = d.vi[0]->height;
	int width = d.vi[0]->width;
    
    temp = !!vsapi->propGetInt(in, "norm", 0, &err);
    if (err)
        norm = true;
	else
	{
    // Let's pretend the only allowed values are 1 or 0...
		if (temp < 0 || temp > 1)
		{
			vsapi->setError(out, "reform:  values allowed for norm are 0 for normalized and 1 for absolute values as coordinates");
			vsapi->freeNode(d.node[0]);
			vsapi->freeNode(d.node[1]);
			return;
		}
	}
	if ( temp == 0)
		norm = false;
	else
		norm = true;

	bool soft = false;
	temp = !!vsapi->propGetInt(in, "soft", 0, &err);
	if (!err)
	{
		if (temp != 0)
			soft = true;
	}
		

	temp = vsapi->propNumElements(in, "rect");

	if (temp == 0)
	{
		d.rect[0][0] = 0;
		d.rect[0][1] = 0;
		if (norm)
		{
			d.rect[1][0] = 1.0f;
			d.rect[1][1] = 1.0f;
		}
		else
		{
			d.rect[1][0] = width - 1;
			d.rect[1][1] = height - 1;
		}
	}

	else
	{
		if(temp != 4)
		{
			vsapi->setError(out, "reform: array rect must have exactly 4 normalized coordinate values corresponding opposite cornrs of rectangle");
			vsapi->freeNode(d.node[0]);
			vsapi->freeNode(d.node[1]);
			return;
		}

		else
		{
			for(int i = 0; i < 2; i ++)
			{
				
				d.rect[i][0] = vsapi->propGetFloat(in, "rect", i + i , 0);
				d.rect[i][1] = vsapi->propGetFloat(in, "rect", i + i + 1 , 0);

				if( norm)
				{					
					d.rect[i][0] *= (width - 1);
					d.rect[i][1] *= (height - 1);
				}
			}
		}
	}
	

	if (d.rect[0][0] == d.rect[1][0] || d.rect[0][1] == d.rect[1][1])
	{		
		vsapi->setError(out, " reform: width or height of rect is zero.");
		vsapi->freeNode(d.node[0]);
		vsapi->freeNode(d.node[1]);
		return;
	}
	

	temp = vsapi->propNumElements(in, "quad");

	if ( temp != 8 )
	{
		vsapi->setError(out, "reform: array quad must have exactly 8 entries corresponding to 4 x, y coordinate pairs in a clockwise direction");
		vsapi->freeNode(d.node[0]);
		vsapi->freeNode(d.node[1]);
		return;
	}
	else
	{
		for(int i = 0; i < 4; i ++)
		{
			
			d.quad[i][0] = vsapi->propGetFloat(in, "quad", i + i + 0, 0);
			d.quad[i][1] = vsapi->propGetFloat(in, "quad", i + i + 1, 0);

			if ( norm)
			{				
				d.quad[i][0] *= width - 1;
				d.quad[i][1] *= height - 1;
			}
			
		}
	}
	
	// check if 3 corners are in one line they are linear

	for (int i = 0; i < 4; i++)
	{
		int a = (i + 1) % 4, b = (i + 2) % 4;
		if (d.quad[i][0] == d.quad[a][0] && d.quad[a][0] == d.quad[b][0]
			|| d.quad[i][1] == d.quad[a][1] && d.quad[a][1] == d.quad[b][1])

		{
			vsapi->setError(out, "reform: three x or y coord are equal and so not a quadrilateral");
			vsapi->freeNode(d.node[0]);
			vsapi->freeNode(d.node[1]);
			return;
		}
	}
	// check if it is convex
	for (int a = 0; a < 4; a++)
	{
		int b = (a + 1) % 4, c = (a + 2) % 4, dd = (a + 3) % 4;

		if (d.quad[a][0] < d.quad[c][0] && d.quad[b][0] > d.quad[c][0] && d.quad[dd][0] > d.quad[c][0]
			|| d.quad[a][1] < d.quad[c][1] && d.quad[b][1] > d.quad[c][1] && d.quad[dd][1] > d.quad[c][1])
		{
			vsapi->setError(out, "reform: x or y coords are resulting in a concave quad");
			vsapi->freeNode(d.node[0]);
			vsapi->freeNode(d.node[1]);
			return;
		}
	}

	d.intq = vsapi->propGetInt(in, "intq", 0, &err);
	if(err)
	{
		d.intq = 2;	//0 nearpt, 1 bilinear, 2 bicubic, 3 Lanczos 6x6 interpolation
	}
	else
	{
		if (d.intq < 0 || d.intq > 3)
		{
			vsapi->setError(out, "reform: invalid value for intq. 0 nearpt, 1 bilinear, 2 bicubic, 3 Lanczos 6x6 interpolation");
			vsapi->freeNode(d.node[0]);
			vsapi->freeNode(d.node[1]);
			return;
		}
	}

	temp = !!vsapi->propGetInt(in, "q2r", 0, &err);
	if (err)
		d.q2r = true;
	else
		d.q2r = temp == 0 ? false : true;
    
    data = (ReformData *) malloc(sizeof(d));

    *data = d;

    // Creates a new filter and returns a reference to it. Always pass on the in and out
    // arguments or unexpected things may happen. The name should be something that's
    
    vsapi->createFilter(in, out, "reform", reformInit, reformGetFrame, reformFree, fmParallel, 0, data, core);
}

//////////////////////////////////////////
// PluginInit
/*
VS_EXTERNAL_API(void) VapourSynthPluginInit(VSConfigPlugin configFunc, VSRegisterFunction registerFunc, VSPlugin *plugin) {
    configFunc("in.vcmohan.move", "reform", "VapourSynth Reform Plugin", VAPOURSYNTH_API_VERSION, 1, plugin);
    registerFunc("reform", "clip:clip;bkg:clip;intq:int:opt;norm:int:opt;rect:float[]:opt;quad:float[];q2r:int:opt;", reformCreate, 0, plugin);
}
*/