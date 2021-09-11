/******************************************************************************
ColorBox filter plugin for vapoursynth by V.C.Mohan
 Creates a color boxes that can change along clip length.
 Author V.C.Mohan
 10 Oct 2020

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
typedef struct {
	
	VSFrameRef* f;
	VSVideoInfo vi;
	bool keep;
	int luma;

	int hIncU, wIncU;		// color inc per pixel in strip
	int hBoxIncU, wBoxIncU; // color increment per box along width and height

	int hIncV, wIncV;		// color inc per pixel in strip
	int hBoxIncV, wBoxIncV; // color increment per box along width and height	
	int wInc, hInc;			// color increment within a box
	int nBoxesW;			// number of boxes along width
	int nBoxesH;			// number of boxes along height
	
	

} ColorBoxData;
//------------------------------------------------------------------------------------------------------

template <typename finc>
void paintBoxes(finc* dp, int pitch, finc col, finc winc, finc hinc, int wd, int ht, int nbits);
//------------------------------------------------------------------------------------------------------

// This function is called immediately after vsapi->createFilter(). This is the only place where the video
// properties may be set. In this case we simply use the same as the input clip. You may pass an array
// of VSVideoInfo if the filter has more than one output, like rgb+alpha as two separate clips.
static void VS_CC colorBoxInit(VSMap* in, VSMap* out, void** instanceData, VSNode* node, VSCore* core, const VSAPI* vsapi) {
	ColorBoxData* d = (ColorBoxData*)*instanceData;
	//vsapi->setVideoInfo(d->vi, 1, node);
	
	vsapi->setVideoInfo(&d->vi, 1, node);
	d->f = NULL;

}
//---------------------------------------------------------------------------------------------------------

// create color box for a plane
template <typename finc>
void paintBoxes(finc* dp, int pitch, finc col, finc winc, finc hinc, int wd, int ht, int nbits)
{
	
	if (nbits == 8)
	{
		uint8_t  andColor = (uint8_t)((1 << nbits) - 1);

		uint8_t wcolNew;
		for (int h = 0; h < ht; h++)
		{
			wcolNew = col;

			for (int w = 0; w < wd; w++)
			{
				dp[w] = wcolNew;
				wcolNew += winc;
				wcolNew &= andColor;
			}

			dp += pitch;
			col += hinc;
		}
	}
	else if (nbits <= 16)
	{
		uint16_t  andColor = (uint16_t)((1 << nbits) - 1);

		uint16_t wcolNew;
		for (int h = 0; h < ht; h++)
		{
			wcolNew = col;

			for (int w = 0; w < wd; w++)
			{
				dp[w] = wcolNew;
				wcolNew += winc;
				wcolNew &= andColor;
			}

			dp += pitch;
			col += hinc;
		}
	}
	else // float 
	{
		uint8_t col8 = (uint8_t)(255 * col);
		uint8_t winc8 = (uint8_t)(255 * winc);
		uint8_t hinc8 = (uint8_t)(255 * hinc);
		uint8_t andColor = (uint8_t)(255);
		uint8_t wcolNew = 0;

		for (int h = 0; h < ht; h++)
		{
			wcolNew = col8;

			for (int w = 0; w < wd; w++)
			{
				dp[w] = (float)(wcolNew - 128) / 255.0;
				wcolNew += winc8;
				wcolNew &= andColor;
			}

			dp += pitch;
			col8 += hinc8;
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
static const VSFrameRef* VS_CC colorBoxGetFrame(int n, int activationReason, void** instanceData,
	void** frameData, VSFrameContext* frameCtx, VSCore* core, const VSAPI* vsapi)
{
	ColorBoxData* d = (ColorBoxData*)*instanceData;

	if (activationReason == arInitial)
	{
		VSFrameRef* dst = NULL;
		if (d->f == NULL) 
		{
			dst = vsapi->newVideoFrame(d->vi.format, d->vi.width, d->vi.height, 0, core);
			const VSFormat* fi = d->vi.format;
			int height = vsapi->getFrameHeight(dst, 0);
			int width = vsapi->getFrameWidth(dst, 0);
			int nbytes = fi->bytesPerSample;
			int nbits = fi->bitsPerSample;
			int pitch = vsapi->getStride(dst, 0) / nbytes;
			uint8_t* dptr = vsapi->getWritePtr(dst, 0);
			float brightness = d->luma / 100.0f;
			int subH = fi->subSamplingH, subW = fi->subSamplingW;
			// Will only color Luma plane
			
			for (int nh = 0; nh < d->nBoxesH; nh++)
			{

				int nht = nh == (d->nBoxesH - 1) ? ((height / d->nBoxesH) + (height % d->nBoxesH)) : (height / d->nBoxesH);
				uint8_t* dpw = dptr;

				for (int nw = 0; nw < d->nBoxesW; nw++)
				{
					int nwd = nw == d->nBoxesW - 1 ? (width / d->nBoxesW) + (width % d->nBoxesW) : (width / d->nBoxesW);

					float bright = brightness -((float)(rand() % (d->luma / 4)) / 100.0f);

					switch (nbytes)
					{
					case 1:
						fillPlaneWithVal(dpw , pitch, nwd, nht, (uint8_t)(bright * (1 << nbits)));
						break;
					case 2:
						fillPlaneWithVal((uint16_t*)dpw, pitch, nwd, nht, (uint16_t)(bright * (1 << nbits)));
						break;
					case 4:
						fillPlaneWithVal((float*)dpw, pitch, nwd, nht, bright);
						break;
					}
					dpw += nwd * nbytes;
					
				}

				dptr += nht * pitch * nbytes;
			}

			for (int plane = 1; plane < fi->numPlanes; plane++)
			{
				uint8_t* dstp = vsapi->getWritePtr(dst, plane);
				int dst_stride = vsapi->getStride(dst, plane);
				int ht = vsapi->getFrameHeight(dst, plane);
				int wd = vsapi->getFrameWidth(dst, plane);
				int pitch = dst_stride / nbytes;

				int boxht = ht / d->nBoxesH, boxwd = wd / d->nBoxesW, hBox, wBox;

				if (nbytes == 1)
				{
					uint8_t col = (uint8_t)(96), winc, hinc;
					uint8_t wBoxInc = plane == 1 ? (uint8_t)d->wBoxIncU : (uint8_t)d->wBoxIncV;
					uint8_t hBoxInc = plane == 1 ? (uint8_t)d->hBoxIncU : (uint8_t)d->hBoxIncV;
					
					uint8_t win = (uint8_t)(d->wInc);
					uint8_t hin = (uint8_t)(d->hInc);
					uint8_t* dp = dstp;

					for (int nh = 0; nh < d->nBoxesH; nh++)
					{
						int nht = nh == d->nBoxesH - 1 ? boxht + (ht % d->nBoxesH) : boxht;
						int w = 0;
						if (nh == d->nBoxesH / 2)
							hinc = hin;
						else
							hinc = 0;

						for (int nw = 0; nw < d->nBoxesW; nw++)
						{
							int nwd = nw == d->nBoxesW - 1 ? boxwd + (wd % d->nBoxesW) : boxwd;

							if (nw == d->nBoxesW / 2)
								winc = win;
							else
								winc = 0;
							paintBoxes(dp + w, pitch, col, winc, hinc, nwd, nht, nbits);

							col += wBoxInc;

							w += boxwd;
						}

						col += hBoxInc;

						dp += pitch * boxht;
					}
				}

				else  if (nbytes == 2)
				{

					uint16_t* dp = (uint16_t*)dstp;
					uint16_t col = (uint16_t)(96 << (nbits - 8)), winc = 0, hinc = 0;
					uint16_t wBoxInc = plane == 1 ? (uint16_t)(d->wBoxIncU << (nbits - 8)) : (uint16_t)(d->wBoxIncV << (nbits - 8));
					uint16_t hBoxInc = plane == 1 ? (uint16_t)(d->hBoxIncU << (nbits - 8)) : (uint16_t)(d->hBoxIncV << (nbits - 8));
					uint16_t  win = (uint16_t)(d->wInc << (nbits - 8));
					uint16_t hin = (uint16_t)(d->hInc << (nbits - 8));

					for (int nh = 0; nh < d->nBoxesH; nh++)
					{
						int nht = nh == d->nBoxesH - 1 ? boxht + (ht % d->nBoxesH) : boxht;
						int w = 0;
						if (nh == d->nBoxesH / 2)
							hinc = hin;
						else
							hinc = 0;

						for (int nw = 0; nw < d->nBoxesW; nw++)
						{
							int nwd = nw == d->nBoxesW - 1 ? boxwd + (wd % d->nBoxesW) : boxwd;

							if (nw == d->nBoxesW / 2)
								winc = win;
							else
								winc = 0;
							paintBoxes(dp + w, pitch, col, winc, hinc, nwd, nht, nbits);

							col += wBoxInc;

							w += boxwd;
						}

						col += hBoxInc;

						dp += pitch * boxht;
					}
				}
				else // float
				{
					float* dp = (float*)dstp;
					
					float col = (float)(96 / 255.0f);
					float win = (float)(d->wInc) / 255.0f;
					float hin = (float)(d->hInc) / 255.0f;
					float wBoxInc = plane == 1 ? (float)(d->wBoxIncU / 255.0f) : (float)(d->wBoxIncV / 255.0f);
					float hBoxInc = plane == 1 ? (float)(d->hBoxIncU / 255.0f) : (float)(d->hBoxIncV / 255.0f);
					float hinc, winc;

					for (int nh = 0; nh < d->nBoxesH; nh++)
					{
						int nht = nh == d->nBoxesH - 1 ? boxht + (ht % d->nBoxesH) : boxht;
						int w = 0;
						if (nh == d->nBoxesH / 2)
							hinc = hin;
						else
							hinc = 0;

						for (int nw = 0; nw < d->nBoxesW; nw++)
						{
							int nwd = nw == d->nBoxesW - 1 ? boxwd + (wd % d->nBoxesW) : boxwd;

							if (nw == d->nBoxesW / 2)
								winc = win;
							else
								winc = 0;

							paintBoxes(dp + w, pitch, col, winc, hinc, nwd, nht, nbits);

							col += wBoxInc;
							w += boxwd;
						}

						col += hBoxInc;

						dp += pitch * boxht;
					}
				}
			}
		

			if (d->vi.fpsNum > 0)
			{
				VSMap* frameProps = vsapi->getFramePropsRW(dst);
				vsapi->propSetInt(frameProps, "_DurationNum", d->vi.fpsDen, paReplace);
				vsapi->propSetInt(frameProps, "_DurationDen", d->vi.fpsNum, paReplace);
			}
		}

		if (d->keep)
		{
			if (dst != NULL)
				d->f = dst;
			return vsapi->cloneFrameRef(d->f);
		}
		else 
		{
			return dst;
		}
	}
	return 0;
}

	
		
// Free all allocated data on filter destruction
static void VS_CC colorBoxFree(void* instanceData, VSCore* core, const VSAPI* vsapi)
{
	ColorBoxData* d = (ColorBoxData*)instanceData;
	vsapi->freeFrame(d->f);
	
	free(d);
}

// This function is responsible for validating arguments and creating a new filter
static void VS_CC colorBoxCreate(const VSMap* in, VSMap* out, void* userData, VSCore* core, const VSAPI* vsapi) {
	ColorBoxData d;
	ColorBoxData* data;

	int err;
	d.vi.width = 640;
	d.vi.height = 480;
	d.vi.fpsNum = 24;
	d.vi.fpsDen = 1;
	
	d.vi.numFrames = int64ToIntS((d.vi.fpsNum * 100) / d.vi.fpsDen);
	
	d.luma = (int)vsapi->propGetInt(in, "luma", 0, &err);
	
	if (err)
		d.luma = 40;
	else
	{
		if (d.luma > 99 || d.luma < 1)
		{
			vsapi->setError(out, "ColorBox: Luma percentage must be between 1 & 99");
			
			return;
		}
	}
	d.nBoxesW = (int)vsapi->propGetInt(in, "nbw", 0, &err);

	if (err)
		d.nBoxesW = 6;
	else
	{
		if (d.nBoxesW > 12 || d.nBoxesW < 2)
		{
			vsapi->setError(out, "ColorBox: nbw must be between 2 and 12");

			return;
		}
	}

	d.nBoxesH = (int)vsapi->propGetInt(in, "nbh", 0, &err);

	if (err)
		d.nBoxesH = 4;
	else
	{
		if (d.nBoxesH > 12 || d.nBoxesH < 2)
		{
			vsapi->setError(out, "ColorBox: nbh must be between 2 and 12");

			return;
		}
	}
	int format = int64ToIntS(vsapi->propGetInt(in, "format", 0, &err));

	if (err)
		d.vi.format = vsapi->getFormatPreset(pfYUV444P8, core);
	else
	{
		d.vi.format = vsapi->getFormatPreset(format, core);
		if ( !d.vi.format)
		{
			vsapi->setError(out, "ColorBox: invalid format");
			
			return;
		}
	}
	if(d.vi.format->colorFamily != cmYUV || (d.vi.format->id == pfYUV444PH) )
	{
		vsapi->setError(out, "ColorBox: YUV integer and single float formats only allowed");
		
		return;
	}
	int subH = d.vi.format->subSamplingH, subW = d.vi.format->subSamplingW;

	
	d.wBoxIncU = 29 << subW;
	d.hBoxIncU = 117 << subH;
	d.wIncU = 64 << subW;
	d.hIncU = 64 << subH;
	d.wBoxIncV = 79 << subW;
	d.hBoxIncV = 47 << subH;
	d.wIncV = 79 << subW;
	d.hIncV = 37 << subH;
	d.wInc = 4 << subW;
	d.hInc = 4 << subH;
	d.keep = true;

	// I usually keep the filter data struct on the stack and don't allocate it
	// until all the input validation is done.
	data = (ColorBoxData*)malloc(sizeof(d));
	*data = d;


	// If your filter is really fast (such as a filter that only resorts frames) you should set the
	// nfNoCache flag to make the caching work smoother.
	vsapi->createFilter(in, out, "colorBox", colorBoxInit, colorBoxGetFrame, colorBoxFree, fmUnordered, nfNoCache, data, core);
}

//////////////////////////////////////////

/*
VS_EXTERNAL_API(void) VapourSynthPluginInit(VSConfigPlugin configFunc, VSRegisterFunction registerFunc, VSPlugin *plugin)
{
	configFunc("com.example.ColorBox", "vcm", "VapourSynth ColorBox", VAPOURSYNTH_API_VERSION, 1, plugin);
	registerFunc("ColorBox", "format:int:opt;luma:int:opt;nbw:int:opt;nbh:int:opt;", colorBoxCreate, 0, plugin);
}
*/