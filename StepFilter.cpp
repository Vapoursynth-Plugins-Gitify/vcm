/******************************************************************************
StepFilter filter plugin for  avisynth+ by V.C.Mohan
 Calculates frame and horizontal line averages and corrects
 additive or multiplicative regular noise

 Thread safe (thread agnostic)

 Author V.C.Mohan
20  April 2021

  Copyright (C) <2021>  <V.C.Mohan>

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, version 3 of the License.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	A copy of the GNU General Public License is at
	see <http://www.gnu.org/licenses/>.
	For details of how to contact author see <http://www.avisynth.nl/users/vcmohan/vcmohan.html>

*************************************************************************************************/  
/*
#include <stdlib.h>
#include "VapourSynth.h"
#include "VSHelper.h"
#include "math.h"
*/
typedef struct {
	VSNodeRef* node;
	const VSVideoInfo* vi;
	float boost; // change reference amplitude by this
	bool add;	// correction is additive or multiplacative
	int segmenthor;	// horizontal correction segment size
	int segmentvert;	// vertical correction segment size
	bool limit; // Do not correct if correction is -ve
	
}StepFilterData;
	
template <typename finc>
float getFrameAverage(const finc* sp, const int pitch,
	float  boost, const int ht, const int wd);
template <typename finc>
void getLineAverageAndCorrect(finc* dp,  int wd,
	float corr, finc min, finc max, bool add, bool limit);
template <typename finc>
float getLineAverage(const finc* sp, int wd);
template <typename finc>
void correctLineAmplitude(finc* sp, int wd,
	float corr, finc min, finc max, bool add, bool limit);
//template <typename finc>
//finc clampVal(float val, finc min, finc max);

template <typename finc>
void getSegLineAverageAndCorrect(finc* dp, const finc* sp,  int wd,
	float corr, finc min, finc max, bool add, bool limit);

template <typename finc>
void fullProcessStepFilter(finc* dp, const finc* sp, int pitch,float boost,
	int wd, int ht, finc min, finc max, bool add, bool limit,
	int segmenthor, int segmentvert);

/***************************************************************/
template <typename finc>
float getLineAverage(const finc* sp, int wd)
{
	float sum = 0.0;

	for (int w = 0; w < wd; w++)
	{
		sum += *(sp + w );

	}
	return sum / wd;
}

template <typename finc>
void correctLineAmplitude(finc* dp, int wd, 
	float corr, finc min, finc max, bool add, bool limit)
{
	// limit is used from class
	if (limit && corr < 0)
		return;
	if (add)
	{

		for (int w = 0; w < wd; w++)
			*(dp + w )
			= clamp(*(dp + w ) + corr, min, max);
	}
	else
	{
		for (int w = 0; w < wd; w++)
			*(dp + w )
			= clamp(*(dp + w ) * corr, min, max);
	}
}
/*template <typename finc>
finc StepFilter::clampVal(float val, finc min, finc max)
{
	return (finc)(val < min ? min : val > max ? max : val);
}*/


//-----------------------------------------------------------------------
template <typename finc>
float getFrameAverage(const finc* sp, const int pitch,
	float  boost, const int ht, const int wd)
{
	float corr = 0.0;

	for (int h = 0; h < ht; h++)
	{
		corr += getLineAverage(sp, wd);

		sp += pitch;
	}

	//stepfilter = sum / ht;
	return boost * corr / ht;
}
/***************************************************************/
template <typename finc>
void getLineAverageAndCorrect(finc* dp,  int wd,
	float corr, finc min, finc max, bool add, bool limit)
{
	float sum = getLineAverage(dp, wd);
	float lineCorr = add ? corr - sum : corr / sum; // > 0 ? corr - sum : 0;
	//float lineCorr = corr / sum;

	correctLineAmplitude(dp, wd, 
		lineCorr, min, max, add, limit);
}

template <typename finc>
void getSegLineAverageAndCorrect(finc* dp, const finc* sp, int wd,
	float corr, finc min, finc max, bool add, bool limit)
{
	float sum = getLineAverage(sp, wd);
	float lineCorr = add ? corr - sum : corr / sum;

	correctLineAmplitude(dp + wd / 2 , 1,
		lineCorr, min, max, add, limit);
}

template <typename finc>
void fullProcessStepFilter( finc* dp, const finc* sp, int pitch, float boost,
	int wd, int ht, finc min, finc max, bool add, bool limit, int segmenthor, int segmentvert)
{
	float frameCorr;

	if (segmentvert == 0)
	{

		frameCorr = getFrameAverage(sp, pitch,
			boost, ht, wd);

		for (int h = 0; h < ht; h++)
		{
			if (segmenthor == 0)
			{
				getLineAverageAndCorrect(dp,
					wd, frameCorr, min, max, add, limit);
			}
			else // if (segmenthor != 0)
			{
				// extend left
				float segCorr = getLineAverage(sp, segmenthor);
				float factor = add ? frameCorr * boost - segCorr : frameCorr * boost / segCorr;

				correctLineAmplitude(dp, segmenthor / 2, 
					factor, min, max, add, limit);

				for (int w = 0; w < wd - segmenthor; w++)
				{
					getSegLineAverageAndCorrect(dp + w , sp + w ,
						segmenthor, frameCorr, min, max, add, limit);
				}
				// extend right
				segCorr = getLineAverage(sp + (wd - segmenthor), segmenthor);

				factor = add ? frameCorr * boost - segCorr : frameCorr * boost / segCorr;

				correctLineAmplitude(dp + (wd - segmenthor / 2) , segmenthor / 2, 
					factor, min, max, add, limit);
			}
			dp += pitch;
			sp += pitch;
		}
	}
	else if (segmentvert != 0)
	{
		for (int h = 0; h < ht - segmentvert; h++)
		{
			frameCorr = getFrameAverage(sp + h * pitch, pitch,
				boost, segmentvert, wd);

			if (segmenthor == 0)
			{
				getLineAverageAndCorrect(dp + (h + segmentvert / 2) * pitch,
					 wd, frameCorr, min, max, add, limit);
			}
			else
			{
				// extend to left margin
				float segCorr = getLineAverage(sp + (h + segmentvert / 2) * pitch,segmenthor);
				float factor = add ? frameCorr * boost - segCorr 
					: frameCorr * boost / segCorr;

				correctLineAmplitude(dp + (h + segmentvert / 2) * pitch, segmenthor / 2,
					factor, min, max, add, limit);
				// process 
				for (int w = 0; w < wd - segmenthor; w++)
					getSegLineAverageAndCorrect(
						dp + (h + segmentvert / 2) * pitch + w,
						sp + (h + segmentvert / 2) * pitch + w,
						segmenthor, frameCorr, min, max, add, limit);

				// extend right margin
				segCorr = getLineAverage(sp + (h + segmentvert / 2) * pitch + (wd - segmenthor), segmenthor);
				factor = add ? frameCorr * boost - segCorr : frameCorr * boost / segCorr;
				correctLineAmplitude(dp + (h + segmentvert / 2) * pitch + (wd - segmenthor / 2), segmenthor / 2,
					factor, min, max, add, limit);
			}
		}
	}
}
// This function is called immediately after vsapi->createFilter(). This is the only place where the video
// properties may be set. In this case we simply use the same as the input limit. You may pass an array
// of VSVideoInfo if the filter has more than one output, like rgb+alpha as two separate limits.
static void VS_CC stepfilterInit(VSMap* in, VSMap* out, void** instanceData,
	VSNode* node, VSCore* core, const VSAPI* vsapi) 
{
	StepFilterData* d = (StepFilterData*)*instanceData;
	vsapi->setVideoInfo(d->vi, 1, node);

}

	//----------------------------------------------------------------------------------------------------------

	// This is the main function that gets called when a frame should be produced. It will, in most cases, get
	// called several times to produce one frame. This state is being kept track of by the value of
	// activationReason. The first call to produce a certain frame n is always arInitial. In this state
	// you should request all the input frames you need. Always do it in ascending order to play nice with the
	// upstream filters.
	// Once all frames are ready, the filter will be called with arAllFramesReady. It is now time to
	// do the actual processing.
static const VSFrameRef* VS_CC stepfilterGetFrame(int n, int activationReason, void** instanceData,
	void** frameData, VSFrameContext* frameCtx, VSCore* core, const VSAPI* vsapi)
{
	StepFilterData* d = (StepFilterData*)*instanceData;

	if (activationReason == arInitial)
	{
		// Request the source frame on the first call
		vsapi->requestFrameFilter(n, d->node, frameCtx);
	}
	else if (activationReason == arAllFramesReady)
	{
		const VSFrameRef* src = vsapi->getFrameFilter(n, d->node, frameCtx);
		// The reason we query this on a per frame basis is because we want our filter
		// to accept limits with varying dimensions. If we reject such content using d->vi
		// would be better.
		const VSFormat* fi = d->vi->format;
		int height = vsapi->getFrameHeight(src, 0);
		int width = vsapi->getFrameWidth(src, 0);
		int nbytes = fi->bytesPerSample;
		int nbits = fi->bitsPerSample;
		VSFrameRef* dst = vsapi->copyFrame(src, core);
		
		int np = fi->colorFamily == cmRGB ? 3 : 1;

		for (int plane = 0; plane < np; plane++)
		{
			const uint8_t* sp = vsapi->getReadPtr(src, plane);
			//int src_stride = vsapi->getStride(src, plane);
			uint8_t* dp = vsapi->getWritePtr(dst, plane);
			int stride = vsapi->getStride(dst, plane);
			int ht = vsapi->getFrameHeight(src, plane);
			int wd = vsapi->getFrameWidth(src, plane);

			int pitch = stride / nbytes;


			if (fi->sampleType == stInteger)
			{
				if (nbytes == 1)
				{
					unsigned char min = 0, max = 255;

					if (fi->colorFamily == cmYUV)
						min = 16, max = 235;

					fullProcessStepFilter(dp,  sp, pitch, d->boost,
										wd, ht, min, max, d->add, d->limit,
										d->segmenthor, d->segmentvert);
				}

				else if (nbytes == 2)
				{
					uint16_t min = 0, max = 255 << (nbits - 8);

					if (fi->colorFamily == cmYUV)
						min = 16 << (nbits - 8), max = 235 << (nbits - 8);

					fullProcessStepFilter((uint16_t*)dp, (const uint16_t*)sp, pitch, d->boost,
											wd, ht, min, max, d->add, d->limit,
											d->segmenthor, d->segmentvert);
				}
			}

			else	// float
			{
				float min = 0, max = 1.0;

				fullProcessStepFilter((float*)dp, (const float*)sp, pitch, d->boost,
										wd, ht, min, max, d->add, d->limit,
										d->segmenthor, d->segmentvert);
			}
		}
			
		vsapi->freeFrame(src);
		return dst;
	}
	return 0;
}

/***************************************************************/
// Free all allocated data on filter destruction
static void VS_CC stepfilterFree(void* instanceData, VSCore* core, const VSAPI* vsapi)
{
	StepFilterData* d = (StepFilterData*)instanceData;
	vsapi->freeNode(d->node);	
	free(d);
}

// This function is responsible for validating arguments and creating a new filter
static void VS_CC stepfilterCreate(const VSMap* in, VSMap* out, void* userData, VSCore* core, const VSAPI* vsapi) {
	StepFilterData d;
	StepFilterData* data;
	int err;

	// Get a limit reference from the input arguments. This must be freed later.
	d.node = vsapi->propGetNode(in, "clip", 0, 0);
	d.vi = vsapi->getVideoInfo(d.node);
	if (d.vi->format->colorFamily != cmRGB && d.vi->format->colorFamily != cmYUV && d.vi->format->colorFamily != cmGray)
	{
		vsapi->setError(out, "StepFilter: RGB, YUV and Gray color formats only for input allowed ");
		vsapi->freeNode(d.node);
		return;
	}
	if (d.vi->format->sampleType == stFloat && d.vi->format->bitsPerSample == 16)
	{
		vsapi->setError(out, "StepFilter: Half float formats not allowed ");
		vsapi->freeNode(d.node);
		return;
	}
	int temp = !!int64ToIntS(vsapi->propGetInt(in, "add", 0, &err));
	if (err)
	{
		d.add = true;
	}
	else if (temp == 0)
	{
		d.add = false;
	}
	else
		d.add = true;

	temp = !!int64ToIntS(vsapi->propGetInt(in, "limit", 0, &err));
	if (err)
	{
		d.limit = false;
	}
	else if (temp == 0)
	{
		d.limit = false;
	}
	else
		d.limit = true;

	d.boost = vsapi->propGetFloat(in, "boost", 0, &err);

	if (err)
	{
		d.boost = 1.0f;
	}
	else
	{
		if (d.boost < 0.5 || d.boost > 5.0)
		{
			vsapi->setError(out, "StepFilter: boost must have a value between 0.5 and 5.0");
			vsapi->freeNode(d.node);
			return;
		}
	}
	d.segmenthor = int64ToIntS(vsapi->propGetInt(in, "segmenthor", 0, &err));
	if (err)
		d.segmenthor = 60;
	else if ( d.segmenthor != 0 && d.segmenthor < 16 && d.segmenthor >= d.vi->width / 2)
	{
		vsapi->setError(out, "StepFilter: segmenthor must be either zero or have a value between 16 and frame width / 2");
		vsapi->freeNode(d.node);
		return;
	}
	d.segmentvert = int64ToIntS(vsapi->propGetInt(in, "segmentvert", 0, &err));
	if (err)
		d.segmentvert = 60;
	else if (d.segmentvert != 0 && d.segmentvert < 16 && d.segmentvert >= d.vi->height / 2)
	{
		vsapi->setError(out, "StepFilter: segmentvert must be either zero or have a value between 16 and frame height / 2");
		vsapi->freeNode(d.node);
		return;
	}

	// I usually keep the filter data struct on the stack and don't allocate it
	// until all the input validation is done.
	data = (StepFilterData*)malloc(sizeof(d));
	*data = d;


	// If your filter is really fast (such as a filter that only resorts frames) you should set the
	// nfNoCache flag to make the caching work smoother.
	vsapi->createFilter(in, out, "StepFilter", stepfilterInit, stepfilterGetFrame, stepfilterFree, fmParallel, 0, data, core);
}

//////////////////////////////////////////

/*
VS_EXTERNAL_API(void) VapourSynthPluginInit(VSConfigPlugin configFunc, VSRegisterFunction registerFunc, VSPlugin *plugin)
{
	configFunc("com.example.StepFilter", "vmod", "VapourSynth StepFilter ", VAPOURSYNTH_API_VERSION, 1, plugin);
	registerFunc("StepFilter", "clip:clip;add:int:opt;boost:float:opt;"
	"segmenthor:int:opt;segmentvert:int:opt;limit:int:opt;", stepfilterCreate, 0, plugin);
}
*/
