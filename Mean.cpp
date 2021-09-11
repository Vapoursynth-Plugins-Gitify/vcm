/******************************************************************************
Mean filter plugin for vapoursynth by V.C.Mohan
Mean value is used if within tol .2 Nov 2020
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
	float tol;	// start %age value of greyness
	int grid;
	
}MeanData;
	
	template <typename finc>
	finc getAverage(finc* dp, int* offsets, int noff);

	template <typename finc>
	void blur(finc* dp, int* offsets, int noff, finc val);

	int setOffsets(int* offsets, int x, int y, int pitch);

/***************************************************************/
template <typename finc>
finc getAverage(finc* dp, int* offsets, int noff)
{
	float sum = 0;

	for (int i = 0; i < noff; i++)

		sum += (int)dp[offsets[i]];

	return (finc)(sum / noff);
}
template <typename finc>
void blur(finc* dp, int* offsets, int noff, finc val)
{
	for (int i = 0; i < noff; i++)
	{
		dp[offsets[i]] = val;
	}
}

int setOffsets(int* offsets, int x, int y, int pitch)
{
	int i = 0;

	for (int h = -y / 2; h <= y / 2; h++)
	{
		for (int w = -x / 2; w <= x / 2; w++)
		{
			if (h == 0 && w == 0) continue;
			offsets[i] = h * pitch + w ;
			i++;
		}
	}
	return i;
}
// This function is called immediately after vsapi->createFilter(). This is the only place where the video
// properties may be set. In this case we simply use the same as the input clip. You may pass an array
// of VSVideoInfo if the filter has more than one output, like rgb+alpha as two separate clips.
static void VS_CC meanInit(VSMap* in, VSMap* out, void** instanceData,
	VSNode* node, VSCore* core, const VSAPI* vsapi) 
{
	MeanData* d = (MeanData*)*instanceData;
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
	static const VSFrameRef* VS_CC meanGetFrame(int n, int activationReason, void** instanceData,
		void** frameData, VSFrameContext * frameCtx, VSCore * core, const VSAPI * vsapi)
	{
		MeanData* d = (MeanData*)*instanceData;

		if (activationReason == arInitial)
		{
			// Request the source frame on the first call
			vsapi->requestFrameFilter(n, d->node, frameCtx);
		}
		else if (activationReason == arAllFramesReady)
		{
			const VSFrameRef* src = vsapi->getFrameFilter(n, d->node, frameCtx);
			// The reason we query this on a per frame basis is because we want our filter
			// to accept clips with varying dimensions. If we reject such content using d->vi
			// would be better.
			const VSFormat* fi = d->vi->format;
			int height = vsapi->getFrameHeight(src, 0);
			int width = vsapi->getFrameWidth(src, 0);
			int nbytes = fi->bytesPerSample;
			int nbits = fi->bitsPerSample;
			VSFrameRef* dst = vsapi->copyFrame(src, core);


			int noff = d->grid * d->grid;	// max number of points
			int* offsets = (int*)vs_aligned_malloc <int>(sizeof(int) * d->grid * d->grid, 32);


			for (int plane = 0; plane < fi->numPlanes; plane++)
			{
				if (plane == 0 || fi->colorFamily == cmRGB
					|| (fi->colorFamily == cmYUV && fi->subSamplingH == 0 && fi->subSamplingW == 0)
					)
				{
					const uint8_t* srcp = vsapi->getReadPtr(src, plane);
					int src_stride = vsapi->getStride(src, plane);
					uint8_t* dstp = vsapi->getWritePtr(dst, plane);
					int stride = vsapi->getStride(dst, plane);
					int ht = vsapi->getFrameHeight(src, plane);
					int wd = vsapi->getFrameWidth(src, plane);

					int pitch = stride / nbytes;

					int noffs = setOffsets(offsets, d->grid, d->grid, pitch);

					srcp += d->grid / 2 * stride;
					dstp += d->grid / 2 * stride;

					for (int h = d->grid / 2; h < ht - d->grid / 2; h++)
					{
						for (int w = d->grid / 2; w < wd - d->grid / 2; w++)
						{
							if (fi->sampleType == stInteger)
							{
								if (nbytes == 1)
								{
									unsigned char mean, toler;

									mean = getAverage(srcp + w, offsets, noffs);
									toler = (unsigned char)(mean * d->tol);

									if (abs(mean - srcp[w]) < toler)
										dstp[w] = mean;
								}

								else if (nbytes == 2)
								{
									uint16_t mean, toler;

									mean = getAverage((uint16_t*)srcp + w, offsets, noffs);
									toler = (uint16_t)(mean * d->tol);

									if (abs(mean - *((uint16_t*)srcp + w)) < toler)
										*((uint16_t*)dstp + w) = mean;

								}

							}

							else	// float
							{
								float mean, toler;

								mean = getAverage((float*)srcp + w, offsets, noffs);
								toler = (float)(mean * d->tol);

								if (abs(mean - *((float*)srcp + w)) < toler)
									*((float*)dstp + w) = mean;

							}

						}
						srcp += stride;
						dstp += stride;
					}
				}
			}

			vs_aligned_free(offsets);
			vsapi->freeFrame(src);
			return dst;
		}
		return 0;
}

/***************************************************************/
// Free all allocated data on filter destruction
static void VS_CC meanFree(void* instanceData, VSCore* core, const VSAPI* vsapi)
{
	MeanData* d = (MeanData*)instanceData;
	vsapi->freeNode(d->node);	
	free(d);
}

// This function is responsible for validating arguments and creating a new filter
static void VS_CC meanCreate(const VSMap* in, VSMap* out, void* userData, VSCore* core, const VSAPI* vsapi) {
	MeanData d;
	MeanData* data;
	int err;

	// Get a clip reference from the input arguments. This must be freed later.
	d.node = vsapi->propGetNode(in, "clip", 0, 0);
	d.vi = vsapi->getVideoInfo(d.node);
	if (d.vi->format->colorFamily != cmRGB && d.vi->format->colorFamily != cmYUV && d.vi->format->colorFamily != cmGray)
	{
		vsapi->setError(out, "Mean: RGB, YUV and Gray color formats only for input allowed ");
		vsapi->freeNode(d.node);
		return;
	}
	if (d.vi->format->sampleType == stFloat && d.vi->format->bitsPerSample == 16)
	{
		vsapi->setError(out, "Mean: Half float formats not allowed ");
		vsapi->freeNode(d.node);
		return;
	}
	d.grid = int64ToIntS(vsapi->propGetInt(in, "grid", 0, &err));
	if (err)
	{
		d.grid = 5;
	}
	else
	{
		if (d.grid < 3 || d.grid > 11 || (d.grid % 2) == 0)
		{
			vsapi->setError(out, "Mean: value of grid need to be an odd number between 3 and 11");
			vsapi->freeNode(d.node);
			return;
		}
	}

	d.tol = (float)vsapi->propGetFloat(in, "tol", 0, &err);

	if (err)
	{
		d.tol = 0.05f;
	}
	else
	{
		if (d.tol < 0.01 || d.tol > 1.0)
		{
			vsapi->setError(out, "Mean: tol must have a value between 0.01 and 1.0");
			vsapi->freeNode(d.node);
			return;
		}
	}
	

	// I usually keep the filter data struct on the stack and don't allocate it
	// until all the input validation is done.
	data = (MeanData*)malloc(sizeof(d));
	*data = d;


	// If your filter is really fast (such as a filter that only resorts frames) you should set the
	// nfNoCache flag to make the caching work smoother.
	vsapi->createFilter(in, out, "Mean", meanInit, meanGetFrame, meanFree, fmParallel, 0, data, core);
}

//////////////////////////////////////////

/*
VS_EXTERNAL_API(void) VapourSynthPluginInit(VSConfigPlugin configFunc, VSRegisterFunction registerFunc, VSPlugin *plugin)
{
	configFunc("com.example.Mean", "vmod", "VapourSynth Mean ", VAPOURSYNTH_API_VERSION, 1, plugin);
	registerFunc("Mean", "clip:clip;grid:int:opt;tol:float:opt;", meanCreate, 0, plugin);
}
*/
