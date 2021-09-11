/******************************************************************************
modFan filter plugin for vapoursynth by V.C.Mohan
Removes diagonal (?) regular freq noise seen sometimes on capture from TV
specially if TV signal strength is low and interference can be any source like cell phone etc
White vertical one pixel wide streaks can  also be removed
  	
Last modified on 26 Aug 2017, 22 Aug 2020 

*************************************************************************************************/  

typedef struct {

					VSNodeRef *node;
					const VSVideoInfo *vi;
					bool edge;	//  pass, edge
					int span;		//  (wavelength  of noise)  pixels along x axis		
					bool yuv;		// Are y  u and v planes also to be processed?		
					float plustol;		// plus side tolerance for 		
					float minustol;		// minus side tolerance 		

				}FanData;
//------------------------------------------------------------------
//Definition of function
  
	template <typename finc>
	void FanFilterPlane( finc * op, int opitch, 			   
		int span, finc plus, finc minus,bool edge, int bwd, int bht, int kb);

//-----------------------------------------------------------------------------------

template <typename finc>
void FanFilterPlane(finc * op, int opitch, 
	int span, finc plus, finc minus,bool edge, int bwd, int bht, int kb)
{
	int span2 = span/2;

	float multsp2 = 1.0f / span2;

	op +=  span2 * kb;

	for (int h = 0; h < bht; h++)
	{

		int kbw = kb * span2;
		// span is odd number
		

		for (int w = span2; w < bwd - span2; w++)
		{		// get the first two sums
			
			// either side sums
			float sum1 = 0;
			float sum2 = 0;		
			float pval = *(op + kbw);	// central value
			int kbww = kb;

			for (int ww = 1; ww <= span2; ww++)
			{
				sum1 += op[kbw - kbww];

				sum2 += op[kbw + kbww];

				kbww += kb;
				
			}
					// get averages of both sides
			sum1 *= multsp2;
			sum2 *= multsp2;

			float avg = (sum1 + sum2) / 2.0f;

			if (sum1 == sum2)
			{
				
				op[kbw] = sum1;
				
			}

			else if ((pval < sum1 + plus &&  pval > sum1 - minus)
				|| (pval < sum2 + plus &&  pval > sum2 - minus))
			{
				// point value within range so no action
				// no action
			}
			else if (edge)
			{
			
				if ( (sum1 - pval) * (sum1 - pval) > (pval - sum2) * (pval - sum2))
					op[kbw] = sum2;	// near to sum2
				else
					op[kbw] = sum1;

			}

			else
			{
				// average value
				op[kbw] = avg;
			}
			kbw += kb;
		}

		op += opitch;
	}			
	
	
		
}
//---------------------------------------------------------------------------------
// This function is called immediately after vsapi->createFilter(). This is the only place where the video
// properties may be set. In this case we simply use the same as the input clip. You may pass an array
// of VSVideoInfo if the filter has more than one output, like rgb+alpha as two separate clips.
static void VS_CC fanInit
(VSMap *in, VSMap *out, void **instanceData,
VSNode *node, VSCore *core, const VSAPI *vsapi)
{
	FanData *d = (FanData *)* instanceData;
	vsapi->setVideoInfo(d->vi, 1, node);
}

//---------------------------------------------------------------------------------	

// This is the main function that gets called when a frame should be produced. It will in most cases get
// called several times to produce one frame. This state is being kept track of by the value of
// activationReason. The first call to produce a certain frame n is always arInitial. In this state
// you should request all input frames you need. Always do it i ascending order to play nice with the
// upstream filters.
// Once all frames are ready the the filter will be called with arAllFramesReady. It is now time to
// do the actual processing.
static const VSFrameRef *VS_CC fanGetFrame
(int n, int activationReason, void **instanceData,
void **frameData, VSFrameContext *frameCtx,
VSCore *core, const VSAPI *vsapi)
{
	FanData *d = (FanData *)* instanceData;

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


		// When creating a new frame for output it is VERY EXTREMELY SUPER IMPORTANT to
		// supply the "dominant" source frame to copy properties from. Frame props
		// are an essential part of the filter chain and you should NEVER break it.
		VSFrameRef * dst = vsapi->copyFrame(src, core);

		// It's processing loop time!
		// Loop over all the planes

		for (int plane = 0; plane < fi->numPlanes; plane++)
		{

			const uint8_t *srcp = vsapi->getReadPtr(src, plane);
			int src_stride = vsapi->getStride(src, plane);
			uint8_t *dp = vsapi->getWritePtr(dst, plane);
			int dst_stride = vsapi->getStride(dst, plane);
			int bht = vsapi->getFrameHeight(src, plane);
			int bwd = vsapi->getFrameWidth(src, plane);
			int samplesize = fi->bytesPerSample;
			int spitch = src_stride / samplesize;
			int dpitch = dst_stride / samplesize;
			int nb = fi->bitsPerSample;
			// copy input frame plane
		//	env->BitBlt(dp, dst_stride, srcp, src_stride, bwd * samplesize, bht);

			// do not process A
			if (plane == 3 || (plane > 0 && ! d->yuv && fi->colorFamily == cmYUV) ) continue;

			if (fi->sampleType == stInteger)
			{
				if (fi->bitsPerSample == 8)
				{
					unsigned char plus = d->plustol * 256, minus = d->minustol * 256;

					FanFilterPlane(dp , dpitch,

						d->span, plus, minus, d->edge, bwd, bht, 1);
				}

				else // if(fi->bitsPerSample > 8)
				{
					uint16_t plus = d->plustol * (1 << nb), minus = d->minustol * (1 << nb);

					FanFilterPlane((uint16_t*)dp, dpitch,

						d->span, plus, minus, d->edge, bwd, bht, 1);
				}
			}  // if stinteger
			else  // if(fi->sampleType == stfloat)
			{

				FanFilterPlane((float*)dp, dpitch,

					d->span, d->plustol, d->minustol, d->edge, bwd, bht, 1);

			}
		}
	
		vsapi->freeFrame(src);
		return dst;		
	}

return 0;
}

//-------------------------------------------------------------------
// Free all allocated data on filter destruction
static void VS_CC fanFree(void *instanceData, VSCore *core, const VSAPI *vsapi)
{
	FanData *d = (FanData *)instanceData;
	vsapi->freeNode(d->node);
	free(d);
}

///-----------------------------------------------------------------
// This function is responsible for validating arguments and creating a new filter
static void VS_CC fanCreate(const VSMap *in,
	VSMap *out, void *userData, VSCore *core,
	const VSAPI *vsapi)
{
	FanData d;
	FanData *data;
	int temp;
	int err;

	// Get a clip reference from the input arguments. This must be freed later.
	d.node = vsapi->propGetNode(in, "clip", 0, 0);
	d.vi = vsapi->getVideoInfo(d.node);

	// Note that
	// vi->format can be 0 if the input clip can change format midstream.
	if (d.vi->format->colorFamily == cmCompat) {
		vsapi->setError(out, "AdaptiveMedian: Compat format input not accepted.");
		vsapi->freeNode(d.node);
		return;
	}


	// If a property read fails for some reason (index out of bounds/wrong type)
	// then err will have flags set to indicate why and 0 will be returned. This
	// can be very useful to know when having optional arguments. Since we have
	// strict checking because of what we wrote in the argument string the only reason
	// this could fail is when the value wasn't set by the user.
	// And when it's not set we want it to default to enabled.

	d.span = vsapi->propGetInt(in, "span", 0, &err);

	if (err)
		d.span = 5;
	else if (d.span < 3 || d.span > 51 || (d.span & 1) == 0)
	{
		vsapi->setError(out, "fan: span value can be an odd number between 3 and 51 only");
		vsapi->freeNode(d.node);
		return;
	}

	temp = vsapi->propGetInt(in, "edge", 0, &err);

	if (err)
		d.edge = true;
	else if (temp < 0 || temp > 1)
	{
		vsapi->setError(out, "fan: edge value can be 0 or 1 only");
		vsapi->freeNode(d.node);
		return;
	}
	else
	{
		d.edge = temp == 0 ? false : true;
	}

	temp = vsapi->propGetInt(in, "uv", 0, &err);

	if (err)
		d.yuv = true;
	else if (temp < 0 || temp > 1)
	{
		vsapi->setError(out, "fan: uv value can be 0 or 1 only");
		vsapi->freeNode(d.node);
		return;
	}
	else
	{
		d.yuv = temp == 0 ? false : true;
	}

	d.plustol = vsapi->propGetFloat(in, "plus", 0, &err);
	if (err)
		d.plustol = 0.02f;
	else if (d.plustol < 0 ||d.plustol > 0.5f)
	{
		vsapi->setError(out, "fan: plustol can have a value between 0.0 and 0.5 only");
		vsapi->freeNode(d.node);
		return;
	}
	d.minustol = vsapi->propGetFloat(in, "minus", 0, &err);
	if (err)
		d.minustol = 0.02f;
	else if (d.minustol < 0 || d.minustol > 0.5f)
	{
		vsapi->setError(out, "fan: minustol can have a value between 0.0 and 0.5 only");
		vsapi->freeNode(d.node);
		return;
	}
	
	data = (FanData *)malloc(sizeof(d));
	*data = d;

	// Create a new filter and returns a reference to it. Always pass on the in and out
	// arguments or unexpected things may happen. The name should be something that's
	// easy to connect to the filter, like its function name.
	// The three function pointers handle initialization, frame processing and filter destruction.
	// The filtermode is very important to get right as it controls how threading of the filter
	// is handled. In general you should only use fmParallel whenever possible. This is if you
	// need to modify no shared data at all when the filter is running.
	// For more complicated filters fmParallelRequests is usually easier to achieve as an
	// be prefetched in parallel but the actual processing is serialized.
	// The others can be considered special cases where fmSerial is useful to source filters and
	// fmUnordered is useful when a filter's state may change even when deciding which frames to
	// prefetch (such as a cache filter).
	// If you filter is really fast (such as a filter that only resorts frames) you should set the
	// nfNoCache flag to make the caching work smoother.
	vsapi->createFilter(in, out, "fan", fanInit,
		fanGetFrame, fanFree,
		fmParallel, 0, data, core);
	return;
}

/*
VS_EXTERNAL_API(void) VapourSynthPluginInit(VSConfigPlugin configFunc, VSRegisterFunction registerFunc, VSPlugin *plugin)
{
configFunc("vc.mohan.dns", "fan", "VapourSynth Fan", VAPOURSYNTH_API_VERSION, 1, plugin);
registerFunc("Fan", "clip:clip;span:int:opt;edge:int:opt;plus:int:opt;minus:int:opt;uv:int:opt;", fanCreate, 0, plugin);
}
*/

