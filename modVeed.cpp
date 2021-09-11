
/******************************************************************************
veed filter plugin for vapoursynth by V.C.Mohan
 Removes Greenish noise or reddish noise. from images
 A gentle filter. Thread safe
 Author V.C.Mohan	
 20 Aug 2020

  Copyright (C) <2006 - 2020>  <V.C.Mohan>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    A copy of the GNU General Public License is at
    see <http://www.gnu.org/licenses/>.
	For details of how to contact author see <http://www.avisynth.org/vcmohan> 
*************************************************************************************************/  
typedef struct {
				VSNodeRef *node;
				const VSVideoInfo *vi;
				int str;				//strength of filter. 1 to 8
				int rad;				// pixels in radius of circle that influence
				int planes[3];			// planes to be processed
				int plimit[3];			// plus limit from ideal
				int mlimit[3];			// minus limit fro ideal
				float * kern;
				int ksize;

				}VeedData;
//-----------------------------------------------------------------------------------------------------
template <typename finc>
void VeedOut(finc * dp, const int dpitch,
				finc * wp, const int wpitch,
				const finc * sp, const int spitch,
				const int ht, const int wd, 
				const float *kr, int span, 
				finc low,  finc high);

template <typename finc>
void VeedOut(finc * dp, const int dpitch,
				finc * wp, const int wpitch,
				const finc * sp, const int spitch,
				const int ht, const int wd, 
				const float *kr, int span, 
				finc low,  finc high)
{

	int span2 = (span >> 1);
	
	int offset[17];
	for (int i = -span2; i < span2 + 1; i ++)
	{ 
		offset[i + span2] = i ;
	}
	// along horizontal
	const finc *spp = sp;
		finc *wpp = wp;
		finc *dpp = dp;

	for (int h = 0; h < ht; h ++)
	{
		
		// leave out left and right border pixels
		for (int w = span2; w < wd - span2; w ++)
		{
			float val = 0;

			for (int i = 0; i < span; i ++)
			{

				val += spp[w + offset[i]] * kr[i];
			}

			wpp[w] = val;			
		}
		spp += spitch;
		wpp += wpitch;
	}

	// vertical direction. Pre calculate offsets to save multiplication
	for (int i = -span2; i < span2 + 1; i ++)
	{
		offset[i + span2] = i * wpitch;
	}	

	for( int w = span2; w < wd - span2; w ++)
	{	
		
		finc *wpp = wp + span2 * wpitch;
		finc *dpp = dp + span2 * dpitch;
		const finc *spp = sp + span2 * spitch;
		// deveed along vertical
		 for (int h = span2; h < ht - span2; h ++)
		 {

			 float val = 0;

			 for (int i = 0; i < span; i ++)
			 {

				 val += wpp[offset[i]] * kr[i ];
			 }
			 
			 if ((spp[w] - val) > low)
			 {

				 dpp[w] = spp[w] - low;
			 }

			 else if ((val - spp[w]) > high)
			 {			

				 dpp[w] = spp[w] + high;
			 }

			 else
			 { 

				 dpp[w] = val;
			}
			 spp += spitch;
			 wpp += wpitch;
			 dpp += dpitch;
		 }		
	} 
}

//----------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------
// This function is called immediately after vsapi->createFilter(). This is the only place where the video
// properties may be set. In this case we simply use the same as the input clip. You may pass an array
// of VSVideoInfo if the filter has more than one output, like rgb+alpha as two separate clips.
static void VS_CC veedInit
				(VSMap *in, VSMap *out, void **instanceData, 
				VSNode *node, VSCore *core, const VSAPI *vsapi) 
{
    VeedData *d = (VeedData *) * instanceData;
    vsapi->setVideoInfo(d->vi, 1, node);	

	d->ksize = 2 * d->rad + 1;	

	d->kern =  vs_aligned_malloc<float>(sizeof(float)*d->ksize, 32);	
		
	for(int i = 0; i < d->ksize; i ++)
	{
		d->kern[i] = pow (2.71828, -(0.5 * (i - d->ksize / 2) * (i - d->ksize / 2) )/(d->str*d->str)) / (d->str * sqrt(6.2831853));
			
	
	}

}

//---------------------------------------------------------------------------------	

	// This is the main function that gets called when a frame should be produced. It will in most cases get
// called several times to produce one frame. This state is being kept track of by the value of
// activationReason. The first call to produce a certain frame n is always arInitial. In this state
// you should request all input frames you need. Always do it i ascending order to play nice with the
// upstream filters.
// Once all frames are ready the the filter will be called with arAllFramesReady. It is now time to
// do the actual processing.
static const VSFrameRef *VS_CC veedGetFrame
				(int n, int activationReason, void **instanceData, 
				void **frameData, VSFrameContext *frameCtx, 
				VSCore *core, const VSAPI *vsapi)
{
    VeedData *d = (VeedData *) * instanceData;

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
		VSFrameRef *work = vsapi->newVideoFrame(fi, width, height, src, core);

        // It's processing loop time!
        // Loop over all the planes
        
        for (int plane = 0; plane < fi->numPlanes; plane++)
		{
			
			const uint8_t *srcp = vsapi->getReadPtr(src, plane);
            int src_stride = vsapi->getStride(src, plane);
            uint8_t *dstp = vsapi->getWritePtr(dst, plane);
			int dst_stride = vsapi->getStride(dst, plane);
			uint8_t *workp = vsapi->getWritePtr(work, plane);
			int work_stride = vsapi->getStride(work, plane);
			int bht = vsapi->getFrameHeight(src, plane);            
            int bwd = vsapi->getFrameWidth(src, plane);
			int samplesize = fi->bytesPerSample;
			int spitch = src_stride / samplesize;
			int dpitch = dst_stride / samplesize;
			int wpitch = work_stride / samplesize;
			int nb = fi->bitsPerSample;
			
			if (d->planes[plane] == 0)
				continue;
			vs_bitblt(workp, work_stride, srcp, src_stride, bwd * samplesize, bht);

			if (fi->sampleType == stInteger)
			{
				if(fi->bitsPerSample == 8)
				{
					// 8 bit integer sample
					uint8_t low = ((1<< nb) * d->mlimit[plane])/ 200;
					uint8_t high = ((1<< nb) * d->plimit[plane])/ 200;
					
					VeedOut(dstp, dpitch, workp, wpitch, srcp, spitch,
					  bht, bwd, d->kern, d->ksize,  low, high);
				}

				else // if(fi->bitsPerSample > 8)
				{
					uint16_t low = ((1<< nb) * d->mlimit[plane])/ 200;
					uint16_t high = ((1<< nb) * d->plimit[plane])/ 200;
					uint16_t * dp = (uint16_t *) dstp;
					uint16_t * wp = (uint16_t *) workp;
					const uint16_t * sp = ( const uint16_t *) srcp;

					VeedOut(dp, dpitch, wp, wpitch, sp, spitch,
					  bht, bwd, d->kern, d->ksize, low, high);
				}

			} // integer

			else // if(fi->sampleType == stfloat)
			{
				float low = (d->mlimit[plane])/ 200;
				float high = (d->plimit[plane])/ 200;
				float * dp = (float *) dstp;
				float * wp = (float *) workp;
				const float * sp = ( const float *) srcp;

				VeedOut(dp, dpitch, wp, wpitch, sp, spitch,
					  bht, bwd, d->kern, d->ksize, low, high);
			}
		} // plane
		 // Release the source frame
       vsapi->freeFrame(src);
	   vsapi->freeFrame(work);
        // A reference is consumed when it is returned so saving the dst ref somewhere
        // and reusing it is not allowed.
		
       return dst;
    }

    return 0;
}
//-------------------------------------------------------------------
// Free all allocated data on filter destruction
static void VS_CC veedFree(void *instanceData, VSCore *core, const VSAPI *vsapi) 
{
   VeedData *d = (VeedData *)instanceData;
    vsapi->freeNode(d->node); 
	
	vs_aligned_free (d->kern);
	free(d);
}

// This function is responsible for validating arguments and creating a new filter
static void VS_CC veedCreate(const VSMap *in, 
					VSMap *out, void *userData, VSCore *core,
					const VSAPI *vsapi) 
{
    VeedData d;
    VeedData *data;
    int temp;
    int err;

    // Get a clip reference from the input arguments. This must be freed later.
    d.node = vsapi->propGetNode(in, "clip", 0, 0);
    d.vi = vsapi->getVideoInfo(d.node);

    // Note that
    // vi->format can be 0 if the input clip can change format midstream.
	if (d.vi->format->colorFamily == cmCompat) {
		vsapi->setError(out, "veed: Compat format not accepted.");
		vsapi->freeNode(d.node);
		return;
	}
	

    // If a property read fails for some reason (index out of bounds/wrong type)
    // then err will have flags set to indicate why and 0 will be returned. This
    // can be very useful to know when having optional arguments. Since we have
    // strict checking because of what we wrote in the argument string the only reason
    // this could fail is when the value wasn't set by the user.
    // And when it's not set we want it to default to enabled.
    
	d.rad = vsapi->propGetInt(in, "rad", 0, &err);

    if (err)
        d.rad = 5;

    // the only allowed 
    if (d.rad < 1 || d.rad > 8 ) 
	{
		vsapi->setError(out,"veed: rad value can be 1 to 8 only");
		vsapi->freeNode(d.node);
		return;
    }
	d.str = vsapi->propGetInt(in, "str", 0, &err);

    if (err)
        d.str = 5;

    // the only allowed 
    if (d.str < 1 || d.str > 8 ) 
	{
		vsapi->setError(out,"veed: str value can be 1 to 8 only");
		vsapi->freeNode(d.node);
		return;
    }

	if ( pow(2.71828, -(0.5 * (d.rad * d.rad) / (d.str*d.str)) / (d.str * sqrt(6.2831853))) < 2.0 / 255)
	{
		vsapi->setError(out,"veed: Either decrease rad or increase str to prevent wasteful processing.");
		vsapi->freeNode(d.node);
		return;
    }

	for ( int i = 0; i < 3; i ++)
	{
		d.planes[i] = !!vsapi->propGetInt(in, "planes", i, &err);

		if (err)
		 d.planes[i] = 1;

			// the only allowed 
		if (d.planes[i] < 0 || d.planes[i] > 1 ) 
		{
			vsapi->setError(out,"veed: planes value can be 0 or 1 only");
			vsapi->freeNode(d.node);
			return;
		}

		d.plimit[i] = vsapi->propGetInt(in, "plimit", i, &err);

		if (err)
		 d.plimit[i] = 3;

			// the only allowed 
		if (d.plimit[i] < 0 || d.plimit[i] > 10 ) 
		{
			vsapi->setError(out,"veed: plimit values can be 0 to 10 only");
			vsapi->freeNode(d.node);
			return;
		}
		d.mlimit[i] = vsapi->propGetInt(in, "mlimit", i, &err);

		if (err)
		 d.mlimit[i] = 3;

			// the only allowed 
		if (d.mlimit[i] < 0 || d.mlimit[i] > 10 ) 
		{
			vsapi->setError(out,"veed: mlimit values can be 0 to 10 only");
			vsapi->freeNode(d.node);
			return;
		}
	}

	if(d.planes[0] == 0 && d.planes[1] == 0 && d.planes[2] == 0 || (d.vi->format->colorFamily == cmGray && d.planes[0] == 0))
	{
			vsapi->setError(out,"veed: values of all planes are zero. At least one should be set to 1");
			vsapi->freeNode(d.node);
			return;
	}

	data =   (VeedData * ) malloc(sizeof(d));
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
    vsapi->createFilter(in, out, "veed", veedInit, 
								veedGetFrame, veedFree, 
								fmParallel, 0, data, core);
    return;
}

/*
VS_EXTERNAL_API(void) VapourSynthPluginInit(VSConfigPlugin configFunc, VSRegisterFunction registerFunc, VSPlugin *plugin)
{
    configFunc("vc.mohan.dns", "veed", "VapourSynth Veed", VAPOURSYNTH_API_VERSION, 1, plugin);
    registerFunc("Veed", "clip:clip;str:int:opt;rad:int:opt;planes:int[]:opt;plimit:int[]:opt;mlimit:int[]:opt;", veedCreate, 0, plugin);
}
*/



