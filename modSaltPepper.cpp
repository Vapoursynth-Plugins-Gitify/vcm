/*------------------------------------------------------------------------------------------

  saltPepper plugin for VapourSynth 

  Removes salt and or Pepper noise from image 

  Technique used in morphing literature is adapted to detect salt and pepper

  Any value outside plus minus tolerance limit is considered noise.

  Thread safe
  Last modified on 20 Aug 2020
  Copyright (C) <2006, 2007, 2009, 2020>  <V.C.Mohan>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    A copy of the GNU General Public License is at
    http://www.gnu.org/licenses/.
	For details of how to contact author see <http://www.avisynth.org/vcmohan/vcmohan.html> 
	

  Author V.C.Mohan

  Date 20 Aug 2020

  ----------------------------------------------------------------------------------------------
*/


	
//-----------------------------------------------------------------------------------	
	template <typename finc>
	void deSalt(finc * dp, int dpitch, const finc *sp, int spitch,
				int ht, int wd, finc tolr, bool average);
	template <typename finc>
	void dePepper(finc * dp, int dpitch, const finc *sp, int spitch,
				int ht, int wd, finc tolr, bool average);

	
/*************************************************
 * The following is the implementation 
 * of the defined functions.
 *************************************************/

//------------------------------------------------------------------------------------
template <typename finc>
void deSalt(finc * dp, int dpitch, const finc *sp, int spitch,
			 int ht, int wd, finc tolr, bool average)
{
	int offset[8];
//	finc tolr = tol;
			// calculate offsets once per frame to avoid multiplications later

	for ( int hh = -1, i = 0; hh <= 1; hh ++)
	{

		for(int w = -1; w <= 1; w ++)
		{

			if ( w == 0 && hh == 0)

				continue;

			else

				offset[ i ++] = hh * spitch + w ;
		}
	}


		// records new value only if it is a salt. Otherwise remains unchanged
	for ( int h = 1,hdpitch = dpitch, hspitch = spitch; h < ht - 1;
				h ++, hdpitch += dpitch, hspitch += spitch)
	{

		for (int w = 1; w < wd - 1; w ++)
		{

			finc val = sp[hspitch + w] - tolr;

			bool okay = true;

			for( int i = 0; i < 8; i ++)
			{
								
						// check if this point is a salt
				if( val <= sp[hspitch + w + offset[i]] )
				{
					
					okay = false;

					break;
				}
			}
			
			if (okay)
			{								// yes it is a salt. 
				
				if( average)
				{
						// replace by average
					float sum = 0;
		
					for( int i = 0; i < 8; i ++)
							
								// sum all points surrounding excluding itself	
							sum += sp[hspitch + w + offset[i]];
						
					
					dp[hdpitch + w] = sum / 8;
				}

				else
				{
						// we will replace by highest value of nearest neighbours
					finc  max = sp[hspitch + w + offset[0]];

					for( int i = 1; i < 8; i ++)
				
						if(max < sp[hspitch + w + offset[i]])
						
							max = sp[hspitch + w + offset[i]];							
				
								
						dp[hdpitch + w] = max;
				}
								
			}
						
		}
	}
	
}

//------------------------------------------------------------------------------------
template <typename finc>
void dePepper(finc * dp, int dpitch, const finc *sp, int spitch,
				int ht, int wd, finc tolr, bool average)
{

	int offset[8];
//	finc tolr = tol;
			// calculate offsets once per frame to avoid multiplications later

	for ( int hh = -1, i = 0; hh <= 1; hh ++)
	{

		for(int w = -1; w <= 1; w ++)
		{

			if ( w == 0 && hh == 0)

				continue;

			else

				offset[ i ++] = hh * spitch + w ;
		}
	}


		// records new value only if it is a pepper. Otherwise remains unchanged
	for ( int h = 1,hdpitch = dpitch, hspitch = spitch; h < ht - 1;
				h ++, hdpitch += dpitch, hspitch += spitch)
	{

		for (int w = 1; w < wd - 1; w ++)
		{

			finc val = sp[hspitch + w] + tolr;

			bool okay = true;

			for( int i = 0; i < 8; i ++)
			{								
						// check if this point is a pepper
				if( val >= sp[hspitch + w + offset[i]] )
				{					
					okay = false;

					break;
				}
			}			

			if(okay)
			{								// yes it is a pepper. so replace by average
				if( average)
				{
						// replace by average
					float sum = 0;
		
					for( int i = 0; i < 8; i ++)
							
								// sum all points surrounding excluding itself	
							sum += sp[hspitch + w + offset[i]];
						
					
					dp[hdpitch + w] = sum / 8;
				}

				else
				{
						// we will replace by lowest value of nearest neighbours
					finc  min = sp[hspitch + w + offset[0]];	// initialize to a value

					for( int i = 1; i < 8; i ++)
				
						if(min > sp[hspitch + w + offset[i]])
						
							min = sp[hspitch + w + offset[i]];							
				
								
						dp[hdpitch + w] = min;
				
				}
			}						
		}
	}	
}
//----------------------------------------------------------------------------------
typedef struct {
				VSNodeRef *node;
				const VSVideoInfo *vi;
				int planes[3];	//1 for de salt,2 for de pepper, 3 for both de salt & pepper
				int tol;			// tolerance				
				bool avg;			// whether average value or min / max value to be used for replacement
			
				}SaltPepperData;
//---------------------------------------------------------------------------------
// This function is called immediately after vsapi->createFilter(). This is the only place where the video
// properties may be set. In this case we simply use the same as the input clip. You may pass an array
// of VSVideoInfo if the filter has more than one output, like rgb+alpha as two separate clips.
static void VS_CC saltpepperInit
				(VSMap *in, VSMap *out, void **instanceData, 
				VSNode *node, VSCore *core, const VSAPI *vsapi) 
{
    SaltPepperData *d = (SaltPepperData *) * instanceData;
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
static const VSFrameRef *VS_CC saltpepperGetFrame
				(int n, int activationReason, void **instanceData, 
				void **frameData, VSFrameContext *frameCtx, 
				VSCore *core, const VSAPI *vsapi)
{
    SaltPepperData *d = (SaltPepperData *) * instanceData;

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
        VSFrameRef *dst = vsapi->copyFrame(src, core);

        // It's processing loop time!
        // Loop over all the planes
        
        for (int plane = 0; plane < fi->numPlanes; plane++)
		{
			

			const uint8_t *srcp = vsapi->getReadPtr(src, plane);
            int src_stride = vsapi->getStride(src, plane);
            uint8_t *dstp = vsapi->getWritePtr(dst, plane);
			int dst_stride = vsapi->getStride(dst, plane);
			int bht = vsapi->getFrameHeight(src, plane);            
            int bwd = vsapi->getFrameWidth(src, plane);
			int samplesize = fi->bytesPerSample;
			int spitch = src_stride / samplesize;
			int dpitch = dst_stride / samplesize;
			int nb = fi->bitsPerSample;
			

			if(  d->planes[plane] == 0)
				continue;

			if (fi->sampleType == stInteger)
			{
				if(fi->bitsPerSample == 8)
				{
					uint8_t min = 0;
					uint8_t max = (1<< nb) - 1;
					uint8_t tol = (max * d->tol) / 200;
					// 8 bit integer samples
			
					if( d->planes[plane] == 1 || d->planes[plane] == 3)
					{				

						deSalt(dstp, dpitch, srcp, spitch, bht, bwd,  tol, d->avg);

					}

					if( d->planes[plane] == 2 || d->planes[plane] == 3)
					{
				
						dePepper(dstp, dpitch, srcp, spitch, bht, bwd,  tol, d->avg);
					}

					
				}// 8 bit sample

				
				else // 16 bit integer samples			
				{
					uint16_t * dp = (uint16_t *) dstp;
					const uint16_t *sp = (uint16_t*)srcp;
					uint16_t max = ( 1 << nb) - 1;
					uint16_t tol = ((max * d->tol) / 200);

					if( d->planes[plane] == 1 || d->planes[plane] == 3)
					{				

						deSalt(dp, dpitch, sp, spitch, bht, bwd, tol, d->avg);

					}

					if( d->planes[plane] == 2 || d->planes[plane] == 3)
					{
				
						dePepper(dp, dpitch, sp, spitch, bht, bwd, tol, d->avg);
					}

					
				}// 16 bit sample
			}// integer samples
			else	// floating point data
			{
				float * dp = (float *) dstp;
				const float *sp = (float *)srcp;
				float min, max; 
				if(	fi ->colorFamily == cmRGB )
				{
					min = 0.0; 
					max = 1.0;
				}
				else if( plane == 0)
				{
					min = 0.0 ;
					max = 1.0;
				}
				else
				{
					min = -0.5;
					max =  0.5;
				}
				// now process
				float tol = (max * d->tol) / 200;

				if( d->planes[plane] == 1 || d->planes[plane] == 3)
				{				

					deSalt(dp, dpitch, sp, spitch, bht, bwd, tol, d->avg);

				}

				if( d->planes[plane] == 2 || d->planes[plane] == 3)
				{
				
					dePepper(dp, dpitch, sp, spitch, bht, bwd, tol, d->avg);
				}

					
			}//floating point data 
		} // plane
	

		 // Release the source frame
       vsapi->freeFrame(src);
        // A reference is consumed when it is returned so saving the dst ref somewhere
        // and reusing it is not allowed.
		
       return dst;
    }

    return 0;
}

//-------------------------------------------------------------------
// Free all allocated data on filter destruction
static void VS_CC saltpepperFree(void *instanceData, VSCore *core, const VSAPI *vsapi) 
{
   SaltPepperData *d = (SaltPepperData *)instanceData;
    vsapi->freeNode(d->node);	
    free(d);
}
//------------------------------------------------------------------------
// This function is responsible for validating arguments and creating a new filter
static void VS_CC saltpepperCreate(const VSMap *in, 
					VSMap *out, void *userData, VSCore *core,
					const VSAPI *vsapi) 
{
    SaltPepperData d;
    SaltPepperData *data;
    int temp;
    int err;

    // Get a clip reference from the input arguments. This must be freed later.
    d.node = vsapi->propGetNode(in, "clip", 0, 0);
    d.vi = vsapi->getVideoInfo(d.node);

    // Note that
    // vi->format can be 0 if the input clip can change format midstream.
	if (d.vi->format->colorFamily == cmCompat)
	{
		vsapi->setError(out, "saltPepper: Compat format not allowed.");
		vsapi->freeNode(d.node);
		return;
	}
	

    // If a property read fails for some reason (index out of bounds/wrong type)
    // then err will have flags set to indicate why and 0 will be returned. This
    // can be very useful to know when having optional arguments. Since we have
    // strict checking because of what we wrote in the argument string the only reason
    // this could fail is when the value wasn't set by the user.
    // And when it's not set we want it to default to enabled.
    for( int i = 0; i < 3; i ++)
	{
		d.planes[i] = vsapi->propGetInt(in, "planes", i, &err);

		if (err)
			d.planes[i] = 3;	// both salt and pepper

		// the only allowed 
		if (d.planes[i] < 0 || d.planes[i] > 3 ) 
		{
			vsapi->setError(out, "saltPepper: planes array  can have values of 0 to 3 only. 0 for no process, 1 for salt only, 2 for pepper only and 3 for both ");
			vsapi->freeNode(d.node);
			return;
		}
	}
	if(d.planes[0] == 0 && d.planes[1] == 0 && d.planes[2] == 0 || (d.vi->format->colorFamily == cmGray && d.planes[0] == 0) )
	{
			vsapi->setError(out, "saltPepper: all planes have values of 0. At least one plane must be non zero with a valid value.");
			vsapi->freeNode(d.node);
			return;
	}

	d.tol = vsapi->propGetInt(in, "tol", 0, &err);

    if (err)
        d.tol = 3;	// both salt and pepper

    // the only allowed 
    if (d.tol < 0 || d.tol > 5 ) 
	{
        vsapi->setError(out, "saltPepper: tol can have value of 0 to 5");
        vsapi->freeNode(d.node);
        return;
	}

	
	temp = !!vsapi->propGetInt(in, "avg", 0, &err);

    if (err)
        d.avg = true;

    // the only allowed values 
    if (temp < 0 || temp > 1 ) 
	{
        vsapi->setError(out, "saltPepper: avg can have a value of 0 or 1 only");
        vsapi->freeNode(d.node);
        return;
	}
	else if (temp == 0)
		d.avg = false;
	else
		d.avg = true;
	    // I usually keep the filter data struct on the stack and don't allocate it
    // until all the input validation is done.
	
    data =   (SaltPepperData * ) malloc(sizeof(d));
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
    vsapi->createFilter(in, out, "saltPepper", saltpepperInit, 
								saltpepperGetFrame, saltpepperFree, 
								fmParallel, 0, data, core);
    return;
}

/*
VS_EXTERNAL_API(void) VapourSynthPluginInit(VSConfigPlugin configFunc, VSRegisterFunction registerFunc, VSPlugin *plugin)
{
    configFunc("vc.mohan.dns", "SaltPepper", "VapourSynth SaltPepper", VAPOURSYNTH_API_VERSION, 1, plugin);
    registerFunc("saltPepper", "clip:clip;planes:int[]:opt;tol:int:opt;avg:int:opt", saltpepperCreate, 0, plugin);
}
*/
