
/* --------------------------------------------------------------------------
Jitter function introduces a random or sinusoidal jitter to the image. 
The jitter  from frame to frame can be constant or vary randomly or linearly
Author V.C.Mohan
Created on 31 july 2020

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
//------------------------------------------------------------------------------ */


typedef struct {
	VSNodeRef *node;
	const VSVideoInfo *vi;
	int type;	//  jitter type 1.rand (random), 2.sine(sinusoid)
	int jmax;	// max amplitude of jitter to insert
	int wl;		// wwave length of sinusoid
	bool stat;	// jitter changes frame to frame
	int speed;	// speed of movement of sinusoid per frame  1.high, 2.med 3.low					
	int dense; // density of random shifts 1.high, 2.med 3.low


	int fmove;			// number of rows sinusoid moves per frame
	int interval;	// in this no shift occurs
	int * shift;		// buffer for shift values
	int modulo;	// mod of row to get shift val
}JitterData;
					
//--------------------------------------------------------------------------------------------
// functions for Jitter
	void makeRandomShiftLUT(int interval, int * lshift,int height, int jmax);
	void makeSineShiftLUT(int * sshift, int jmax, int wl);
	int	getShift( const int* lshift, const int row, int modulo);
	template <typename finc>
	void shiftRow( const finc * sp, finc * rp, int shiftval, int width);
	template <typename finc>
	void fillColor(finc * rp, int shiftval, finc color);

//--------------------------------------------------------------------------------------
	static void VS_CC jitterInit(VSMap *in, VSMap *out, void **instanceData, VSNode *node, VSCore *core, const VSAPI *vsapi)
	{
		JitterData *d = (JitterData *)* instanceData;
		vsapi->setVideoInfo(d->vi, 1, node);

		d->shift = NULL;

		if (d->type == 1)

		{
			d->shift = (int*)vs_aligned_malloc<int>(sizeof(int) * d->vi->height, 32);

			d->modulo = d->vi->height;

			// intervals with no jitter
			if (d->dense == 1)
			{

				d->interval = d->vi->height / 32;
			}

			else if (d->dense == 2)
			{
				d->interval = d->vi->height / 16;
			}
			else
			{
				d->interval = d->vi->height / 8;
			}

			makeRandomShiftLUT( d->interval, d->shift, d->vi->height, d->jmax);
		}
		else if (d->type == 2)
		{
			d->shift = (int *)vs_aligned_malloc<int>(sizeof(int) * d->wl, 32);

			d->modulo = d->wl;

			makeSineShiftLUT(d->shift, d->jmax, d->wl);

			if (d->speed == 1)

				d->fmove = d->wl / 16;

			else if (d->speed == 2)

				d->fmove = d->wl / 8;

			else
				d->fmove = d->wl / 4;
		}
		
	}

//----------------------------------------------------------------------
	// Free all allocated data on filter destruction
	static void VS_CC jitterFree(void *instanceData, VSCore *core, const VSAPI *vsapi)
	{
		JitterData *d = (JitterData *)instanceData;
		vsapi->freeNode(d->node);
		free(d);
		vs_aligned_free(d->shift);
	}
//------------------------------------------------------------------

int getShift(const int* lshift, const int row, int modulo)
{
	// calculate shift value for this row
	 return lshift[row % modulo];
	
}
//----------------------------------------------------------------------------
void makeSineShiftLUT(int * sshift, int jmax, int wl)
{
	for (int h = 0; h < wl; h++)
	{
		sshift[h] = int( jmax * ( ( 1.0f + sin(2 * M_PI * h / wl) )/ 2)  );
	}
}
//-----------------------------------------------------------------------------
void makeRandomShiftLUT( int interval, int * lshift, int height,  int jmax)
{
	// creates random shift at random row
	// initialize shift buffer to zero
	for (int k = 0; k < height; k++)
	{
		lshift[k] = 0;
	}
	// first rows with zero shift
	int initialZeroShift = rand() % interval;

	for (int k = initialZeroShift; k < height; )
	{

		lshift[k] = rand() % jmax;	// shift val for this row

		k += 1 + (rand() % interval);		// rows with zero shift
	}
}

//--------------------------------------------------------------------------
template <typename finc>
void shiftRow(const finc * sp, finc * rp,  int shiftval, int width)
{

	for(int i = shiftval; i < width; i++)
	{
		*(rp + i ) = *(sp + (i -shiftval ) );
	}
}
//-------------------------------------------------------------------------------------
template <typename finc>
void fillColor(finc * rp, int shiftval, finc color)
{
	for (int i = 0; i < shiftval; i++)
	{
		*(rp + i ) = color;
	}
}



/***************************************************************/
// This is the main function that gets called when a frame should be produced. It will, in most cases, get
// called several times to produce one frame. This state is being kept track of by the value of
// activationReason. The first call to produce a certain frame n is always arInitial. In this state
// you should request all the input frames you need. Always do it in ascending order to play nice with the
// upstream filters.
// Once all frames are ready, the filter will be called with arAllFramesReady. It is now time to
// do the actual processing.
static const VSFrameRef *VS_CC jitterGetFrame(int n, int activationReason, void **instanceData, void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi)
{
	JitterData *d = (JitterData *)* instanceData;

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
		const unsigned char* fp[] = { NULL, NULL, NULL, NULL };
		unsigned char * dp[] = { NULL, NULL, NULL, NULL };
		int pitch[4];
		
		int pht[4], pwd[4];
		for (int p = 0; p < nplanes; p++)
		{
			pitch[p] = vsapi->getStride(dst, p) / nbytes;
			dp[p] = vsapi->getWritePtr(dst, p);
			fp[p] = vsapi->getReadPtr(src, p);
			pht[p] = vsapi->getFrameHeight(src, p);
			pwd[p] = vsapi->getFrameWidth(src, p);
		}


		int * localShift = d->shift;	// will point to a fixed buf or a locally created buffer.
		int frameShift;

		if (d->type == 1)
		{
			// random shift
			frameShift = rand() % d->modulo;

		}
		else //sine
		{
			// either stationary or random buf
			frameShift = n * d->fmove;
		}

		int col[] = { 0, fi->colorFamily == cmYUV ? 127 : 0, fi->colorFamily == cmYUV ? 127 : 0 };

		
		for (int h = 0; h < height; h++)
		{
			int shiftValue =  getShift(localShift, h + frameShift, d->modulo);

			if (shiftValue != 0)
			{				

				for (int p = 0; p < nplanes; p++)
				{
					if (nbytes == 1)
					{

						shiftRow(fp[p], dp[p], shiftValue, width);
						fillColor(dp[p], shiftValue, (uint8_t)col[p]);
					}

					else if (nbytes == 2) // 16 bit samples
					{
						shiftRow((uint16_t*)fp[p], (uint16_t*)dp[p], shiftValue, width);

						fillColor((uint16_t*)dp[p], shiftValue, (uint16_t) (col[p] << (nbits - 8) ));
					}

					else
					{	// float
						shiftRow((float*)fp[p], (float*)dp[p], shiftValue, width);

						fillColor((float*)dp[p], shiftValue, 0.0f);
					}

				}		// for plane		

			}	// if shiftval != 0	
			for (int p = 0; p < nplanes; p ++)
			{
				fp[p] += pitch[p] * nbytes;
				dp[p] += pitch[p] * nbytes;
			}

		}	// for h =0;	

		return dst;
	}
	return 0;
}
	
// This function is responsible for validating arguments and creating a new filter
static void VS_CC jitterCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi)
{
	JitterData d;
	//   JitterData *data;
	int err;
	int temp;
	
	// Get a clip reference from the input arguments. This must be freed later.
	d.node = vsapi->propGetNode(in, "clip", 0, 0);
	d.vi = vsapi->getVideoInfo(d.node);

	// In this integer and float. Note that
	// vi->format can be 0 if the input clip can change format midstream.
	if (!isConstantFormat(d.vi)) // || d.vi->format->sampleType != stInteger || d.vi->format->bitsPerSample != 8)
	{
		vsapi->setError(out, "Jitter: format of clip must be constant");
		vsapi->freeNode(d.node);
		return;
	}

	if (d.vi->format->colorFamily != cmRGB && d.vi->format->colorFamily != cmYUV  && d.vi->format->colorFamily != cmGray)
	{
		vsapi->setError(out, "Jitter: only RGB or YUV or Gray color formats allowed");
		vsapi->freeNode(d.node);
		return;
	}

	if (d.vi->format->subSamplingH != 0 || d.vi->format->subSamplingW != 0)
	{
		vsapi->setError(out, "Jitter: color planes should have no subsampling");
		vsapi->freeNode(d.node);
		return;
	}
	

	temp = vsapi->propGetInt(in, "type", 0, &err);

	if (err)
		d.type = 1;
	else if ( temp < 1 || temp > 2)
	{
		vsapi->setError(out, "Jitter: type can have a value of 1 for random or 2 for sinuoidal jitter only");
		vsapi->freeNode(d.node);
		return;
	}
		d.type = temp;

	if (d.type == 2)
	{
		d.wl = vsapi->propGetInt(in, "wl", 0, &err);

		if (err)
			d.wl = d.vi->height / 8;

		else if (d.wl < 8 || d.wl > d.vi->height)
		{
			vsapi->setError(out, "Jitter: wavelength wl must be between 8 and frame height ");
			vsapi->freeNode(d.node);
			return;
		}
	}

	d.jmax = vsapi->propGetInt(in, "jmax", 0, &err);

	if (err)
		d.jmax = d.vi->width / 16;

	else if (d.jmax < 8 || d.jmax > d.vi->width / 4)
	{
		vsapi->setError(out, "Jitter: jmax the maximum amplitude of jitter should be between 8 and quarter of frame width ");
		vsapi->freeNode(d.node);
		return;
	}
	
	const char *den  = vsapi->propGetData(in, "dense", 0, &err);
	
	if (err)
		d.dense = 2;
	else if (strcmp(den,"high") == 0)		
		d.dense = 1;
	else if (strcmp(den,"med") == 0)		
		d.dense = 2;
	else if (strcmp(den,"low") == 0)		
		d.dense = 3;
	else 
	{
		vsapi->setError(out, "Jitter: dense can be high, med, or low/ only ");
		vsapi->freeNode(d.node);
		return;
	}
	
	temp = !!vsapi->propGetInt(in, "stat", 0, &err);
	if (err)
		d.stat = false;
	else if (temp == 0)
		d.stat = false;
	else
		d.stat = true;

	if (d.stat)
	{
		const char * speed = vsapi->propGetData(in, "speed", 0, &err);
		
		if (err)
			d.speed = 2;
		else if (strcmp(speed,"high") == 0)		
			d.speed = 1;
		else if (strcmp(speed,"med") == 0)		
			d.speed = 2;
		else if (strcmp(speed,"low") == 0)		
			d.speed = 3;
		else
		{
			vsapi->setError(out, "Jitter: speed can be high, med, or low only ");
			vsapi->freeNode(d.node);
			return;
		}
		
	}
	

	JitterData * data = (JitterData *)malloc(sizeof(d));

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
	vsapi->createFilter(in, out, "Jitter", jitterInit, jitterGetFrame, jitterFree, fmParallel, 0, data, core);
}

//////////////////////////////////////////
// Initial entry

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
configFunc("in.vcmohan.Jitter", "vcm", "VapourSynth Jitter plugin", VAPOURSYNTH_API_VERSION, 1, plugin);
registerFunc("Jitter", "clip:clip;type:int:opt;jmax:int:opt;dense:data:opt;stat:int:opt;speed:data:opt;", JitterCreate, 0, plugin);
}
*/
