/******************************************************************************
FQSharp filter plugin for vapoursynth by V.C.Mohan
This filter operates in freq domain (2d) and blurs or unblurs
linear (motion) or circular (focus) styles

Author V.C.Mohan. 
June 2015

  Copyright (C) < 2014>  <V.C.Mohan>

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

********************************************************************************/


	void ApplyFilter(fftwf_complex * fout, float * Filter, int wd, int ht);

	template <typename finc>
	void getRealInput(float *data,const finc * fptr, int pitch,
						int wd, int ht, int wpad, int hpad, bool centered);

	
	template <typename finc>
	void getRealOutput(float *data, finc * fptr, int pitch,
						int wd, int ht, int wpad,  bool cent, finc min, finc max);

	


//	int getSign ( int h, int w);

	int getSign(int i);

	int DrawPSF(float *psf, bool linear, int xval, int yval, int bestx, int besty, float spike = 0.0);

	void DesignInverse(fftwf_complex * fout, float * Filter, float wn,
							int bestx, int besty,  float scale);

	
//-------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------
			//Draws the PSF linear or circular at center of the psf buffer 
		// draws psf in float data
int  DrawPSF(float *psf, bool linear, int xval, int yval, int bestx, int besty, float spike)
{

	int count;

				// zero psf area
	for(int h = 0; h < bestx * besty; h ++)
			
		psf[h] = 0.0;

	if( linear)
	{
		int length = abs( yval) > xval? abs(yval) : xval;
			
		count = 2 * length + 1;
		
					// draw the blur line at the center of frames best sizes
		
		if( abs(yval) > xval)
		{
			if(yval < 0)
			{
				yval = -yval;

				xval = - xval;
			}
			
			
			for(int h = - abs( yval); h <= abs( yval) ; h ++)
			{
			
				int w = (h * xval  ) / yval;		// get nearest integer

				int fraction = abs(h * xval  ) % abs(yval);		// get fractional part

					// distribute amp in ratio of fraction
			
				psf[ (besty/2 + h ) * bestx + ( bestx/2 + w)] 
					= (1.0f * (abs(yval) - fraction)) / (count * abs(yval));

				if( h * xval  > 0)

					psf[ (besty/2 + h ) * bestx + ( bestx/2 + w + 1)] 
						= (1.0f * fraction) / (count * abs(yval));

				else	// if( h  < 0)

					psf[ (besty/2 + h ) * bestx + ( bestx/2 + w - 1)] = (1.0f * fraction) / (count * abs(yval));				

			}
		}

			
		else	// xval is greater than yval

		{
				//xval is always a positive


			for(int w = - xval; w <= xval; w ++)
			{
			
				int h = (w * yval) / xval;
				
				int fraction = abs((w * yval  )) % xval;		// get fractional part

				psf[ (besty/2 + h ) * bestx + ( bestx/2 + w)] 
					= (1.0f * ( xval - fraction)) / (count * xval);

				if( w * yval > 0)

					psf[ (besty/2 + h + 1 ) * bestx + ( bestx/2 + w )] 
						= (1.0f * fraction) / (count * xval);

				else	// if( w * yval < 0)

					psf[ (besty/2 + h - 1 ) * bestx + ( bestx/2 + w )] 
						= (1.0f * fraction) / (count * xval);
							
			}

		}

	}

	else		// circular

	{
	
					// draw the blur circle at the center of frame

		count = 0;

		for(int h = - xval; h <= xval; h ++)

			for( int w = - xval; w <= xval; w ++)

				if( h * h + w * w <= xval * xval)

					count ++;
		
			// incase of deblur there is a spike added to center value
			// actually all values except center are reduced.

		for(int h = - xval; h <= xval; h ++)

			for( int w = - xval; w <= xval; w ++)

				if( h * h + w * w <= xval * xval)
									
					psf[ (besty/2 + h ) * bestx + ( bestx/2 + w)] = 1.0f / count;
				
					// center value will have spike to provide white noise and reduce instability
		if(spike > 0.001)

			psf[ (besty/2 ) * bestx + ( bestx/2 )] = (1.0f + spike)/ count;
	}

	return count;
}

//------------------------------------------------------------------------------------------------------------
void DesignInverse(fftwf_complex * fout, float * Filter, 
							  float wn,int bestx, int besty,
							   float scale)
{
	// the forward transform of PSF is in fout. Only real  positive values
				// get max value
	float mval = fout[0][0];
	

	for(int h = 0; h < bestx * besty; h ++)

		
				// get max value. Possibly the zeroth val is max but lets find
		if( (fout[ h ][0]) > mval)

				mval = fout[ h ][0];
		
	mval *= wn;		// this is the min value we will accept for inversion
	
					//  apply freq mask, and wn of inversion
	// PSF during forward transform is scaled up by sqrt(bestx * besty)
	// while inverting it it is 1/ this value. so no sscaling for transform is needed.
	float scaler = scale * (1.0 + wn)/ (bestx * besty ) /(1.0f - wn); // This trial/error derived approximation for scaler
									// as the higher wn is frequencies are inverted less
									// and we lose amplitude
					// we add white noise and  scaler . 
	for(int h = 0; h < bestx * besty; h ++)
	
		
		Filter[ h ] =  scaler / (fout[ h][0] + mval);

		

}

//------------------------------------------------------------------------------------------------------------
	
//---------------------------------------------------------------------------------------------
	// use this when data type is  complex 
int getSign(int i)
{
	return (i   & 1) == 0 ? 1 : -1;

}
/*// use this for real  data
int getSign(int h, int w)
{
	
	return ((h + (w >> 1)) & 1) == 0 ? 1 : -1;
}
*/

//---------------------------------------------------------------------------------------------
template <typename finc>
		// convert input unsigned char data to  float type
void getRealInput(float *data,const finc * fptr, 
						  int pitch, int wd, int ht,int bytes,
						  int wpad, int hpad, bool centered)
{
	// parameter bytes is for avisynth. not needed here
				// convert frame y values to float and keep in data  buffer
	float *data1 = data;
	const finc * fptr1 = fptr;
	
	if(centered)
	{
				// values multiplied by -1^(x+y) i.e sign to get the spectogram centered in frame
		for(int h = 0; h < ht; h ++)
		{

			for(int w = 0; w < wd; w ++)
			{
				
				data1[ w] = getSign( h, w) * fptr1[w];
			
			}
			data1 += wpad;
			fptr1 += pitch;
		}

	}

	else		// not centered. So as it is 
	{

		for(int h = 0; h < ht; h ++)
		{
			for(int w = 0;w < wd; w ++)
			{
			
				data1[ w] = fptr1[w];
			}
			data1 += wpad;
			fptr1 += pitch;
		}
			
	}

	data1 = data + ht * wpad;
			
		// fill with zeroes rest of data buffer
	for(int h = ht; h < hpad; h ++)
	{
		for(int w = 0; w < wpad;w ++)
		{

			data1[ w] = 0.0;
		}
		data1 += wpad;
	}

	data1 = data;

	for(int h = 0; h < hpad; h ++)
	{
		for(int w = wd;w < wpad; w ++)
		{

			data1[ w] = 0.0;
		}
		data1 += wpad;
	}

}


//---------------------------------------------------------------------------------------------
template <typename finc>
	// convert  float type output to unsigned char  data 
void getRealOutput(float *data, finc* fptr, int pitch,int wd, int ht,int wpad,int bytes, bool cent, finc min, finc max)
{
	// vapoursynth does not require bytes parameter
				// as values after fft may become -ve or go beyond 255 check and clip
	float *data1= data;
	finc * fptr1 = fptr;
	if(cent )
	{
		for(size_t h = 0; h < ht; h ++)
		{
			for(size_t w = 0; w < wd; w ++)
			{
				int val = data1[ w] * getSign(h, w);
				
				if(val < min)

					fptr1[w] = min;

				else if(val > max)

					fptr1[  w] = max;

				else

					fptr1[ w] = val;
			}
			data1 += wpad;
			fptr1 += pitch;
		}
	}

	else
	{
		for(size_t h = 0; h < ht; h ++)
		{
			for(size_t w = 0; w < wd; w ++)
			{
				
				if(data1[ w] < min)

					fptr1[ w] = min;

				else if(data1[ w] > max)

					fptr1[ w] = max;

				else

					fptr1[  w] = data1[ w];
			}
			data1 += wpad;
			fptr1 += pitch;
		}
	}
				
}
//---------------------------------------------------------------------------------------------


void hammingWindowing(float * cosBell, int pitch, int width, int height, int rfilt)
{
			
			// design a cosine bell hamming windowing function
	int rfiltsq = rfilt * rfilt;
//	float rfilt = sqrt(rfiltsq );	

	for(int h = 0; h < height/2; h ++)
	{
			for (int w = 0; w < width/2; w ++)
			{
				if( h * h + w * w <= rfiltsq)
				{
					float radial = sqrt((float) h * h + w * w );
					float bell =  0.46 + 0.54 * cos( (3.1416 * radial) / rfilt);

					for ( int hh = -1; hh <= 1; hh += 2)
					{
						for ( int ww = - 1; ww <= 1; ww += 2)
						{

							cosBell[ (height/2 + h * hh) * pitch + width/2 + w * ww ] *= bell;
						}
					}
				}
				
				else
				{
					for ( int hh = -1; hh <= 1; hh += 2)
					{
						for ( int ww = - 1; ww <= 1; ww += 2)
						{

							cosBell[ (height/2 + h * hh) * pitch + width/2 + w * ww ] *= 0.0;
						}
					}			
				}			
			}
	}
}

//-------------------------------------------------------------------------------------------------------


void ApplyFilter(fftwf_complex* fout, float * Filter, int wd, int ht )
{
			// applies the designed filter.in freq domain. Scalar multiply
			// of freq response of designed filter with freq transform of input

	//float scale = 1.0 / (hbest * wbest );
	int nval = ht * wd;
	for(int h = 0; h < nval; h++)
		
	{
		fout[h ][0] *= Filter[h ]; 
		fout[h ][1] *= Filter[h ]; 

	}
		
}
//-------------------------------------------------------------------------------------------------------------------------


