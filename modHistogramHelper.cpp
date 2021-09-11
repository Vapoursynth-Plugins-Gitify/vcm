//---------------------------------------------------------------------
#ifndef HISTOGRAMADJUSTHELPER_CPP_V_C_MOHAN
#define HISTOGRAMADJUSTHELPER_CPP_V_C_MOHAN
// declarations
void getHistFromTable (const int* tbl,const int npoints, int nv, float * buf);
int  getHistCummTable(int* tbl, int npoints, int nv, float * buf);
void sigmaHist(float *buf, float  *buf2, int nval);
template <typename finc>
finc clamper(float val, finc min, finc max);


// as some values are bits per sample dependant could not use template
void getHistFromWindow8 ( const uint8_t * fp,int fpitch, 
				int nb, int wd, int ht, float * buf);

void getHistFromWindow16 ( const uint16_t * fp,int fpitch, 
				int nb, int wd, int ht, float * buf);

void getHistFromWindowf ( const float * fp,int fpitch, 
				int nb, int wd, int ht, float * buf);
//-----------------------------------------------------------------
void fillAdjustedValues8(const uint8_t * sp, uint8_t * dp, int pitch,
						int wd, int ht, int nb, float * buf, int limit );

void fillAdjustedValues16(const uint16_t * sp, uint16_t * dp, int pitch,
						int wd, int ht, int nb, float * buf, int limit );

void fillAdjustedValuesf(const float * sp, float * dp, int pitch,
						int wd, int ht, int nb, float * buf, int limit );
void fillMatchedValues8(const uint8_t * sp, uint8_t * dp, int pitch,
						int wd, int ht, int nb, float * hbuf,float * mbuf, int limit );
void fillMatchedValues16(const uint16_t * sp, uint16_t * dp, int pitch,
						int wd, int ht, int nb, float * hbuf,float * mbuf, int limit );
void fillMatchedValuesf(const float * sp, float * dp, int pitch,
						int wd, int ht, int nb, float * hbuf,float * mbuf, int limit );
int getMatchedValue( const float * mbuf, float val, int max);

template <typename finc>
finc clamper( float val, finc min, finc max);
	
/***************************************************************/
//------------------------------------------------------------------------------
int  getHistCummTable(int* tbl, int npoints, int nv, float * buf)
{
	// costructs histogram from table entries// pairs of %age level and  % agecummulative population 
	// level was restricted to 4096 maximum at input
	float cmult = 0.01f;	//value corresponding to a population of 1 
	float index = nv / 100.0f;

	int level = tbl[0] * index;
	int cpop = tbl[1] * cmult;

	int count = 0;

	if (level > 0)
	{
		int level = tbl[0] * index;
		int cpop = tbl[1] * cmult;
		// progressive sum of numbers is n * (n + 1) / 2
		float inc = 2.0 * cpop / (level + level * level);

		for (int j = 0; j < level; j++)
		{
			//
			buf[j] = j * inc;
			count++;
		}
	}

	else
	{
		buf[0] = level;

		count++;
	}

	for (int i = 2; i < npoints; i += 2)
	{
		// linearly interpolate %age population values for all inbetween luma levels = tbl[i - 2] * index;
		int level0 = tbl[i - 2] * index;
		int level1 = tbl[i] * index;
		float cpop0 = tbl[i - 1] * cmult;
		float cpop1 = tbl[1 + 1] * cmult;

		for (int j = level0; j < level1; j++)
		{
			// y = y1 + ((y2 - y1) * (x - x1)) / (x2 - x1)
			buf[j] = cpop0 + (j - level0) * (cpop1 - cpop0) / (level1 - level0);
			count++;

		}

	}

	if (tbl[npoints - 2] < 100 * index)
	{
		if (tbl[npoints - 1] < 100 * cmult)
		{
			int level0 = tbl[npoints - 2] * index;
			int level1 = nv;
			float cpop0 = tbl[npoints - 1] * cmult;
			float cpop1 = 100 * cmult;

			for (int j = level0; j < level1; j++)
			{
				// y = y1 + ((y2 - y1) * (x - x1)) / (x2 - x1)
				buf[j] = cpop0 + (j - level0) * (cpop1 - cpop0) / (level1 - level0);

				count++;

			}
		}

		else
		{
			// last point has full cummulative value;
			int level0 = tbl[npoints - 2] * index;
			int level1 = 100 * index;
			float cpop0 = tbl[npoints - 1] * cmult;


			for (int j = level0; j < nv; j++)
			{

				buf[j] = cpop0;
				count++;

			}
		}

	}

	return count;


}

//..............................................................................

void getHistFromWindow8 ( const uint8_t * fp,int fpitch, 
				int nb, int wd, int ht, float * buf)
{
		// populates buf with a normalized histogram. 
	
		// normalized increment for each count
	float norminc = 1.0f / ( wd * ht );
		// zero out buffer
	for(int h = 0; h < 1 << nb; h ++)
		
			buf[h] = 0;		
	
		// populate with count * norminc for each intensity value
	for(int h = 0; h < ht; h ++)
	{
		for(int w = 0; w < wd; w ++)
		{

			buf[ fp[ w] ] += norminc;
		}

		fp += fpitch;
	}
		

}
//....................................................................
void getHistFromWindow16 ( const uint16_t * fp, int fpitch, 
						int nb, int wd, int ht, float * buf)
{
		// populates buf with a normalized histogram.
	// normalized increment for each count
	float norminc = 1.0f / ( wd * ht );
	int shift =  nb - 12;
	if( nb <= 12)
	{		
		// zero out buffer
		for(int n = 0; n < 1 << nb; n ++)
		
			buf[n] = 0;	
	
		// populate with count * norminc for each intensity value
		for(int h = 0; h < ht; h ++)
		{
			for(int w = 0; w < wd; w ++)
			{

				buf[ fp[ w]  ] += norminc;
			}

			fp += fpitch;
		}
	}
	else
	{		
		// zero out buffer
		for(int n = 0; n < 4096; n ++)
		
			buf[n] = 0;	
	
		// populate with count * norminc for each intensity value
		for(int h = 0; h < ht; h ++)
		{
			for(int w = 0; w < wd; w ++)
			{

				buf[ fp[ w] >> shift] += norminc;
			}

			fp += fpitch;
		}
	}

}
//..............................................................................
void getHistFromWindowf ( const float * fp, int fpitch, 
				int nb, int wd, int ht, float * buf)
{
		// populates buf with a normalized histogram. 
	int max = (1 << nb) - 1;
		// normalized increment for each count
	float norminc = 1.0f / ( wd * ht );
		// zero out buffer
	for(int h = 0; h < (1 << nb); h ++)
		
			buf[h] = 0.0f;
		
//	int index;;
		// populate with count * norminc for each intensity value
	for(int h = 0; h < ht; h ++)
	{
		for(int w = 0; w < wd; w ++)
		{
		//	index = fp[ w] * max;	// as Y values are 0 to 1.0, we should get 0 to 4096
		
			buf[(int)(fp[ w] * max) ] += norminc; 
		}

		fp += fpitch;
	}		

}

//-------------------------------------------------------------------------------------
			// fill in front end
void getHistFromTable (const int* tbl,const int npoints,const int nv, float * buf)
{
		// costructs histogram from table entries// pairs of index and population values percentages
		// make index 4096
	float norminc = 0.1f / nv ;	//value corresponding to a population of 1 
	float index = nv/100.0f;
	
	for ( int i	= 0; i < npoints; i += 2)
	{
		if(i == 0)
		{
			// fill in the front part of buffer by extending first table value
			if( tbl[0] > 0)
			{

				for( int j = 0; j < (tbl[0] * index); j++)
				{

					buf[j]= tbl[1] * norminc;
				}
			}
		}
		else	// linearly interpolate population values for all inbetween luma values
		{
			for( int j = (tbl[i-2] * index); j < (tbl[i] * index); j ++)
			{
				// y = y1 + ((y2 - y1) * (x - x1)) / (x2 - x1)
				buf[j] = norminc * (tbl[i-1] + ((tbl[i+1] - tbl[i-1]) * (j - tbl[i-2]))/(tbl[i] - tbl[i-2]) ) ;
			}
		}
	}


	if(tbl[npoints-2] < 100)
	{
			// fill last part of buffer by extending last value
		for(int j = (tbl[npoints - 2]* index); j < nv; j++)
		{

			buf[j] = tbl[npoints - 1] * norminc;

			
		}
	}
	//normalize
	float sum = 0.0;

	for ( int i = 0; i< nv; i ++)

		sum += buf[i];

	for ( int i = 0; i < nv; i ++)

		buf[i] /= sum;


	
}
					
//..............................................................................

		
void sigmaHist(float *buf, float  *buf2, int index)
{
	// makes a progressive sum 	and normalizes
//	float sum = buf[0];
//	int index = nb == 0? 4096 : 1 << nb;
	buf2[0] = buf[0];

	for( int i = 1; i < index; i++)
	{
		
	
		buf2[i] = buf2[i -1] + buf[i]; 


	}

	// normalize buffer

	for(int i = 0; i < index; i ++)
	{
		buf2[i] /= buf2[index - 1];
	}

}
//-----------------------------------------------------------------------------------
template <typename finc>
finc clamper( float val, finc min, finc max)
{
	return val < min? 0 : val > max? max : val;
}
//-----------------------------------------------------------------------------------
void fillAdjustedValues8(const uint8_t * sp, uint8_t * dp, int pitch,
						int wd, int ht, int nb, float * buf, int limit )
{
	for( int h = 0; h <  ht; h ++)
	{
		for ( int w =  0; w <  wd; w ++)
		{
			uint8_t temp = buf[sp[w]] * 255;
			dp[w] = sp[w] + ((temp - sp[w] ) * (100 - limit)) / 100; 
			
		}

		dp += pitch;
		sp += pitch;
	}

}


//------------------------------------------------------------------------------
void fillAdjustedValues16(const uint16_t * sp, uint16_t * dp, int pitch,
						int wd, int ht, int nb, float * buf, int limit )
{
	int max = nb <= 12 ? (1 << nb) - 1 : 4095;
	if (nb <= 12)
	{
		for( int h = 0; h <  ht; h ++)
		{
			for ( int w =  0; w <  wd; w ++)
			{
				
				uint16_t temp = buf[sp[w]] * max;
				dp[w] = sp[w] + ((temp - sp[w] ) * (100 - limit)) / 100; 
				
			}

			dp += pitch;
			sp += pitch;
		}
	}

	else
	{
		for( int h = 0; h <  ht; h ++)
		{
			for ( int w =  0; w <  wd; w ++)
			{
				uint16_t temp = buf[sp[w] >> (nb - 12)] * (1  << (nb ));
				
				dp[w] = sp[w] + ((temp - sp[w] ) * (100 - limit)) / 100; 
				
			}

			dp += pitch;
			sp += pitch;
		}
	}

}

//------------------------------------------------------------------------------
void fillAdjustedValuesf(const float * sp, float * dp, int pitch,
						int wd, int ht, int nb, float * buf, int limit )
{
	int max = (1 << nb) - 1;
//	buf += 2048;
	float mult = (100.0f - limit) / 100.0f;
	for( int h = 0; h <  ht; h ++)
	{
		for ( int w =  0; w <  wd; w ++)
		{
			
			float temp = buf[int (sp[w] * max)] ;
			dp[w] = sp[w] + (temp - sp[w] ) * mult;
			
		}

		dp += pitch;
		sp += pitch;
	}

}

//-------------------------------------------------------------------------------------------------------------------------

int getMatchedValue( const float * mbuf, float val, int max)
{
	int start = 1.0 / val;

	if( mbuf[start] == val)
		return start;

	for ( int i = start; i < max; i ++)
	{
		if ( i > start && mbuf[i] >= val )

			return i;

		else if( i == start && mbuf[i] > val)
			// we are searching in wrong direction
			break;
	}

	for ( int i = start; i >= 0; i --)
	{
		if ( mbuf[i] <= val)

			return i;
	}

	return 0;
}

//--------------------------------------------------------------------------------------
void fillMatchedValues8(const uint8_t * sp, uint8_t * dp, int pitch,
						int wd, int ht, int nb, float * hbuf,float * mbuf, int limit )
{
	int max = 256;
	
	float temp;
	if ( nb <= 12)
	{
		for( int h = 0; h <  ht; h ++)
		{
			for ( int w =  0; w <  wd; w ++)
			{
				temp = hbuf[ sp[w] ];
				uint8_t val = getMatchedValue( mbuf, temp, max) ;
		
				dp[w] = sp[w] + ((val -sp[w]  ) * (100 - limit)) / 100;
			}

			dp += pitch;
			sp += pitch;
		}
	}
}
//--------------------------------------------------------------------------------------
void fillMatchedValues16(const uint16_t * sp, uint16_t * dp, int pitch,
						int wd, int ht, int nb, float * hbuf,float * mbuf, int limit )
{
	int max = nb <= 12 ? 1 << nb : 4096;
	float mult = 1.0f / max;
	float temp;
	if ( nb <= 12)
	{
		for( int h = 0; h <  ht; h ++)
		{
			for ( int w =  0; w <  wd; w ++)
			{
				temp = hbuf[ sp[w] ];
				uint16_t val = getMatchedValue( mbuf, temp, max) ;
		
				dp[w] = sp[w] + (( val -sp[w]  ) * (100 - limit)) / 100;
			}

			dp += pitch;
			sp += pitch;
		}
	}

	else
	{
		for( int h = 0; h <  ht; h ++)
		{
			for ( int w =  0; w <  wd; w ++)
			{
				temp = hbuf[sp[w] >> (nb - 12) ];
				uint16_t val = getMatchedValue( mbuf, temp, max) << (nb - 12) ;
		
				dp[w] =  sp[w] + (( val - sp[w] ) * (100 - limit)) / 100;
			}

			dp += pitch;
			sp += pitch;
		}
	}

}
//--------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
void fillMatchedValuesf(const float * sp, float * dp, int pitch,
						int wd, int ht, int nb, float * hbuf,float * mbuf, int limit )
{
	int max = 4096;
	float mult = 1.0f / max;
	float temp;

	for( int h = 0; h <  ht; h ++)
	{
		for ( int w =  0; w <  wd; w ++)
		{
			temp = hbuf[int ( sp[w] * max) ];
			float val = (getMatchedValue( mbuf,temp,max) * mult);
		
			dp[w] = sp[w] + ((val - sp[w]  ) * (100 - limit)) / 100;
		}

		dp += pitch;
		sp += pitch;
	}

}
//--------------------------------------------------------------------------------------------------------
#endif