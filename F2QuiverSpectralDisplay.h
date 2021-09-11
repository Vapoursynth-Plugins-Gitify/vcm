// F2Quiver Frequency spectrum a nd Filter Display functions
#ifndef F2QUIVER_FREQSPECTRAL_AND_FILTER_DISPLAY_VCMOHAN
#define F2QUIVER_FREQSPECTRAL_AND_FILTER_DISPLAY_VCMOHAN
//-------------------------------------------------------------------
// declarations
	template <typename finc>
	void drawHorizontalRuler2D(finc * dp , const int pitch, const int  wd ,int frqwidth, finc max, const int nyq);
	template <typename finc>
	void drawVerticalRuler2D(finc * dp , const int pitch, const int  ht , int hbest, finc max, const int nyq);
	template <typename finc>
	void displayPowerSpectrum2D(finc * dp, int pitch, float * powerSpectrum, int wd, int ht, int frqwidth,finc max, float gamma);
	template <typename finc>
	void displayFilter2D(finc * dp, int dpitch, float * filter,int fpitch, int wd , int ht, finc max, int nfft );
	void getPowerSpectrum2D( float * powerSpectrum, fftwf_complex * out, int nval);
//----------------------------------------------------------------------------
	// definitions

template <typename finc>
void displayFilter2D(finc * dp, int dpitch, float * filter,int fpitch, int wd , int ht, finc max, int nfft )
{
	for ( int h = 0; h < ht; h ++)
	{
		for (int w = 0; w < wd; w ++)
		{
			dp[w] = max * filter[w]; // *nfft;
		}
		filter += fpitch;
		dp += dpitch;
	}
}
		
	
//-------------------------------------------------------------------------------------------
template <typename finc>
void displayPowerSpectrum2D(finc * dp,int pitch, float * powerSpectrum, int wd, int ht, int frqwidth, finc max, float gamma)
{
	for ( int h = 0; h < ht; h ++)
	{
		for (int w = 0; w < wd ; w ++)
		{
			dp[w] = pow(powerSpectrum[w] , gamma) * max;
		}
		dp += pitch;
		powerSpectrum += frqwidth;
	}


}

//-------------------------------------------------------------------------------------------
void getPowerSpectrum2D( float * powerSpectrum, fftwf_complex * out, int nval)
{
	float pmax = 0.0f;

	for ( int i = 0; i < nval; i ++)
	{
		powerSpectrum[i] = out[i][0] * out[i][0] + out[i][1] * out[i][1];

		if ( pmax < powerSpectrum[i])

			pmax = powerSpectrum[i];
	}
	
	pmax = 1.0f / pmax;

	for ( int i = 0; i < nval; i ++)
	{
		powerSpectrum[i] *=  pmax;
	}

}
//---------------------------------------------------------------------------------------------------

template <typename finc>
void drawHorizontalRuler2D(finc * dstp , const int pitch, const int  wd , int frqwidth, finc max, const int nyq)
{
	if ( frqwidth/2 > wd )
		// can not display even half of ruler
		return;
	finc * dp = dstp + frqwidth / 2 ;	// ruler center point
	int limitr = wd - frqwidth / 2;
	int limitl = frqwidth / 2;

	for ( int fq = 0; fq < nyq/2 ; fq += 10)
	{
		int w = (fq * frqwidth) / (  nyq);	// ok. only ratio

	//	if( w >= wd / 2) continue;
		

		if ( (fq % 100) == 0)
		{
			// every 100 th freq
			if ( wd > 25)
			{
				for ( int h = 0; h < 10; h ++)
				{
					if ( w < limitr)
						dp[h * pitch + w] = max;
					if( w < limitl)
						dp[h * pitch - w] = max;
				}
			}
		}

		else if ( (fq % 50) == 0 && wd >= 50 &&  ( wd < 100 || wd >= 250) )
		{
			// every 50th freq
			for ( int h = 0; h < 8; h ++)
			{
				if ( w < limitr)
					dp[h * pitch + w] = (4 * max )/ 5;
				if ( w < limitl)
					dp[h * pitch - w] = ( 4 * max) / 5;
				
			}
		}

		else if ( (fq % 20) == 0 && wd >= 100 && wd <  250)
		{
			// every 20th freq if 50th not drawn
			for ( int h = 0; h < 8; h ++)
			{
				if ( w < limitr)
					dp[h * pitch + w] = (4 * max )/ 5;
				if ( w < limitl)
					dp[h * pitch - w] = ( 4 * max) / 5;
				
			}
		}

		else if( (fq % 10 ) == 0)
		{
			for ( int h = 0; h < 4; h ++)
			{
				if ( w < limitr)
					dp[h * pitch + w] = (3 * max )/ 4;
				if ( w < limitl)
					dp[h * pitch - w] = ( 3 * max) / 4;
				
			}
		}
	}
	// draw an arrow head to indicate center

		for ( int xy = 0; xy < 3; xy ++)
		{
	//		if ( xy < limitr)
				dp[ (10 - xy ) * pitch + xy] = max;
	//		if (  xy < limitl)
				dp[ (10 - xy ) * pitch - xy] = max;
		}
	
}

//-------------------------------------------------------------------------------
template <typename finc>
void drawVerticalRuler2D(finc * dsp , const int pitch, const int  ht, int hbest ,finc max, const int nyq)
{
	if ( hbest / 2 > ht)
		return;	// even center and half of ruler will not be visible
	finc * dp = dsp + hbest / 2 * pitch;	// center of ruler
	int limith = ht - hbest / 2;
	int limitl = hbest / 2;	// may not be reqd

	for ( int fq = 0; fq < nyq; fq += 10)
	{
		int h = (fq * hbest) / (  nyq);
	//	if ( h >= ht/2) continue;
		if ( (fq % 100) == 0)
		{
			// every 100 th freq
			if ( ht > 25)
			{
				for ( int w = 0; w < 10; w ++)
				{
					if ( h < limith)
						dp[h * pitch + w] = max;
					if( h < limitl)
						dp[ - h * pitch + w] = max;
				}
			}
		}

		else if ( (fq % 50) == 0 && ht >= 50 &&  ( ht < 100 || ht >= 250) )
		{
			// every 50th freq
			for ( int w = 0; w < 8; w ++)
			{
				if ( h < limith)
					dp[h * pitch + w] = (4 * max )/ 5;
				if ( h < limitl)
					dp[ - h * pitch + w] = ( 4 * max) / 5;
				
			}
		}

		else if ( (fq % 20) == 0 && ht >= 100 && ht <  250)
		{
			// every 20th freq if 50th not drawn
			for ( int w = 0; w < 8; w ++)
			{
				if ( h < limith)
					dp[h * pitch + w] = (4 * max )/ 5;
				if ( h < limitl)
					dp[ -h * pitch + w] = ( 4 * max) / 5;	
			}
		}

		else if( (fq % 10 ) == 0)
		{
			for ( int w = 0; w < 4; w ++)
			{
				if ( h < limith)
					dp[h * pitch + w] = (3 * max )/ 4;
				if ( h < limitl)
					dp[ -h * pitch + w] = ( 3 * max) / 4;
				
			}
		}
	}
	// draw an arrow head to indicate center

		for ( int xy = 0; xy < 3; xy ++)
		{
			if( xy < limith)
				dp[ (10 - xy ) + pitch * xy] = max;
			if ( xy < limitl)
			dp[ (10 - xy ) - pitch * xy] = max;
		}
	
}
//------------------------------------------------------------------------------------------------------

#endif