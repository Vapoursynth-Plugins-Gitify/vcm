  
	// declarations
	template < typename finc>
	void drawVBoldLine( finc * dp, int length, int pitch, finc val);
	template < typename finc>
	void drawDashedLine( finc * dp, int length, int pitch, finc val, int ldash);
	template < typename finc>
	void drawDottedLine( finc * dp, int length, int pitch, finc val);
	template < typename finc>
	void DrawGridLines(finc * dp, int dpitch, int wd, int ht, int kb,
						int modbold, int moddash, int ginterval,
						finc bold, finc dash, finc dot,
						int ldash, int ldot);
	template < typename finc>
	void DrawCenteredGridLines(finc * dp, int dpitch, int wd, int ht, int kb,
		int modbold, int moddash, int ginterval,
		finc bold, finc dash, finc dot,
		int ldash, int ldot);

	template < typename finc>
	void DrawAxisScale( finc * dp, int dpitch, int wd, int ht, int kb,			  
							int modbold, int moddash, int ginterval,
							finc bold, finc dash, finc dot,
							int lbold, int ldash, int ldot);

	
/*************************************************
 * The following is the implementation 
 * of the defined functions.
 *************************************************/
//................................................
template < typename finc>
void drawVBoldLine( finc * dp, int length, int pitch, finc val)
{

	for ( int i = 0; i < length; i ++, dp += pitch)

		dp[0] = val;
}

//........................................................................................

template < typename finc>
void drawDashedLine( finc * dp, int length, int pitch, finc val, int ldash)
{

	for ( int l = 0, lpitch = 0; l < length; l += ldash + ldash)

		for( int m = 0, lpitch = l * pitch; m < ldash; m ++, lpitch += pitch)

			dp[lpitch] = val;
}

//.................................................................................................

template < typename finc>
void DrawGridLines(finc * dp, int dpitch, int wd, int ht, int kb,
						int modbold, int moddash, int ginterval,
						finc bold, finc dash, finc dot,
						int ldash, int ldot)
{
	int gridpitch = ginterval * dpitch;

	for(int h = 0, hdpitch = 0 ; h < ht; h += ginterval, hdpitch += gridpitch)
	{
	
		if((h % modbold) == 0)
		{
					// draw full lines 

			drawVBoldLine(dp + hdpitch, wd, kb, bold);	// luma
	
			continue;
		}

		else if((h % moddash) == 0)
		{
						// draw dashed line 4 pixels dash, 4 pixels gap

			drawDashedLine(dp + hdpitch, wd, kb, dash, ldash);
		
			continue;
		}
		
		else
					
					// dotted line
			drawDashedLine( dp + hdpitch, wd, kb, dot, ldot);
					
	}
}

//---------------------------------------------------------------------------------------------
template < typename finc>
void DrawCenteredGridLines(finc * dp, int dpitch, int wd, int ht, int kb,
	int modbold, int moddash, int ginterval,
	finc bold, finc dash, finc dot,
	int ldash, int ldot)
{
	// dp is placed on start of central grid line, ht is half height 
	int gridpitch = ginterval * dpitch;

	for (int h = 0, hdpitch = 0; h < ht; h += ginterval, hdpitch += gridpitch)
	{

		if ((h % modbold) == 0)
		{
			// draw full lines 

			drawVBoldLine(dp + hdpitch, wd, kb, bold);	// luma
			drawVBoldLine(dp - hdpitch, wd, kb, bold);	// luma

			continue;
		}

		else if ((h % moddash) == 0)
		{
			// draw dashed line 4 pixels dash, 4 pixels gap

			drawDashedLine(dp + hdpitch, wd, kb, dash, ldash);
			drawDashedLine(dp - hdpitch, wd, kb, dash, ldash);

			continue;
		}

		else

			// dotted line
			drawDashedLine(dp + hdpitch, wd, kb, dot, ldot);
			drawDashedLine(dp - hdpitch, wd, kb, dot, ldot);

	}
}

//.............................................................................................

template < typename finc>
void DrawAxisScale( finc * dp, int dpitch, int wd, int ht, int kb,			  
							int modbold, int moddash, int ginterval,
							finc bold, finc dash, finc dot,
							int lbold, int ldash, int ldot)
{
	int gridpitch = ginterval * dpitch;

	int center = ht / 2 * dpitch;

	for(int h = 0, hdpitch = 0; h < ht / 2; h += ginterval, hdpitch += gridpitch)
					
	if((h % modbold) == 0)
	{		
								// draw a horizontal line of length lbold
		drawVBoldLine(dp + center + hdpitch  - lbold /2  * kb, lbold, kb, bold);

		drawVBoldLine(dp + center - hdpitch  - lbold /2  * kb, lbold, kb, bold);

								
			continue;
	}
			
	else if((h % moddash)==0)
	{
		drawVBoldLine(dp + center + hdpitch  - ldash/2 * kb, ldash, kb, dash);

		drawVBoldLine(dp + center - hdpitch  - ldash/2 * kb, ldash, kb, dash);

		continue;

	}

	else
	{
					
		drawVBoldLine(dp + center + hdpitch  - ldot/2 * kb, ldot, kb, dot);

		drawVBoldLine(dp + center - hdpitch  - ldot/2 * kb, ldot, kb, dot);
	}
	
}

//...........................................................................................
void getComponentsRGB( uint8_t * cp, int RGB)
{
	cp[2] = RGB & 255;			// blue
	cp[1] = (RGB >> 8) & 255;	// green
	cp[0] = (RGB >> 16) & 255;	//red
}

void getYUVfromRGB( uint8_t * yp, int RGB)
{
	
	uint8_t blue = RGB & 255;	//blue
	uint8_t green = (RGB >> 8) & 255;	// green
	uint8_t red = (RGB >> 16) & 255;	// red

	yp[0] = 0.299*red+0.587*green+0.114*blue;	// y	
	yp[1] = -0.169*red-0.332*green+0.5*blue+128;	// u			
	yp[2] = 0.5*red-green*0.419-blue*0.0813+128;	// v
}
	
