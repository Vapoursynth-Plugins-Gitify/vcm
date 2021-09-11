#pragma once
#ifndef STATS_AND_OFFSETS_LOOK_UP_TABLES_V_C_MOHAN
#define STATS_AND_OFFSETS_LOOK_UP_TABLES_V_C_MOHAN
//------------------------------------------------------------
//Look up tables for rectangular and circular grids
//Variance, average of grid
//-----------------------------------------------------------------
int makeLinearLUT(int* offsets, int pitch, int xcoord, int ycoord );
int makeRectGridLUT(int* offsets, int pitch, int xgrid, int ygrid = 0, int coffset = 0);
int makeCircularLUT(int* offsets, int pitch, int rad, int coffset = 0);
template <typename finc>
float getMeanValue(const finc* sp, int* offsets, int noff);
template <typename finc>
float getVariance(finc* sp, int* Offsets, int noff, float avg);


//---------------------------------------------------------------------

int makeLinearLUT(int* offsets, int pitch, int xcoord, int ycoord)
{
	int count = 0, npoints = 0;
	
	int absx = abs(xcoord), absy = abs(ycoord);

	if (absx < absy)
	{
		npoints = 2 * absy + 1;

		if (xcoord == 0)
		{
			// xcoord = 0. So avoid divide by zero
			for (int i = 0; i < npoints; i++)
			{
				offsets[i] = (i - npoints / 2) * pitch;
				count++;
			}
		}

		else
		{
			for (int i = 0; i < npoints; i++)
			{
				offsets[i] = (absy / ycoord) * (i - npoints / 2) * pitch
					+ absx / xcoord * (((i - npoints / 2) * absx
						+ (absy / 2)) / absy);
				count++;
			}
		}

	}
	else
	{
		npoints = 2 * absx + 1;
		if (ycoord == 0)	// avoid div by zero
		{

			for (int i = 0; i < npoints; i++)
			{
				offsets[i] = (i - npoints / 2);
				count++;
			}
		}

		else
		{
			for (int i = 0; i < npoints; i++)
			{
				offsets[i] = (absx / xcoord) * (i - npoints / 2)
					+ (absy / ycoord) * pitch * (((i - npoints / 2) * absy
						+ (absx / 2)) / absx);
				count++;
			}
		}
	}

	return count;
}
		//definitions
// makes offsets table from left top corner of grid and adds coffset to bring them to any point required
int makeRectGridLUT(int* offsets, int pitch, int xgrid, int ygrid, int coffset)
{
	if (ygrid == 0)
		ygrid = xgrid;

	int i = 0;

	for (int h = 0; h < ygrid; h++)
	{
		for (int w = 0; w < xgrid; w++)
		{
			offsets[i] = h * pitch + w  + coffset;
			i++;
		}
	}
	return i;
}
// makes offset table from center of circle and adds required offset to move ref to any point
int makeCircularLUT(int* cOff, int pitch, int rad, int coffset)
{
	int count = 0;
	int rsq = rad * rad;

	for (int h = -rad; h <= rad; h++)
	{
		int hsq = h * h;
		for (int w = -rad; w <= rad; w++)
		{
			if (hsq + (w * w) <= rsq)
			{
				cOff[count] = h * pitch + w  + coffset;
				count++;
			}

		}
	}
	return count;
}

// get mean value of a grid/ disc given offsets LUT pointer and number offsets
template <typename finc>
float getMeanValue(const finc* sp, int* offsets, int noff)
{
	float sum = 0;

	for (int i = 0; i < noff; i++)

		sum += (float)sp[offsets[i]];

	return (sum / noff);
}
// gives variance of given grid/ area with offsetsLUT, n offsets and average of grid / area
template <typename finc>
float getVariance(finc* sp, int* offsets, int noff, float avg)
{
	float sum = 0.0f;
	for (int i = 0; i < noff; i++)
	{
		sum += (avg - sp[offsets[i]]) * (avg - sp[offsets[i]]);
	}
	return sum / noff;
}
#endif
