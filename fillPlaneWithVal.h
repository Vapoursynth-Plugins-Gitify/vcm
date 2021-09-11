#ifndef FILLPLANEWITHVALUE_V_C_MOHAN
#define FILLPLANEWITHVALUE_V_C_MOHAN

template <typename finc>
void fillPlaneWithVal(finc * dp, int pitch, int wd, int ht, finc val);

template <typename finc>
void fillPlaneWithVal(finc * dp, int pitch, int wd, int ht, finc val)
{
	for (int h = 0; h < ht; h++)
	{
		for (int w = 0; w < wd; w++)
		{
			dp[w] = val;
		}

		dp += pitch;
	}
}
#endif