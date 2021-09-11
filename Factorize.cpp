#ifndef FACTORIZE_V_C_MOHAN
#define FACTORIZE_V_C_MOHAN

void GetFactors(int n,int * facbuf);			
int getBestDim(int dimension, int *facbuf, int nfact=64);		
		

void GetFactors(int n,int * facbuf)
{
	// finds factors and fills facbuf with factor, remainder. maximum allowed or 32 pairs
	// factors are 4,2,3,5,....  remaining primes. 
	int p=4;
	double floor_sqrt;
	floor_sqrt = floor( sqrt((double)n) );

		 /*factor out powers of 4, powers of 2, then any remaining primes */
	 do 
	 {
		while (n % p)
		{
			switch (p) 
			{
				case 4: p = 2; break;
				case 2: p = 3; break;
				default: p += 2; break;
			}
			if (p > floor_sqrt)
				p = n;          /* no more factors, skip to end */
		}
		n /= p;
		*facbuf++ = p;
		*facbuf++ = n;
	} while (n > 1);
}

int getBestDim(int dimension, int *facbuf, int nfact)
{
	// returns nearest larger value having factors limited  2, 3, and 5
	int n = dimension;
	int largest = 7;
	size_t i;

	while (largest >5)
	{
		GetFactors(n,facbuf);
					// check value of largest factor
		for(i =0; i<nfact; i+=2)
			if(facbuf[i+1]==1)	// 
			{
				largest = facbuf[i];
				break;
			}

			if(largest >5)
				n += 4; 
	}
	return (n);
}

#endif 