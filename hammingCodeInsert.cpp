// code for inserting to apply hamming filter
// 	// go back to spactial domain after centering

fftwf_complex* outp = d->outBuf;
float* fqFilter = d->FreqFilter;
int bwdR = d->frqwidth;
// transfer to outBuf and apply -1^n to center output after fft
int start = 1;

for (int h = 0; h < d->hbest; h++)
{
	int wstart = start;

	for (int w = 0; w < bwdR; w++)
	{
		outp[w][0] = wstart * fqFilter[w];

		outp[w][1] = 0.0;

		wstart = -wstart;
	}
	outp += bwdR;
	fqFilter += bwdR;
	start = -start;

}

d->fftwf_execute(d->pinv);

// hamming window as per filter radius in spatial domain
F2QhammingWindowing(d->inBuf, d->wbest, d->wbest, d->hbest, d->frad);

// forward transform into freq domain

d->fftwf_execute(d->pf);

// zero phase and also scale as it has gone twice transform

float scaler = 1.0f / (d->hbest * d->wbest);

for (int i = 0; i < bwdR * d->hbest; i++)

	d->FreqFilter[i] = scaler * sqrt(getAmpSquareOfComplex( d->outBuf+i));


// now atlast we are ready to process. 
