#include "mex.h"

void columnG_fast(mwSize p, double * G_small, double w, mwSize Nsen, mwSize Nsrc, double * G_col);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if(nrhs != 3)
    	mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "Three inputs required.");
    mwSize p = mxGetScalar(prhs[0]);
 	double * G_small = mxGetPr(prhs[1]);
 	mwSize Nsen = mxGetM(prhs[1]);
 	mwSize Nsrc = mxGetN(prhs[1]);
 	double w = mxGetScalar(prhs[2]);
 	double * G_col;
 	plhs[0] = mxCreateDoubleMatrix(Nsen * Nsen,1,mxREAL);
 	G_col = mxGetPr(plhs[0]);
 	columnG_fast(p, G_small, w, Nsen, Nsrc, G_col);
}

void columnG_fast(mwSize p, double * G_small, double w, mwSize Nsen, mwSize Nsrc, double * G_col)
{
	mwSize r = (p - 1) % 4;
	mwSize s = (p - 1) / 4;
	mwSize Nsites = Nsrc / 2; 
	mwSize is = s % Nsites;
	
	mwSize js = (s - is) / Nsites;
	mwSize i, j;
	mwSize k, l;
	if(r == 0) 
	{
		i = 2 * is;
		j = 2 * js;
	}
	else if(r == 1)
	{
		i = 2 * is + 1;
		j = 2 * js;
	}
	else if(r == 2)
	{
		i = 2 * is;
		j = 2 * js + 1;
	}
	else if(r == 3)
	{
		i = 2 * is + 1;
		j = 2 * js + 1;
	}
	for (k = 0; k < Nsen; ++k)
		for (l = 0; l < Nsen; ++l)
			G_col[k + Nsen * l] = G_small[k + Nsen * i] * G_small[l + Nsen * j] * w;
}
