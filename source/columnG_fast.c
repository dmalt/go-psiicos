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
 	plhs[0] = mxCreateDoubleMatrix((Nsen * Nsen + Nsen) / 2,1,mxREAL);
 	G_col = mxGetPr(plhs[0]);
 	columnG_fast(p, G_small, w, Nsen, Nsrc, G_col);
}

void columnG_fast(mwSize p, double * G_small, double w, mwSize Nsen, mwSize Nsrc, double * G_col)
{
	mwSize i = p % Nsrc;
	if(!i)
		i = Nsrc;
	i --;
	mwSize j = (p - i) / Nsrc;
	mwSize k, l, s = 0;
	for (k = 0; k < Nsen; ++k)
		for (l = k; l < Nsen; ++l)
		{
			G_col[s] = G_small[k + Nsen * i] * G_small[l + Nsen * j] * w;
			s++;
		}
}
