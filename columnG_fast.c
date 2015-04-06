#include "mex.h"

void columnG_fast(mwSize p, double * G_small, double * w, mwSize Nsen, mwSize Nsrc, double * G_col);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if(nrhs != 3)
    	mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "Three inputs required.");
    mwSize p = mxGetScalar(prhs[0]);
 	double * G_small = mxGetPr(prhs[1]);
 	mwSize Nsen = mxGetM(prhs[1]);
 	mwSize Nsrc = mxGetN(prhs[1]);
 	double * w = mxGetPr(prhs[2]);
 	double * G_col;
 	plhs[0] = mxCreateDoubleMatrix(Nsen * Nsen * 2,1,mxREAL);
 	G_col = mxGetPr(plhs[0]);
 	columnG_fast(p, G_small, w, Nsen, Nsrc, G_col);
}

void columnG_fast(mwSize p, double * G_small, double * w, mwSize Nsen, mwSize Nsrc, double * G_col)
{
	mwSize q = p % (Nsrc * Nsrc);
	mwSize half = p / (Nsrc * Nsrc); /* We have columns of two types: those which end with zeros and those which start with zeros. half determines the type*/
	mwSize left = 0, right = 1;
	if(!q)
		q = Nsrc * Nsrc;
	mwSize i = q % Nsrc;
	if(!i)
		i = Nsrc;
	i --;
	mwSize j = (q - i) / Nsrc;
	double w_p = w[p-1];
	mwSize k, l;
	for (k = 0; k < Nsen; ++k)
		for (l = 0; l < Nsen; ++l)
		{
			if(half == left)
			{
				G_col[k + Nsen * l] = G_small[k + Nsen * i] * G_small[l + Nsen * j] * w_p;
				G_col[Nsen * Nsen + k + Nsen * l] = 0.;
			}
			else if(half == right)
			{
				G_col[k + Nsen * l] = 0.;
				G_col[Nsen * Nsen + k + Nsen * l] = G_small[k + Nsen * i] * G_small[l + Nsen * j] * w_p;
			}
		}
}
