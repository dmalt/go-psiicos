#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if(nrhs != 5)
    	mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "Five inputs required.");
    mwSize Nsen_sq = mxGetM(prhs[0]);
 	mwSize Ntime = mxGetN(prhs[0]);
 	double * M = mxGetPr(prhs[0]);

    mwSize Nsen_sq = mxGetM(prhs[1]);
 	mwSize Ntime = mxGetN(prhs[1]);
 	double * G = mxGetPr(prhs[1]);

    mwSize Nsen_sq = mxGetM(prhs[2]);
 	mwSize Ntime = mxGetN(prhs[2]);
 	double * X = mxGetPr(prhs[2]);

 	mwSize Nsrc_sq = mxGetN(prhs[1]);
 	double * G_small = mxGetPr(prhs[1]);
 	double w = mxGetScalar(prhs[2]);
 	double * G_col;
 	plhs[0] = mxCreateDoubleMatrix(Nsen * Nsen,1,mxREAL);
 	G_col = mxGetPr(plhs[0]);
 	columnG_fast(p, G_small, w, Nsen, Nsrc, G_col);
}