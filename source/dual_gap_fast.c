#include "mex.h"
/*extern "C"
{
       
}*/
	#include <cblas.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if(nrhs != 6)
    	mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "Six inputs required.");

    mwSize Nsen_sq = mxGetM(prhs[0]);
 	mwSize Ntime = mxGetN(prhs[0]);
 	double * M = mxGetPr(prhs[0]);

 	mwSize Nsrc_sq = mxGetN(prhs[1]);
 	double * G = mxGetPr(prhs[1]);

 	double * X = mxGetPr(prhs[2]);

 	mwSize lambda = mxGetScalar(prhs[3]);

 	mwSize S = mxGetScalar(prhs[4]);

 	double * R = mxGetPr(prhs[5]);
}


CalcDualGap(mwSize Nsrc_sq, mwSize Ntime, mwSize Nsen_sq, double * M, double * G, double * X, double * R)
{
	double * temp = malloc(sizeof(double));
	cblas_dgemv(CblasRowMajor, CblasNoTrans, Ch, Ch, 1, colRmat[t], Ch, G[2 * i], 1, 0., temp, 1);
	return 0;
}
