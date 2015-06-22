#include "mex.h"
/*extern "C"
{
       
}*/
	#include <cblas.h>

void CalcDualGap(mwSize, mwSize, mwSize, double *, double *, double *, double *, double *);


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
 	plhs[0] = mxCreateDoubleMatrix(Nsen_sq,Ntime,mxREAL);
 	double * temp = malloc(sizeof(double));
 	temp = mxGetPr(plhs[0]);
 	CalcDualGap(Nsrc_sq, Ntime, Nsen_sq, M, G, X, R, temp);
}

void CalcDualGap(mwSize Nsrc_sq, mwSize Ntime, mwSize Nsen_sq, double * M, double * G, double * X, double * R, double * temp)
{
	/* cblas_dgemv(CblasColMajor, CblasNoTrans, Nsen_sq, Nsrc_sq, 1, G, Nsen_sq, X, 1, 0., temp, 1); */
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, Nsen_sq, Ntime, Nsrc_sq, 1, G, Nsen_sq, X, Nsrc_sq, 0, temp, Nsen_sq);
}

