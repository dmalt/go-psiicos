#include "mex.h"
#include <cblas.h>

void CalcDualGap(mwSize, mwSize, mwSize, mwSize, double *, double *, double *, double *, double *);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if(nrhs != 6)
    	mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "Six inputs required.");

    mwSize Nsen_sq  = mxGetM(prhs[0]);
 	mwSize Ntime  	= mxGetN(prhs[0]);
 	double * M  	= mxGetPr(prhs[0]);
 	mwSize Nsrc_sq 	= mxGetN(prhs[1]);
 	double * G 		= mxGetPr(prhs[1]);
 	double * X 		= mxGetPr(prhs[2]);
 	mwSize lambda 	= mxGetScalar(prhs[3]);
 	mwSize S 		= mxGetScalar(prhs[4]);
 	double * R 		= mxGetPr(prhs[5]);

 	plhs[0] = mxCreateDoubleMatrix(Nsen_sq, Ntime, mxREAL);
 	double * temp = malloc(Nsen_sq * Ntime * sizeof(double));
 	temp = mxGetPr(plhs[0]);
 	CalcDualGap(Nsrc_sq, Ntime, Nsen_sq, S, M, G, X, R, temp);
}

void CalcDualGap(mwSize Nsrc_sq, mwSize Ntime, mwSize Nsen_sq, mwSize S,\
				 double * M, double * G, double * X, double * R, double * temp)
{
	mwSize s;
	double sigma = 0.;
	double * B = malloc(sizeof(double) * S);
	double * G_s = malloc(sizeof(double) * Nsen_sq * 4); /* G_s corresponds to topography of two interacting sites on cortex, 
	i.e. two interacting pairs of dipoles*/
	double * X_s = malloc(sizeof(double) * Ntime * 4);
	/*double * temp =  */
	for (s = 0; s < S; ++s)
	{
		G_s = &G[Nsen_sq * s * 4];
		X_s = &X[Ntime * s * 4];
		cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, Nsen_sq, Ntime, 4, 1, G_s, Nsen_sq, R, Nsen_sq, 0, temp, Nsen_sq);
		B[s] = cblas_dnrm2(Nsen_sq * Ntime, temp, 1);
		sigma += cblas_dnrm2(4 * Ntime, X_s, 1 ); 
		mexPrintf("%f\n",B[s]);
	}
	mexPrintf("%f\n",sigma);

	/*cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, Nsen_sq, Ntime, Nsrc_sq, 1, G, Nsen_sq, X, Nsrc_sq, 0, temp, Nsen_sq);*/
	/* cblas_dgemv(CblasColMajor, CblasNoTrans, Nsen_sq, Nsrc_sq, 1, G, Nsen_sq, X, 1, 0., temp, 1); */
}

