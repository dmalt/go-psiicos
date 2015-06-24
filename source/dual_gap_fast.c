#include "mex.h"
#include <cblas.h>
#include <stdio.h>

double max_array(double * a, mwSize num_elements);
void print_array(int a[], int num_elements);
void CalcDualGap(mwSize, mwSize, mwSize, mwSize, double, double *, double *, double *, double *, double *);

double max_array(double * a, mwSize num_elements)
{
   mwSize i;
   double max=-32000.;
   for (i=0; i<num_elements; i++)
	 if (a[i]>max)
	    max=a[i];
   return(max);
}

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
 	CalcDualGap(Nsrc_sq, Ntime, Nsen_sq, S, lambda, M, G, X, R, temp);
}

void CalcDualGap(mwSize Nsrc_sq, mwSize Ntime, mwSize Nsen_sq, mwSize S, double lambda,\
				 double * M, double * G, double * X, double * R, double * temp)
{
	mwSize src, sen, t;
	double sigma = 0.;
	double * B = malloc(sizeof(double) * S);
	double * R_ = malloc(sizeof(double) * Nsen_sq * Ntime);
	double C = 0;
	double * G_s = malloc(sizeof(double) * Nsen_sq * 4); /* G_s corresponds to topography of two interacting sites on cortex, 
	i.e. two interacting pairs of dipoles*/
	double maxC_1 = 1;
	double * X_s = malloc(sizeof(double) * Ntime * 4);
	/*double * temp =  */
	for (src = 0; src < S; ++src)
	{
		G_s = &G[Nsen_sq * src * 4];
		X_s = &X[Ntime * src * 4];
		/* G(:,range)' * R: */
		cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, Nsen_sq, Ntime, 4, 1, G_s, Nsen_sq, R, Nsen_sq, 0, temp, Nsen_sq);
		B[src] = cblas_dnrm2(Nsen_sq * Ntime, temp, 1);
		sigma += cblas_dnrm2(4 * Ntime, X_s, 1); 
		mexPrintf("%f\n",B[src]);
	}
	C = max_array(B, S) / lambda;
	if(C > 1)
		maxC_1 = C;
	cblas_dcopy(Nsen_sq * Ntime, R, 1, R_, 1);
	cblas_dscal(Nsen_sq * Ntime, 1. / maxC_1, R_, 1);
	mexPrintf("%f\n",C);

	/*cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, Nsen_sq, Ntime, Nsrc_sq, 1, G, Nsen_sq, X, Nsrc_sq, 0, temp, Nsen_sq);*/
	/* cblas_dgemv(CblasColMajor, CblasNoTrans, Nsen_sq, Nsrc_sq, 1, G, Nsen_sq, X, 1, 0., temp, 1); */
}

