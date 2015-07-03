#include "mex.h"
#include <cblas.h>
#include <stdio.h>

float max_array(float* a, mwSize num_elements);
void print_array(int a[], int num_elements);
float CalcDualGap(mwSize, mwSize, mwSize, mwSize, float, float*, float*, float*, float*, float*);
mwSize BlockCoorDescent(mwSize, mwSize, mwSize, float* G, float* Y_p, float* M_, float lambda, float eps, float mu);

float max_array(float* a, mwSize num_elements)
{
   mwSize i;
   float max=-32000.;
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
    mwSize Nsrc_sq  = mxGetN(prhs[0]);
	float* G 		= mxGetPr(prhs[0]);
 
 	mwSize Ntime  	= mxGetN(prhs[1]);
 	float* X 		= mxGetPr(prhs[1]);
 	float* M_  	= mxGetPr(prhs[2]);
 	float lambda 	= mxGetScalar(prhs[3]);
 	float epsilon 	= mxGetScalar(prhs[4]);
 	float mu 		= mxGetScalar(prhs[5]);

 	plhs[0] = mxCreateDoubleScalar(mxREAL);
 	plhs[1] = mxCreateDoubleMatrix(Nsrc_sq, Nsen_sq, mxREAL);
 	float* eta = mxGetPr(plhs[0]);
 	*eta = 1;
 	float* pr = mxGetPr(plhs[1]);
 	pr[0] = 2; pr[1] = 4;
 	BlockCoorDescent(Nsrc_sq, Nsen_sq, Ntime, G, X, M_, lambda, epsilon, mu);
}



mwSize BlockCoorDescent(mwSize Nsrc_sq, mwSize Nsen_sq, mwSize Ntime,\
						 float* G, float* Y_p, float* M_, float lambda, float eps, float mu)
{
	mwSize maxIter = 1000000; /* One million */
	mwSize S = Nsrc_sq / 4; /*Need to check that Nsr_sq is divisible by 4*/
	float* Y_n = mxCalloc(Nsrc_sq * Ntime, sizeof(float));
	cblas_scopy(Nsrc_sq * Ntime, Y_p, 1, Y_n, 1);
	float* R = mxCalloc(Nsen_sq * Ntime, sizeof(float));
	cblas_scopy(Nsen_sq * Ntime, M_, 1, R, 1);
	cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans,\
				 Nsen_sq, Ntime, Nsrc_sq , -1., G, Nsen_sq, Y_p, Nsrc_sq, 1., R, Nsen_sq);
	mwSize i;
	for (i = 0; i < Nsen_sq * Ntime; ++i)
	{
		mexPrintf("%f\n", R[i]);
	}

	mexPrintf("BCD\n");
	mwSize iter, src;
	float* G_s = (float*)mxMalloc(sizeof(float) * Nsen_sq * 4);	
	float* Y_n_s = (float*)mxMalloc(sizeof(float) * Ntime * 4);
	for (iter = 0; iter < maxIter; ++iter)
	{
		for (src = 0; src < S; ++src)
		{
			cblas_scopy(4 * Ntime, Y_n, 1, Y_n_s, 1);
		}
	}
	return 0;
}
