#include "mex.h"
#include <cblas.h>
#include <stdio.h>

double max_array(double * a, mwSize num_elements);
void print_array(int a[], int num_elements);
double CalcDualGap(mwSize, mwSize, mwSize, mwSize, double, double *, double *, double *, double *, double *);

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
    mwSize Nsrc_sq  = mxGetM(prhs[0]);
    mwSize S = Nsrc_sq / 4; /*Need to check that Nsr_sq is divisible by 4*/
	double * G 		= mxGetPr(prhs[0]);
 
 	mwSize Ntime  	= mxGetN(prhs[1]);
 	double * X 		= mxGetPr(prhs[1]);
 	double * M_  	= mxGetPr(prhs[2]);
 	double lambda 	= mxGetScalar(prhs[3]);
 	double epsilon 	= mxGetScalar(prhs[4]);
 	double mu 		= mxGetScalar(prhs[5]);

 	plhs[0] = mxCreateDoubleScalar(mxREAL);
 	plhs[1] = mxCreateDoubleMatrix(Nsrc_sq, Nsen_sq, mxREAL);
 	double * eta = mxGetPr(plhs[0]);
 	*eta = 1;
 	double * pr = mxGetPr(plhs[1]);
 	pr[0] = 2; pr[1] = 4;
}
