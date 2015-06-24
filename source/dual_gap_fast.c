#include "mex.h"
#include <cblas.h>
#include <stdio.h>

double max_array(double * a, mwSize num_elements);
void print_array(int a[], int num_elements);
double CalcDualGap(mwSize, mwSize, mwSize, mwSize, double, double *, double *, double *, double *);

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
 	double lambda 	= mxGetScalar(prhs[3]);
 	mwSize S 		= mxGetScalar(prhs[4]);
 	double * R 		= mxGetPr(prhs[5]);

 	plhs[0] = mxCreateDoubleScalar(mxREAL);
 	
 	double * eta = mxGetPr(plhs[0]);
 	*eta = CalcDualGap(Nsrc_sq, Ntime, Nsen_sq, S, lambda, M, G, X, R);
}

double CalcDualGap(mwSize Nsrc_sq, mwSize Ntime, mwSize Nsen_sq, mwSize S, double lambda,\
				 double * M, double * G, double * X, double * R)
{
	mwSize src, sen, t;
	double sigma = 0.;
	double * B = (double *)mxMalloc(sizeof(double) * S);
	double * R_ = (double *)mxMalloc(sizeof(double) * Nsen_sq * Ntime);
	double * temp = (double *)mxMalloc(Nsen_sq * Ntime * sizeof(double));
	double C = 0.;
	double * G_s = (double *)mxMalloc(sizeof(double) * Nsen_sq * 4); /* G_s corresponds to topography of two interacting sites on cortex, 
	i.e. two interacting pairs of dipoles*/
	double maxC_1 = 1.;
	double * X_s = (double *)mxMalloc(sizeof(double) * Ntime * 4);
	/*double * temp =  */
	for (src = 0; src < S; ++src)
	{
		G_s = &G[Nsen_sq * src * 4];
		X_s = &X[Ntime * src * 4];
		/* G(:,range)' * R: */
		cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 4, Ntime, Nsen_sq , 1, G_s, Nsen_sq, R, Nsen_sq, 0, temp, 4);
		B[src] = cblas_dnrm2(4 * Ntime, temp, 1);
		/*mexPrintf("%f\n", B[src]);*/
		sigma += cblas_dnrm2(4 * Ntime, X_s, 1); 
	}
	C = max_array(B, S) / lambda;
	/*mexPrintf("%f\n", C);
	mexPrintf("lambda = %f\n", lambda);
	mexPrintf("%f\n", max_array(B,S));
*/
	if(C > 1.)
		maxC_1 = C;
	cblas_dcopy(Nsen_sq * Ntime, R, 1, R_, 1);
	cblas_dscal(Nsen_sq * Ntime, 1. / maxC_1, R_, 1);

	double normR  = cblas_dnrm2(Nsen_sq * Ntime, R, 1);

	double normR_ = normR / maxC_1; 
	double R_dotM = cblas_ddot(Nsen_sq * Ntime, R_, 1, M, 1); /* = sum(sum(R_ .* M )) */
	/*mexPrintf("R_dotM = %f\n", R_dotM);
	mexPrintf("normR = %f\n", normR);
	mexPrintf("normR_ = %f\n", normR_);
	mexPrintf("sigma = %f\n", sigma);*/
	double eta = 0.5 * normR * normR + lambda * sigma + 0.5 * normR_ * normR_ - R_dotM;
	return eta;
	/*cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, Nsen_sq, Ntime, Nsrc_sq, 1, G, Nsen_sq, X, Nsrc_sq, 0, temp, Nsen_sq);*/
	/* cblas_dgemv(CblasColMajor, CblasNoTrans, Nsen_sq, Nsrc_sq, 1, G, Nsen_sq, X, 1, 0., temp, 1); */
}

