#include "mex.h"
#include <cblas.h>
#include <stdio.h>
#include <stdlib.h>
#include "math.h"

double max_array(double* a, mwSize num_elements);
void print_array(int a[], int num_elements);
double CalcDualGap(mwSize, mwSize, mwSize, double, double *, double *, double *, double *);
double BlockCoorDescent(mwSize, mwSize, mwSize,\
 double* G, double* rY_p, double* iY_p, double* rM_, double* iM_, double lambda, double eps, double mu,\
  double* rYout, double* iYout);

double max_array(double* a, mwSize num_elements)
{
   mwSize i;
   double max=-32000.;
   for (i=0; i<num_elements; i++)
	 if (a[i]>max)
	    max=a[i];
   return(max);
}

double fmax(double a, double b)
{
	if(a>=b)
		return a;
	else
		return b;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if(nrhs != 6)
    	mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "Six inputs required.");
    mwSize Nsen_sq  = mxGetM(prhs[0]);
    mwSize Nsrc_sq  = mxGetN(prhs[0]);
	double* G 		= mxGetPr(prhs[0]);
 
 	mwSize Ntime  	= mxGetM(prhs[1]);
 	double* reY 	= mxGetPr(prhs[1]);
 	double* imY 	= mxGetPi(prhs[1]);
 	/*mwSize i;
 	for (i = 0; i < 4 * Ntime; ++i)
	{
		mexPrintf("cblas_dcopy: imY = %f\n", imY[i]);
	}*/
 	double* reM_  	= mxGetPr(prhs[2]);
 	double* imM_  	= mxGetPi(prhs[2]);
 	double lambda 	= mxGetScalar(prhs[3]);
 	double epsilon 	= mxGetScalar(prhs[4]);
 	double mu 		= mxGetScalar(prhs[5]);

 	double* rYout;
 	double* iYout;
 	plhs[0] = mxCreateDoubleMatrix(Ntime, Nsrc_sq, mxCOMPLEX);
 	rYout = mxGetPr(plhs[0]);
 	iYout = mxGetPi(plhs[0]);

 	plhs[1] = mxCreateDoubleScalar(mxREAL);
 	double* iter = mxGetPr(plhs[1]);
 	/*double* pr = mxGetPr(plhs[1]);
 	pr[0] = 2; pr[1] = 4;*/
 	*iter = BlockCoorDescent(Nsrc_sq, Nsen_sq, Ntime, G, reY, imY, reM_, imM_, lambda, epsilon, mu, rYout, iYout);
}



double BlockCoorDescent(mwSize Nsrc_sq, mwSize Nsen_sq, mwSize Ntime,\
						 double* G, double* rY_p, double* iY_p, double* rM_, double* iM_, double lambda, double eps, double mu,\
						 double* rYout, double* iYout)
{
	mwSize i;
	mwSize maxIter = 1e6;  /*One million */
	mwSize S = Nsrc_sq / 4; /*Need to check that Nsr_sq is divisible by 4*/

	double* Y_n = mxCalloc(Nsrc_sq * Ntime * 2, sizeof(double));
	double* Y_p = mxCalloc(Nsrc_sq * Ntime * 2, sizeof(double));
	
	if(iY_p == NULL)
	{
		iY_p = mxCalloc(Nsrc_sq * Ntime, sizeof(double));
		for ( i = 0; i < Ntime *  Nsrc_sq; ++i)
			iY_p[i] = 0.;
	}
/* Concatenate real and imag parts of Y; need a loop because we input Y transposed and it gets tricky to cat real and imag parts*/
	for (i = 0; i < Nsrc_sq; ++i)
	{		
		cblas_dcopy(Ntime, &rY_p[i * Ntime], 1, &Y_n[2 * i * Ntime], 1);
		cblas_dcopy(Ntime, &iY_p[i * Ntime], 1, &Y_n[(2 * i + 1) * Ntime], 1);
		cblas_dcopy(Ntime, &rY_p[i * Ntime], 1, &Y_p[2 * i * Ntime], 1);
		cblas_dcopy(Ntime, &iY_p[i * Ntime], 1, &Y_p[(2 * i + 1) * Ntime], 1);
	}

	double* M_ = mxCalloc(Nsen_sq * Ntime * 2, sizeof(double));
	cblas_dcopy(Nsen_sq * Ntime, rM_, 1, &M_[0], 1);
	cblas_dcopy(Nsen_sq * Ntime, iM_, 1, &M_[Nsen_sq * Ntime], 1);

	double* R = mxCalloc(Nsen_sq * Ntime * 2, sizeof(double));
	cblas_dcopy(Nsen_sq * Ntime, rM_, 1, &R[0], 1);
	cblas_dcopy(Nsen_sq * Ntime, iM_, 1, &R[Nsen_sq * Ntime], 1);

	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, Nsen_sq, Ntime * 2, Nsrc_sq , -1.,  G, Nsen_sq, Y_p, Ntime * 2,  1., R, Nsen_sq);
	
	mexPrintf("BCD\n");
	mwSize iter, src;
	double* G_s = (double*)mxMalloc(sizeof(double) * Nsen_sq * 4);	
	double* Y_n_s = (double*)mxMalloc(sizeof(double) * 2 * Ntime * 4);
	double* Y_p_s = (double*)mxMalloc(sizeof(double) * 2 * Ntime * 4);
	double* delta_Y = (double*)mxMalloc(sizeof(double) * 2 * Ntime * 4);
	double eta = 2.;
	double Y_n_s_norm = 0.;
	double scale = 0.;
	for (iter = 0; iter < maxIter; ++iter)
	{
		src = iter % S;	
		cblas_dcopy(4 * Ntime * 2, &Y_n[4 * Ntime * 2 * src], 1, &Y_n_s[0], 1);
		cblas_dcopy(4 * Ntime * 2, &Y_n[4 * Ntime * 2 * src], 1, &Y_p_s[0], 1);
		cblas_dcopy(4 * Nsen_sq, &G[4 * Nsen_sq * src], 1, G_s, 1);

		cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, Ntime * 2, 4, Nsen_sq, mu,  R, Nsen_sq, G_s, Nsen_sq,  1., Y_n_s, Ntime * 2);

		Y_n_s_norm = cblas_dnrm2(2 * Ntime * 4, Y_n_s, 1);
		scale = fmax(1 - mu * lambda / Y_n_s_norm, 0);

		cblas_dscal(2 * Ntime * 4, scale, Y_n_s, 1);

		cblas_dcopy(4 * Ntime * 2, Y_n_s, 1, delta_Y, 1);

		cblas_daxpy(2 * Ntime * 4, -1., Y_p_s, 1, delta_Y, 1);

		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, Nsen_sq, Ntime * 2, 4 , -1.,  G_s, Nsen_sq, delta_Y, Ntime * 2,  1., R, Nsen_sq);

		cblas_dcopy(4 * Ntime * 2, Y_n_s, 1, &Y_n[4 * Ntime * 2 * src], 1);

		if(!(iter % 20 * S))
			eta  = CalcDualGap( Ntime * 2, Nsen_sq, S, lambda, M_, G, Y_n, R);	

		if(eta < eps)
			break;
	}
	for (i = 0; i < Nsrc_sq; ++i)
	{		
		cblas_dcopy(Ntime, &Y_n[2 * i * Ntime], 1, &rYout[i * Ntime], 1);
		cblas_dcopy(Ntime, &Y_n[(2 * i + 1) * Ntime], 1, &iYout[i * Ntime], 1);
	}
	return iter;
} /* --- BlockCoorDescent --- */


double CalcDualGap(mwSize Ntime, mwSize Nsen_sq, mwSize S, double lambda,\
				 double * M, double * G, double * X, double * R)
{
	/*mexPrintf("Ntime = %d, Nsen_sq = %d, S = %d, lambda = %f ", Ntime, Nsen_sq, S, lambda );*/
	mwSize i;
	mwSize src, sen, t;
	double sigma = 0.;
	double * B = (double *)mxMalloc(sizeof(double) * S);
	double * R_ = (double *)mxMalloc(sizeof(double) * Nsen_sq * Ntime);
	double * temp = (double *)mxMalloc(Nsen_sq * Ntime * sizeof(double));
	double C = 0.;
	double * dG_s; /*= (double *)mxMalloc(sizeof(double) * Nsen_sq * 4);*/ /* dG_s corresponds to topography of two interacting sites on cortex, 
	i.e. two interacting pairs of dipoles*/
	double maxC_1 = 1.;
	double * X_s;/*= (double *)mxMalloc(sizeof(double) * Ntime * 4);*/
	/*double * temp =  */
	for (src = 0; src < S; ++src)
	{
		dG_s = &G[Nsen_sq * src * 4];
		X_s = &X[Ntime * src * 4];
		/* G(:,range)' * R: */
		cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 4, Ntime, Nsen_sq , 1., dG_s, Nsen_sq, R, Nsen_sq, 0., temp, 4);
		B[src] = cblas_dnrm2(4 * Ntime, temp, 1);
		sigma += cblas_dnrm2(4 * Ntime, X_s, 1); 
	}
	C = max_array(B, S) / lambda;
	if(C > 1.)
		maxC_1 = C;
	cblas_dcopy(Nsen_sq * Ntime, R, 1, R_, 1);
	cblas_dscal(Nsen_sq * Ntime, 1. / maxC_1, R_, 1);

	double normR  = cblas_dnrm2(Nsen_sq * Ntime, R, 1);

	double normR_ = normR / maxC_1; 
	double R_dotM = cblas_ddot(Nsen_sq * Ntime, R_, 1, M, 1); /* = sum(sum(R_ .* M )) */
	double eta = 0.5 * normR * normR + lambda * sigma + 0.5 * normR_ * normR_ - R_dotM;
	mxFree(B);
	mxFree(R_);
	mxFree(temp);

	return eta;
	/*cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, Nsen_sq, Ntime, Nsrc_sq, 1, G, Nsen_sq, X, Nsrc_sq, 0, temp, Nsen_sq);*/
	/* cblas_dgemv(CblasColMajor, CblasNoTrans, Nsen_sq, Nsrc_sq, 1, G, Nsen_sq, X, 1, 0., temp, 1); */
}


