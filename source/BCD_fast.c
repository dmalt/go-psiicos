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

	/*mexPrintf("Hello!\n");*/
	openblas_set_num_threads(8);
	mwSize i;
	mwSize maxIter = 1e6;  /*One million */
	mwSize S = Nsrc_sq / 4; /*Need to check that Nsr_sq is divisible by 4*/
	double* rY_n = mxCalloc(Nsrc_sq * Ntime, sizeof(double));
	double* iY_n = mxCalloc(Nsrc_sq * Ntime, sizeof(double));
	double* Y_n = mxCalloc(Nsrc_sq * Ntime * 2, sizeof(double));
	/*mexPrintf("Hello!\n");*/
	cblas_dcopy(Nsrc_sq * Ntime, rY_p, 1, rY_n, 1);
	/*mexPrintf("Hello!\n");*/
	
	if(iY_p == NULL)
	{
		iY_p = mxCalloc(Nsrc_sq * Ntime, sizeof(double));
		for ( i = 0; i < Ntime *  Nsrc_sq; ++i)
		{
			/*mexPrintf("iY_p[i] = %f\n", iY_p[i]);*/
			iY_p[i] = 0.;
		}
	}
	cblas_dcopy(Nsrc_sq * Ntime, iY_p, 1, iY_n, 1);
/*	mexPrintf("Hello!\n");
*/
	double* rR = mxCalloc(Nsen_sq * Ntime, sizeof(double));
	double* iR = mxCalloc(Nsen_sq * Ntime, sizeof(double));
	double* R = mxCalloc(Nsen_sq * Ntime * 2, sizeof(double));
	double* M_ = mxCalloc(Nsen_sq * Ntime * 2, sizeof(double));
	cblas_dcopy(Nsen_sq * Ntime, rM_, 1, &M_[0], 1);
	cblas_dcopy(Nsen_sq * Ntime, iM_, 1, &M_[Nsen_sq * Ntime], 1);


	cblas_dcopy(Nsen_sq * Ntime, rM_, 1, rR, 1);
	cblas_dcopy(Nsen_sq * Ntime, iM_, 1, iR, 1);
/*	mexPrintf("Hello!\n")*/;

	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, Nsen_sq, Ntime, Nsrc_sq , -1.,  G, Nsen_sq, rY_p, Ntime,  1., rR, Nsen_sq);
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, Nsen_sq, Ntime, Nsrc_sq , -1.,  G, Nsen_sq, iY_p, Ntime,  1., iR, Nsen_sq);
	
/*	for (i = 0; i < Nsen_sq * Ntime; ++i)
	{
		mexPrintf("iR[i] = %f\n", iR[i]);
		mexPrintf("iM_ = %f\n", iM_[i]);
	}
*/
	mexPrintf("BCD\n");
	/*mexPrintf("randmax = %d, rand = %d\n", RAND_MAX, rand() % S);*/
	mwSize iter, src;
	double* G_s = (double*)mxMalloc(sizeof(double) * Nsen_sq * 4);	
	double* rY_n_s = (double*)mxMalloc(sizeof(double) * Ntime * 4);
	double* iY_n_s = (double*)mxMalloc(sizeof(double) * Ntime * 4);                   
	double* rY_p_s = (double*)mxMalloc(sizeof(double) * Ntime * 4);
	double* iY_p_s = (double*)mxMalloc(sizeof(double) * Ntime * 4);
	double* rdelta_Y = (double*)mxMalloc(sizeof(double) * Ntime * 4);
	double* idelta_Y = (double*)mxMalloc(sizeof(double) * Ntime * 4);
	double eta = 2.;
	double rY_n_s_norm = 0.,  iY_n_s_norm = 0., Y_n_s_norm = 0.;
	double scale = 0.;
	for (iter = 0; iter < maxIter; ++iter)
	{
		src = iter % S;	
		/*mexPrintf("src = %d", src);*/
		cblas_dcopy(4 * Ntime, &rY_n[4 * Ntime * src], 1, rY_n_s, 1);
		cblas_dcopy(4 * Ntime, &iY_n[4 * Ntime * src], 1, iY_n_s, 1);
		cblas_dcopy(4 * Ntime, &rY_n[4 * Ntime * src], 1, rY_p_s, 1);
		cblas_dcopy(4 * Ntime, &iY_n[4 * Ntime * src], 1, iY_p_s, 1);
		cblas_dcopy(4 * Nsen_sq, &G[4 * Nsen_sq * src], 1, G_s, 1);
		/*for (i = 0; i < 4 * Ntime; ++i)
		{
			mexPrintf("cblas_dcopy: rY_n_s = %f, rY_p_s = %f, rdelta_Y = %f\n", rY_n_s[i], rY_p_s[i], rdelta_Y[i]);
		}*/
		/*for (i = 0; i < 4 * Nsen_sq; ++i)
		{
			mexPrintf("cblas_dcopy: G_s = %f\n",G_s[i]);
		}*/
		cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, Ntime, 4, Nsen_sq, mu,  rR, Nsen_sq, G_s, Nsen_sq,  1., rY_n_s, Ntime);
		cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, Ntime, 4, Nsen_sq, mu,  iR, Nsen_sq, G_s, Nsen_sq,  1., iY_n_s, Ntime);
		/*for (i = 0; i < 4 * Ntime; ++i)
		{
			mexPrintf("cblas_dcopy: rY_n_s = %f\n", rY_n_s[i]);
		}*/
		rY_n_s_norm = cblas_dnrm2(Ntime * 4, rY_n_s, 1);
		iY_n_s_norm = cblas_dnrm2(Ntime * 4, iY_n_s, 1);
		Y_n_s_norm = sqrt(rY_n_s_norm * rY_n_s_norm + iY_n_s_norm * iY_n_s_norm);
		/*mexPrintf("norm = %f\n", Y_n_s_norm);*/
		scale = fmax(1 - mu * lambda / Y_n_s_norm, 0);
		/*mexPrintf("scale = %f\n", scale);*/

		cblas_dscal(Ntime * 4, scale, rY_n_s, 1);
		cblas_dscal(Ntime * 4, scale, iY_n_s, 1);

		cblas_dcopy(4 * Ntime, rY_n_s, 1, rdelta_Y, 1);
		cblas_dcopy(4 * Ntime, iY_n_s, 1, idelta_Y, 1);

		cblas_daxpy(Ntime * 4, -1., rY_p_s, 1, rdelta_Y, 1);
		cblas_daxpy(Ntime * 4, -1., iY_p_s, 1, idelta_Y, 1);
		/*for (i = 0; i < Nsen_sq * Ntime ; ++i)
		{
			mexPrintf("iR[i] = %f\n", iR[i]);
		}*/
		/*for ( i = 0; i < Ntime *  4; ++i)
			{
				mexPrintf("idelta_Y[i] = %f\n", idelta_Y[i]);
				mexPrintf("i_Y_n_s[i] = %f\n", iY_n_s[i]);
				mexPrintf("rdelta_Y[i] = %f\n", rdelta_Y[i]);
			}*/
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, Nsen_sq, Ntime, 4 , -1.,  G_s, Nsen_sq, rdelta_Y, Ntime,  1., rR, Nsen_sq);
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, Nsen_sq, Ntime, 4 , -1.,  G_s, Nsen_sq, idelta_Y, Ntime,  1., iR, Nsen_sq);

		cblas_dcopy(4 * Ntime, rY_n_s, 1, &rY_n[4 * Ntime * src], 1);
		cblas_dcopy(4 * Ntime, iY_n_s, 1, &iY_n[4 * Ntime * src], 1);
		/*for (i = 0; i < Nsen_sq * Ntime; ++i)
		{
			mexPrintf("rR[i] = %f\n", rR[i]);
		}*/
			/*for ( i = 0; i < Ntime *  Nsrc_sq; ++i)
			{
				mexPrintf("rY_n[i] = %f\n", rY_n[i]);
			}*/
		for (i = 0; i < Nsrc_sq; ++i)
		{		
			cblas_dcopy(Ntime, &rY_n[i * Ntime], 1, &Y_n[2 * i * Ntime], 1);
			cblas_dcopy(Ntime, &iY_n[i * Ntime], 1, &Y_n[(2 * i + 1) * Ntime], 1);
		}
		cblas_dcopy(Nsen_sq * Ntime, rR, 1, &R[0], 1);
		cblas_dcopy(Nsen_sq * Ntime, iR, 1, &R[Nsen_sq * Ntime], 1);
		/*for ( i = 0; i < 2* Ntime *  Nsrc_sq; ++i)
			{
				mexPrintf("Y_n[i] = %f\n", Y_n[i]);
			}*/
		/*for ( i = 0; i < Ntime *  Nsen_sq * 2; ++i)
			{
				mexPrintf("M_[i] = %f\n", M_[i]);
			}*/
	/*	for (i = 0; i < Nsen_sq * Ntime * 2; ++i)
		{
			mexPrintf("R[i] = %f\n", R[i]);
		}*/
			/*	for (i = 0; i < Nsrc_sq * Nsen_sq; ++i)
		{
			mexPrintf("G = %f\n",G[i]);
		}*/
		/*	mexPrintf("Ntime = %d, Nsen_sq = %d, S = %d, lambda = %f ", Ntime, Nsen_sq, S, lambda );*/
		if(!(iter % 20 * S))
		{
			eta  = CalcDualGap( Ntime * 2, Nsen_sq, S, lambda, M_, G, Y_n, R);
			
		}
		/*if(!(iter % 10000))
			mexPrintf("eta = %f\n", eta);*/
		if(eta < eps)
			break;
	}
	cblas_dcopy(Nsrc_sq * Ntime, rY_n, 1, rYout, 1);
	cblas_dcopy(Nsrc_sq * Ntime, iY_n, 1, iYout, 1);
	return iter;
} /* --- BlockCoorDescent --- */


double CalcDualGap(mwSize Ntime, mwSize Nsen_sq, mwSize S, double lambda,\
				 double * M, double * G, double * X, double * R)
{
	/*mexPrintf("Ntime = %d, Nsen_sq = %d, S = %d, lambda = %f ", Ntime, Nsen_sq, S, lambda );*/
	mwSize i;
	/*for ( i = 0; i < Ntime *  Nsen_sq; ++i)
			{
				mexPrintf("M[i] = %f\n", M[i]);
			}*/
	/*for (i = 0; i < S*4 * Nsen_sq; ++i)
	{
		mexPrintf("G = %f\n",G[i]);
	}*/
	/*for ( i = 0; i < Ntime *  S * 4; ++i)
	{
		mexPrintf("X[i] = %f\n", X[i]);
	}*/
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
		/*mexPrintf("%f\n", B[src]);*/
		/*for (i = 0; i < Ntime * 4; ++i)
		{
			mexPrintf("X_s[i] = %0.6f\n", X_s[i]);
		}
		mexPrintf("\n");*/
		sigma += cblas_dnrm2(4 * Ntime, X_s, 1); 
		/*mexPrintf("src = %d,sigma = %f\n", src, sigma);*/
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
	mxFree(B);
	mxFree(R_);
	mxFree(temp);

	return eta;
	/*cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, Nsen_sq, Ntime, Nsrc_sq, 1, G, Nsen_sq, X, Nsrc_sq, 0, temp, Nsen_sq);*/
	/* cblas_dgemv(CblasColMajor, CblasNoTrans, Nsen_sq, Nsrc_sq, 1, G, Nsen_sq, X, 1, 0., temp, 1); */
}


