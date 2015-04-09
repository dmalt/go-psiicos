#include "mex.h"
#include <cblas.h>
#include <omp.h>
#include <stdio.h>

void columnG_fast(mwSize p, double * G_small, double * w, mwSize Nsen, mwSize Nsrc, double * G_col);
void calc_violations(double * G_small, double * w, mwSize Nsen, mwSize Nsrc, double lambda, double * R, mwSize T, double * violations);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if(nrhs != 4)
    	mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "Four inputs required.");
    double * G_small = mxGetPr(prhs[0]);
 	mwSize Nsen = mxGetM(prhs[0]);
 	mwSize Nsrc = mxGetN(prhs[0]);
 	double * w = mxGetPr(prhs[1]);
 	double lambda = mxGetScalar(prhs[2]);
 	double * R = mxGetPr(prhs[3]);
 	mwSize T = mxGetN(prhs[3]);

 	double * violations;
 	plhs[0] = mxCreateDoubleMatrix(Nsrc * Nsrc * 2, 1, mxREAL);
	violations = mxGetPr(plhs[0]);
	calc_violations(G_small, w, Nsen, Nsrc, lambda, R, T, violations);
}

void calc_violations(double * G_small, double * w, mwSize Nsen, mwSize Nsrc, double lambda, double * R, mwSize T, double * violations)
{
	mwSize Sr = Nsrc * Nsrc * 2;
	mwSize Sn = Nsen * Nsen * 2; /* Number of elements in a column */
	mwSize s, i, nthreads;
	double * column = mxGetPr(mxCreateDoubleMatrix(Sn, 1, mxREAL));
	double * result = mxGetPr(mxCreateDoubleMatrix(T, 1, mxREAL));
	double * RtimeMat = mxGetPr(mxCreateDoubleMatrix(10, 1, mxREAL));
	double norm = 0.;
	/*mexPrintf("T = %d\n", T);*/
/*	#pragma omp parallel firstprivate(column, result, norm) shared(s, Sr, Sn, T, R, lambda, violations) private(nthreads) num_threads(8)
	{*/	
		/*nthreads = omp_get_num_threads();
    	printf("Number of threads = %d\n", nthreads);*/
		#pragma omp parallel for private(s) num_threads(8)/*ordered */
		for (s = 0; s < 10; ++s)
		{
			columnG_fast(s+1, G_small, w, Nsen, Nsrc, column);
			cblas_dgemv(101, 111, T, Sn, 1, R, Sn, column, 1, 0., result, 1);
			/*for(i = 0; i < T; i++)
				mexPrintf("%f,", result[i]);*/
			/*mexPrintf("\n");*/
			norm = cblas_dnrm2(T, result, 1);
			violations[s] = norm - lambda;
			printf("s = %d\n", s);
		}
	/*}*//*	mxDestroyArray(column);
*/}



void columnG_fast(mwSize p, double * G_small, double * w, mwSize Nsen, mwSize Nsrc, double * G_col)
{
	mwSize q = p % (Nsrc * Nsrc);
	mwSize half = (p - 1) / (Nsrc * Nsrc); /* We have columns of two types: those which end with zeros and those which start with zeros. half determines the type*/
	mwSize left = 0, right = 1;
	if(!q)
		q = Nsrc * Nsrc;
	mwSize i = q % Nsrc;
	if(!i)
		i = Nsrc;
	i --;
	mwSize j = (q - i) / Nsrc;
	double w_p = w[p-1];
	mwSize k, l;
	for (k = 0; k < Nsen; ++k)
		for (l = 0; l < Nsen; ++l)
		{
			if(half == left)
			{
				G_col[k + Nsen * l] = G_small[k + Nsen * i] * G_small[l + Nsen * j] * w_p;
				G_col[Nsen * Nsen + k + Nsen * l] = 0.;
			}
			else if(half == right)
			{
				G_col[k + Nsen * l] = 0.;
				G_col[Nsen * Nsen + k + Nsen * l] = G_small[k + Nsen * i] * G_small[l + Nsen * j] * w_p;
			}
		}
}
