#include "mex.h"

void columnG_fast(mwSize p, double * G_small, double * w, mwSize Nsen, mwSize Nsrc, double * G_col);
void calc_violations(double * G_small, double * w, mwSize Nsen, mwSize Nsrc, double lambda, double * R, mwSize T, double * violations);
double calc_product_norm(double * col, double * R, mwSize Sn, mwSize T);


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
	mwSize s;
	double * column = mxGetPr(mxCreateDoubleMatrix(Sn,1,mxREAL));
	double norm = 0.;
	for (s = 0; s < 1e6; ++s)
	{
		columnG_fast(s+1, G_small, w, Nsen, Nsrc, column);
		norm = calc_product_norm(column, R, Sn, T);
		violations[s] = norm - lambda;
		/*if(!(s % 10000))
			mexPrintf("s = %d\n", s);*/
	}
/*	mxDestroyArray(column);
*/}

double calc_product_norm(double * col, double * R, mwSize Sn, mwSize T)
{
	double * product = mxGetPr(mxCreateDoubleMatrix(T,1,mxREAL));
	mwSize t, sn;
	double norm = 0;
	for (t = 0; t < 1; ++t)
	{
		product[t] = 0;
		for (sn = 0; sn < Sn; ++sn)
			product[t] += col[sn] * R[sn + Sn * t];
		norm += product[t] * product[t];
	}
	/*mxDestroyArray(product);*/
	norm = sqrt(norm);
	return norm;
}

void columnG_fast(mwSize p, double * G_small, double * w, mwSize Nsen, mwSize Nsrc, double * G_col)
{
	mwSize q = p % (Nsrc * Nsrc);
	mwSize half = p / (Nsrc * Nsrc); /* We have columns of two types: those which end with zeros and those which start with zeros. half determines the type*/
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