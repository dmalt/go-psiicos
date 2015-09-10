#include "mex.h"
#include <cblas.h>

void columnG_fast(mwSize p, double * G_small, double w, mwSize Nsen, mwSize Nsrc, double * G_col);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if(nrhs != 3)
    	mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "Three inputs required.");
    mwSize p = mxGetScalar(prhs[0]);
 	double * G_small = mxGetPr(prhs[1]);
 	mwSize Nsen = mxGetM(prhs[1]);
 	mwSize Nsrc = mxGetN(prhs[1]);
 	double w = mxGetScalar(prhs[2]);
 	double * G_col;
 	plhs[0] = mxCreateDoubleMatrix(Nsen * Nsen,1,mxREAL);
 	G_col = mxGetPr(plhs[0]);
 	columnG_fast(p, G_small, w, Nsen, Nsrc, G_col);
}

void columnG_fast(mwSize p, double * G_small, double w, mwSize Nsen, mwSize Nsrc, double * G_col)
{
	mwSize r = (p - 1) % 4;
	mwSize s = (p - 1) / 4;
	mwSize Nsites = Nsrc / 2; 
	mwSize is = s % Nsites;
	
	mwSize js = (s - is) / Nsites;
	mwSize i, j;
	mwSize k, l;

	mwSize tcount;
	bool cblasErrFlag = 0;
	if(r == 0) 
	{
		i = 2 * is;
		j = 2 * js;
	}
	else if(r == 1)
	{
		i = 2 * is + 1;
		j = 2 * js;
	}
	else if(r == 2)
	{
		i = 2 * is;
		j = 2 * js + 1;
	}
	else if(r == 3)
	{
		i = 2 * is + 1;
		j = 2 * js + 1;
	}
	/*for (k = 0; k < Nsen; ++k)
		for (l = 0; l < Nsen; ++l)
			G_col[k + Nsen * l] = G_small[k + Nsen * i] * G_small[l + Nsen * j] * w;*/
	/*double * G_col_ = mxCalloc(Nsen * Nsen, sizeof(double));*/
	cblas_dger(CblasColMajor, /* you’re using row-major storage */
           Nsen,           /* the matrix X has dx1 rows ...  */
           Nsen,           /*  ... and dx2 columns.          */
           w,           /* scale factor to apply to x1x2' */
           &G_small[Nsen * i],
           1,             /* stride between elements of x1. */
           &G_small[Nsen * j],
           1,             /* stride between elements of x2. */
           G_col,
           Nsen);          /* leading dimension of matrix */

    /*for (tcount = 0; tcount < Nsen * Nsen; ++tcount)
    {
    	G_col_[tcount] -=G_col[tcount];
    	if(G_col_[tcount] > 1e-6)
    	{
    		cblasErrFlag = 1;
    		mexPrintf("Tes failed\n");
    		break;
    	}
    }
    if(!cblasErrFlag){
    	mexPrintf("Test OK!\n");
    }*/
}
