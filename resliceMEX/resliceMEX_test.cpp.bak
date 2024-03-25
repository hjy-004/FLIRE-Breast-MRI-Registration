#ifdef STAND_ALONE
// For StandAlone App.
#include "../volReslice/resliceMEX_external.h"
#else
// For MEXing
#include <mex.h>
#endif

#include "getvxlval.h"

#ifdef STAND_ALONE
void resliceMEX(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[])
#else
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
#endif
{
	double *pinfo_i, *psi, *pimgsi, *pso, *pmat;
	int	interpm, padding;
	int size_i[3], size_o[3];
	int i,j,k;
	pimgsi	= mxGetPr(prhs[0]);
	psi		= mxGetPr(prhs[1]);
	for (i=0;i<3;i++)
		size_i[i]=(int)psi[i];
	pmat	= mxGetPr(prhs[2]);
	pso		= mxGetPr(prhs[3]);
	for (i=0;i<3;i++)
		size_o[i]=(int)pso[i];
	interpm = (int)mxGetScalar(prhs[4]);
    	padding = (int)mxGetScalar(prhs[5]);

	plhs[0]=mxCreateNumericArray(3, size_o, mxDOUBLE_CLASS, mxREAL);
	double *pimgs_o=mxGetPr(plhs[0]);    
        
	int indo;
	double vxl[4];
	for (k=1;k<=size_o[2];k++)
	for (i=1;i<=size_o[0];i++)
	for (j=1;j<=size_o[1];j++)
	{
		indo=index(size_o,i,j,k);
		vxl[0]=i;vxl[1]=j;vxl[2]=k;vxl[3]=1.;
		getvxlval(vxl, size_i, pimgsi, pmat, interpm, padding, pimgs_o[indo]);
    }
	return;

}
