
#include "mex.h"
#include "matrix.h"
#include <math.h>

/* The gateway routine.  */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  const mwSize  *dims1,*dims2;
  mwSize *dimsout;
  mxArray  *mout;
  mwSize   ndim1,ndim2;
  double    *pdata1,*pdata2;   	  
  double    *pdataout;
  mwSize  i,j,k,l,n,p,q,nc,nnc;

  /* Check proper input and output */
  if (nrhs != 2) {
    mexErrMsgTxt("Two inputs required.");
  } else if (nlhs != 1) {
    mexErrMsgTxt("One output required.");
  } else if ((mxGetClassID(prhs[0])!=mxDOUBLE_CLASS) || (mxGetClassID(prhs[1])!=mxDOUBLE_CLASS)) {
    mexErrMsgTxt("Expect DOUBLE.");
  }
  
  /* Get input arguments */
  ndim1 = mxGetNumberOfDimensions(prhs[0]);
  dims1 = mxGetDimensions(prhs[0]);
  ndim2 = mxGetNumberOfDimensions(prhs[1]);
  dims2 = mxGetDimensions(prhs[1]);

  if ((ndim1 > 3) || (ndim2 > 3)) {
    mexErrMsgTxt("Expect at most 3 dim as inputs.");		     
  } else if ((dims1[1]!= dims2[1]) || (dims1[0]!= dims2[0])){
    mexErrMsgTxt("Expect MxNxP MxNxQ data matrices as inputs.");		     
  }

  p = (ndim1==3)?dims1[2]:1;
  q = (ndim2==3)?dims2[2]:1;
  
  pdata1 = mxGetData(prhs[0]); 
  pdata2 = mxGetData(prhs[1]); 

  /* create output matrix */
  dimsout = mxCalloc(2, sizeof(mwSize));      
  dimsout[0] = p;
  dimsout[1] = q;
  plhs[0] = mxCreateNumericArray(2, dimsout, mxDOUBLE_CLASS , mxREAL);
  pdataout = mxGetData(plhs[0]);
    
  /* calculate the Euclidean distance*/
  n = dims1[1];
  nc = dims1[0]; /* number of centers*/
  nnc = n*nc;
  
  for (i=0;i<p;i++) {
    for (j=0;j<q;j++) {
      double y = 0;
      double yrev = 0;
      for (l=0;l<nc;l++) { /* n-center dim (first)*/
	double x = 0;
	double xrev = 0;
	for (k=0;k<n;k++) { /* euclid dim (second)*/
	  double d1 = pdata1[i*nnc + k*nc + l];
	  double d = d1 - pdata2[j*nnc + k*nc + l];
	  x += d*d;
	  d = d1 - pdata2[j*nnc + k*nc + (nc-l-1)]; /*reverse*/
	  xrev += d*d;
	} 
	y += sqrt(x);
	yrev += sqrt(xrev);
      }
      pdataout[j*p + i] = ((y<yrev)?y:yrev)/nc; /*min mean dist*/
    }
  }
  mxFree(dimsout);
}
