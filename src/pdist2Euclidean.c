
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
  mwSize  i,j,k,n;

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

  if ((ndim1 != 2) || (ndim2 != 2)) {
    mexErrMsgTxt("Expect matrices as inputs.");		     
  } else if (dims1[1]!= dims2[1]){
    mexErrMsgTxt("Expect MXxN MYxN data matrices as inputs.");		     
  }

  pdata1 = mxGetData(prhs[0]); 
  pdata2 = mxGetData(prhs[1]); 

  /* create output matrix */
  dimsout = mxCalloc(2, sizeof(mwSize));      
  dimsout[0] = dims1[0];
  dimsout[1] = dims2[0];
  plhs[0] = mxCreateNumericArray(2, dimsout, mxDOUBLE_CLASS , mxREAL);
  pdataout = mxGetData(plhs[0]);
    
  /* calculate the Euclidean distance*/
  n = dims1[1];
  for (i=0;i<dimsout[0];i++) {
    for (j=0;j<dimsout[1];j++) {
      double x = 0;
      for (k=0;k<n;k++) {
	double d = pdata1[i + k*dims1[0]] - pdata2[j + k*dims2[0]];
	x += d*d;
      } 
      pdataout[j*dimsout[0] + i] = sqrt(x);
    }
  }
  mxFree(dimsout);
}
