
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
  double    *pdagIdx,*pdagPos,*pobjidx,*ptrace,*pifish,*pmobj;   	  
  mwSize  i,j,k,n,m,m2,t,mt,nhyp,nfish,nele1,nele2,tidx,stepsback;

  /* usage: otrace = backtrace_(dagIdx,dagPos,mt,ifish,mobj)*/
  
  /* Check proper input and output */
  if (nrhs != 6 ) {
    mexErrMsgTxt("6 inputs required. usage: backtrace_(dagIdx,dagPos,ifish,mobj,mt,stepsback)");
  } else if (nlhs != 2) {
    mexErrMsgTxt("Two outputs required.");
  } else if ((mxGetClassID(prhs[0])!=mxDOUBLE_CLASS) || (mxGetClassID(prhs[1])!=mxDOUBLE_CLASS)) {
    mexErrMsgTxt("Expect DOUBLE.");
  }

  
  /* Get input arguments */
  ndim1 = mxGetNumberOfDimensions(prhs[0]);
  dims1 = mxGetDimensions(prhs[0]);
  ndim2 = mxGetNumberOfDimensions(prhs[1]);
  dims2 = mxGetDimensions(prhs[1]);

  if ((ndim1 != 3) || (ndim2 != 3)) {
    mexErrMsgTxt("Expect 3D matrices as inputs.");		     
  } else if ((dims1[1]!= dims2[1]) || (dims1[2]!= dims2[2])){
    mexErrMsgTxt("Expect 2. and 3 dimensions to be identical.");		     
  }

  nele1 = mxGetNumberOfElements(prhs[2]);
  nele2 = mxGetNumberOfElements(prhs[3]);
  if (nele1 != nele2){
    mexErrMsgTxt("Expect IFISH and MOBJ to be of the same length");		     
  }

  mt = (mwSize) mxGetScalar(prhs[4]);
  mt = mt > dims1[2] ? dims1[2] : mt;

  stepsback = (mwSize) mxGetScalar(prhs[5]);

  if (stepsback>mt) {
    mexErrMsgTxt("stepsback cannot be larger than mt");
  }
  
  n = dims2[0]; /*pos dim*/
  
  dimsout = mxCalloc(3, sizeof(mwSize));      
  dimsout[0] = n; /*pos dim*/
  dimsout[1] = stepsback; /* */
  dimsout[2] = nele1;
  
  plhs[0] = mxCreateNumericArray(3, dimsout, mxDOUBLE_CLASS , mxREAL);
  ptrace = mxGetData(plhs[0]);
  mxFree(dimsout);

  dimsout = mxCalloc(2, sizeof(mwSize));      
  dimsout[0] = stepsback; /* */  
  dimsout[1] = nele1;
  plhs[1] = mxCreateNumericArray(2, dimsout, mxDOUBLE_CLASS , mxREAL);
  pobjidx = mxGetData(plhs[1]);
  mxFree(dimsout);

  pdagIdx = (double *)  mxGetData(prhs[0]);
  pdagPos = (double *)  mxGetData(prhs[1]);

  pifish = (double *) mxGetData(prhs[2]);
  pmobj = (double *)  mxGetData(prhs[3]);

  nfish = dims1[0];
  nhyp = dims1[1];

  /*populate array*/
  for (i=0;i<nele1;i++) {

    tidx = (stepsback-1) + i*stepsback;    

    pobjidx[tidx] = pmobj[i];
    for (j=0;j<n;j++) {
      ptrace[j + n*tidx] =  pdagPos[(mwSize) ( j + n*(pobjidx[tidx]-1 + nhyp*(mt-1)))];
    }

    for (t=2;t<=stepsback;t++) {

      tidx = (stepsback-t) + i*stepsback;    
      pobjidx[tidx] = pdagIdx[(mwSize) (pifish[i]-1 + nfish*(pobjidx[tidx+1]-1 +  nhyp*(mt-t+1)))];

      m = (mwSize) n*(pobjidx[tidx]-1+nhyp*(mt-t));
      for (j=0;j<n;j++) {
	ptrace[j + n*tidx] =  pdagPos[ j + m];
      }
    }
  }

}
