
#include <math.h>
#include "munkres.h"
#include <vector>
#include "mex.h"

/* The gateway routine.  */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  const mwSize  *dims1,*dims2,*dims3;
  mwSize *dimsout;
  mwSize   ndim1,ndim2,ndim3, rows,cols;
  double    *pdata1,*pdata2,*pdata3;   	  
  double    *pdataout;
  double   nonCostTrack,nonCostDetect;
  mwSize  i,j,n;

  constexpr auto infinity = std::numeric_limits<double>::infinity();
  
/* Check proper input and output */
  if ((nrhs != 2) && (nrhs != 3))  {
    mexErrMsgTxt("Two or three inputs required (cost,costUnassignedTracks,costUnassignedDetections).");
  } else if ((nlhs <1 ) || (nlhs > 4)){
    mexErrMsgTxt("Outputs required (assignments,unmatchedTracks,unmatchedDetections).");
  } else if ((mxGetClassID(prhs[0])!=mxDOUBLE_CLASS) || (mxGetClassID(prhs[1])!=mxDOUBLE_CLASS)) {
    mexErrMsgTxt("Expect DOUBLE.");
  }


 
  /* Get input arguments */
  ndim1 = mxGetNumberOfDimensions(prhs[0]);
  dims1 = mxGetDimensions(prhs[0]);

  
  if (ndim1 != 2) {
    mexErrMsgTxt("Expect cost matrix as input.");		     
  }
  
  dimsout = (mwSize *) mxCalloc(2, sizeof(mwSize));      

  n =dims1[1]+dims1[0];
  rows = dims1[0];
  cols = dims1[1];
  
  Matrix<double> matrix(n, n);
  pdata1 = (double *) mxGetData(prhs[0]); 
  /* make padded data */
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
      matrix(i,j) = infinity;
    }
  }

  for (i=0;i<rows;i++) {
    for (j=0;j<cols;j++) {
      matrix(i,j) = pdata1[i + j*rows];
    }
  }

  for (j=cols;j<n;j++) {
    for (i=rows;i<n;i++) {
      matrix(i,j) = 0;
    }
  }

  /* Non assignment costs*/
  ndim2 = mxGetNumberOfDimensions(prhs[1]);
  dims2 = mxGetDimensions(prhs[1]);

  if ((ndim2 != 2) || ((dims2[1]!=1) && (dims2[0]!=1))) {
    mexErrMsgTxt("Expect scalar or vector as costUnassignedTracks");		     
  }
  else if ((dims2[1]==dims2[0]) && ((dims2[0]==1) || (dims2[1]==1))){
    /* Scalar*/
    pdata2 = (double *) mxGetData(prhs[1]); 
    nonCostTrack = pdata2[0];

    if (infinity==nonCostTrack)
      mexErrMsgTxt("Cost cannot be inf.");		     
    for (i=0;i<rows;i++) {
      matrix(i,i+cols) = nonCostTrack;
    }
  }
  else if ((dims2[1]!=dims2[0]) && ((dims2[0]==rows) || (dims2[1]==rows) )) {
    pdata2 = (double *) mxGetData(prhs[1]); 

    for (i=0;i<rows;i++) {
      if (infinity==pdata2[i])
	mexErrMsgTxt("Cost cannot be inf.");		     
     
      matrix(i,i+cols) = pdata2[i];
    }
  } else {
    mexErrMsgTxt("Expect scalar or vector as  costUnassignedTracks");
  }
  

  if (nrhs==3) {
    if (mxGetClassID(prhs[2])!=mxDOUBLE_CLASS) {
      mexErrMsgTxt("Expect DOUBLE.");
    }

    ndim3 = mxGetNumberOfDimensions(prhs[2]);
    dims3 = mxGetDimensions(prhs[2]);
    if ((ndim3 != 2) || ((dims3[1]!=1) && (dims3[0]!=1))) {
      mexErrMsgTxt("Expect scalar or vector as costUnassignedDetections");		     
    }
    else if ((dims3[1]==dims3[0]) && ((dims3[0]==1) || (dims3[1]==1))){
      /* Scalar*/
      pdata3 = (double *) mxGetData(prhs[2]); 
      nonCostDetect =  pdata3[0];
      if (infinity==nonCostDetect)
	mexErrMsgTxt("Cost cannot be inf.");		     

      for (j=0;j<cols;j++) {
	matrix(rows+j,j) = nonCostDetect;
      }
    }  else if ((dims3[1]!=dims3[0]) && ((dims3[0]==cols) || (dims3[1]==cols) )) {
      /*VEctor cost*/
      pdata3 = (double *) mxGetData(prhs[2]); 
      for (j=0;j<cols;j++) {
	if (infinity== pdata3[j])
	  mexErrMsgTxt("Cost cannot be inf.");		     
	matrix(rows+j,j) = pdata3[j];
      }
      
    } else {
      mexErrMsgTxt("Expect scalar or vector as  costUnassignedDetections");
    }
  }
  else if ((dims2[1]==dims2[0]) && (dims2[0]==1)) {
    for (j=0;j<cols;j++) {
      matrix(rows+j,j) = nonCostTrack; /* take same cost as for tracks*/ 
    }
  } else {
    /* need scalar detection if vector in tracknoncost*/
    mexErrMsgTxt("Need costUnassignedDetections for vectorial costUnassignedTracks.");
  }
    
    

  

  /* apply munkres */
  Munkres<double> m;
  m.solve(matrix);
  
  /*get results*/
  std::vector<unsigned int> matchesRow(0);
  std::vector<unsigned int> matchesCol(0);
  std::vector<unsigned int> unmatchedTracks(0);
  std::vector<unsigned int> unmatchedDetections(0);


  for (j=0;j<n;j++) {
    for (i=0;i<n;i++) {
      if (matrix(i,j)==0){
	if ((i<rows) && (j<cols)) {
	  
	  matchesRow.push_back((unsigned int) i+1);
	  matchesCol.push_back((unsigned int) j+1);
	  
	} else if ((i>=rows) && (j<cols)) {
	  unmatchedDetections.push_back((unsigned int) j+1);
	} else if ((i<rows) && (j>=cols)) {
	  unmatchedTracks.push_back((unsigned int) i+1);
	}
      }
    }
  }


  
  /* create output matrix */
  dimsout[0] = matchesRow.size();
  dimsout[1] = 2;
  plhs[0] = mxCreateNumericArray(2, dimsout, mxDOUBLE_CLASS , mxREAL);
  pdataout = (double * ) mxGetData(plhs[0]);
  for (i=0;i<dimsout[0];i++) {
    pdataout[i ] =  (double) matchesRow[i];
    pdataout[i + matchesRow.size()] = (double) matchesCol[i];
  }

  if (nlhs>1)  {
    dimsout[0] = unmatchedTracks.size();
    dimsout[1] = 1;
    plhs[1] = mxCreateNumericArray(2, dimsout, mxDOUBLE_CLASS , mxREAL);
    pdataout = (double * ) mxGetData(plhs[1]);
    for (i=0;i<dimsout[0];i++) {
      pdataout[i ] = (double) unmatchedTracks[i];
    }
  }

  if (nlhs>2) {
    dimsout[0] = unmatchedDetections.size();
    dimsout[1] = 1;
    plhs[2] = mxCreateNumericArray(2, dimsout, mxDOUBLE_CLASS , mxREAL);
    pdataout =  (double * )mxGetData(plhs[2]);
    for (i=0;i<dimsout[0];i++) {
      pdataout[i ] = (double) unmatchedDetections[i];
    }
  }
  
  if (nlhs>3) {
    dimsout[0] = n;
    dimsout[1] = n;
    plhs[3] = mxCreateNumericArray(2, dimsout, mxDOUBLE_CLASS , mxREAL);
    pdataout =  (double * )mxGetData(plhs[3]);
    for (i=0;i<n;i++) {
      for (j=0;j<n;j++) {
	pdataout[i + n*j] = (double) matrix(i,j);
      }
    }
  }
  mxFree(dimsout);
}
