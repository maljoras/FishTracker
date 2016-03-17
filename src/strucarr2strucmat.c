/*
 * changes array of stuct into a struc of array-fields for vastly improved 
 * performance for large array sizes. 
 *
 * malte.rasch@bnu.edu.cn
 * 
 * inspried by the phonebook.c example of Matlab. 
 */
#include "mex.h"
#include "string.h"

#define MAXCHARS 80   /* max length of string contained in each 
                         field */

/* The gateway routine.  */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  const char **fnames;       /* pointers to field names */
  const size_t  *dims,*cdims;
  mxArray    *tmp, *fout;
  char       *pdata;
  int        ifield, jstruct, *classIDflags;
  int        NStructElems, nfields, ndim;
  size_t     cndim;
  int        s,*nofield,*nelements,*firstnonzero;
  int        i,alldim,snfields;
  size_t     *ccdims;
	  
  /* Check proper input and output */
  if (nrhs != 1)
    mexErrMsgTxt("One input required.");
  else if (nlhs > 1)
    mexErrMsgTxt("Too many output arguments.");
  else if (!mxIsStruct(prhs[0]))
    mexErrMsgTxt("Input must be a structure.");

  /* Get input arguments */
  nfields = mxGetNumberOfFields(prhs[0]);
  NStructElems = mxGetNumberOfElements(prhs[0]);

  /* Allocate memory  for storing classIDflags */
  classIDflags = mxCalloc(nfields, sizeof(int));

  nofield = mxCalloc(nfields, sizeof(int));
  nelements = mxCalloc(nfields, sizeof(int));
  firstnonzero = mxCalloc(nfields, sizeof(int));
  
  /* Check empty field, proper data type, and data type
      consistency; get classID for each field. */
  for (ifield = 0; ifield < nfields; ifield++) {
    nofield[ifield] = 0;
    nelements[ifield] = 0;
    firstnonzero[ifield] = 0;
    int s = 0;
    for (jstruct = 0; jstruct < NStructElems; jstruct++) {
      tmp = mxGetFieldByNumber(prhs[0], jstruct, ifield);
      int n;
      if (tmp == NULL)
	n = 0;
      else
	n  = mxGetNumberOfElements(tmp);

      if ((nelements[ifield]!=0) && (nelements[ifield]!=n) && (n!=0)){
	
	nofield[ifield] = 1;
	mexPrintf("%s%d\t%s%d\n", 
		  "FIELD:", ifield+1, "STRUCT INDEX :", jstruct+1);
	mexPrintf("Above field dimension mismatch: n=%d != %d = nelements\n",n,nelements[ifield]);
	break;
      }
      if (n!=0) {
	nelements[ifield] = n;
	s = s+1;
      } else
	continue;
      
	
      if (s == 1) {
        if ( /*(!mxIsChar(tmp) && (!(mxIsNumeric(tmp) || mxIsLogical(tmp)))) || */
            mxIsSparse(tmp)) { 
          mexPrintf("%s%d\t%s%d\n", 
              "FIELD:", ifield+1, "STRUCT INDEX :", jstruct+1);
          mexPrintf("Above field must have either "
              "string or numeric non-sparse data.\n");
	  nofield[ifield] = 1;
	  break;
        } 
        classIDflags[ifield] = mxGetClassID(tmp);
	firstnonzero[ifield] = jstruct; 
      } else {
        if (mxGetClassID(tmp) != classIDflags[ifield]) {
          mexPrintf("%s%d\t%s%d\n", 
              "FIELD:", ifield+1, "STRUCT INDEX :", jstruct+1);
          mexPrintf("Inconsistent data type in above field!\n");
	  nofield[ifield] = 1;
	  break;
        } 
        else if (!mxIsChar(tmp) && ((mxIsComplex(tmp))))  {
          mexPrintf("%s%d\t%s%d\n", 
              "FIELD:", ifield+1, "STRUCT INDEX :", jstruct+1);
          mexPrintf("Numeric data in above field "
              "must be scalar and noncomplex!\n");
	  nofield[ifield] = 1;
	  break;
        } 
      }
    }
  }
  /* Allocate memory  for storing pointers */
  fnames = mxCalloc(nfields, sizeof(*fnames));

  /* Get field name pointers */
  s = 0;
  for (ifield = 0; ifield < nfields; ifield++) {

    if (nofield[ifield]|| nelements[ifield]==0) {continue;};
    fnames[s] = mxGetFieldNameByNumber(prhs[0],ifield);
    s++;
  }

  /* Create a 1x1 struct matrix for output */
  plhs[0] = mxCreateStructMatrix(1, 1, s, fnames);
  mxFree(fnames);
  ndim = mxGetNumberOfDimensions(prhs[0]);
  dims = mxGetDimensions(prhs[0]);
  alldim = 1;
  for (i = 0; i < ndim; i++) {
    alldim *= dims[i];
  }


  s = 0;
  for (ifield = 0; ifield < nfields; ifield++) {

    if ((nofield[ifield]==1) || (nelements[ifield]==0)) {continue;};

    /* Create cell/numeric array */
    if ((classIDflags[ifield] == mxCHAR_CLASS) || (classIDflags[ifield] == mxSTRUCT_CLASS) ) {
      fout = mxCreateCellArray(ndim, dims);
  
    }  else {

      tmp = mxGetFieldByNumber(prhs[0],firstnonzero[ifield],ifield);
      cdims = mxGetDimensions(tmp);
      cndim = mxGetNumberOfDimensions(tmp);
      ccdims = mxCalloc(cndim+1, sizeof(size_t));      
      for (i = 0; i < cndim; i++) { ccdims[i] = cdims[i];}
      ccdims[cndim] = alldim;  /*cat to new dimension*/
      fout = mxCreateNumericArray(cndim+1, ccdims, 
				  classIDflags[ifield], mxREAL);
      pdata = mxGetData(fout); 
      mxFree(ccdims);
    }

    /* Copy data from input structure array */
    for (jstruct = 0; jstruct < NStructElems; jstruct++) {
      tmp = mxGetFieldByNumber(prhs[0],jstruct,ifield);
      if ((tmp!=NULL) && (mxIsChar(tmp) || (mxIsStruct(tmp))) ) {
        mxSetCell(fout, jstruct, mxDuplicateArray(tmp));
      }
      else {
        size_t     sizebuf;
        sizebuf = mxGetElementSize(fout)*nelements[ifield];
	
	if ((tmp==NULL) || (0==mxGetNumberOfElements(tmp))) {

	  if (classIDflags[ifield] == mxDOUBLE_CLASS) {
	    double * pddata = mxGetPr(fout); 
	    double nanVal = mxGetNaN();
	    for (i = 0; i<nelements[ifield]; i++) {
	      pddata[i +nelements[ifield]*jstruct] = nanVal;
	    }
	  } else if (classIDflags[ifield] == mxSINGLE_CLASS) {
	    float * pr = (float *) pdata;
	    float nanVal = mxGetNaN();
	    for (i = 0; i<nelements[ifield]; i++) {
	      pr[i +nelements[ifield]*jstruct] = nanVal;
	    }
	  }
	  else if (classIDflags[ifield] == mxINT32_CLASS) {
	    int * pr = (int *) pdata;
	    for (i = 0; i<nelements[ifield]; i++) {
	      pr[i +nelements[ifield]*jstruct] = 0;
	    }
	  }
	  else if (classIDflags[ifield] == mxUINT8_CLASS) {
	    unsigned char * pr = (unsigned char *) pdata;
	    for (i = 0; i<nelements[ifield]; i++) {
	      pr[i +nelements[ifield]*jstruct] = 0;
	    }
	  }
	  else {
	    mexPrintf("%d",classIDflags[ifield]);
	    mexErrMsgTxt("Unsupported data type:");
	  }
	
	} else {
	  memcpy(pdata, mxGetData(tmp),sizebuf);
	}
	pdata += sizebuf;
      }
    }

    /* Set each field in output structure */
    mxSetFieldByNumber(plhs[0], 0, s, fout);  
    s++;
  }
  mxFree(classIDflags);
  mxFree(nofield);
  mxFree(nelements);
  mxFree(firstnonzero);
  return;
}
