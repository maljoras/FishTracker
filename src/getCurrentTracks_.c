
#include "mex.h"
#include "string.h"

#define MAXCHARS 80   /* max length of string contained in each field */

/*  the gateway routine.  */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    const char **fnamesin;
    const char **fnamesout;
    mxArray    *tmp;
    int        ifield, nfields;
    mwIndex    jstruct;
    mwSize     nStructElems;
    double     nfish;
    double    *saveFieldSegIf;
    const mxArray   *saveFieldsIn, *saveFieldsOut;


    
    /* INPUTS is : trackInfo = getCurrentTracks_(nfish,self.tracks,saveFieldIn,saveFieldOut,saveFieldSegIf)*/
    
    /* check proper input and output */
    if(nrhs!=5)
        mexErrMsgIdAndTxt( "FISHTRACKER:GetCurrentTracks:invalidNumInputs",
                "5 inputs required.");
    else if(nlhs != 1)
        mexErrMsgIdAndTxt( "FISHTRACKER:GetCurrentTracks:invalidNumOutputs",
                "1 output required.");
    else if(!mxIsDouble(prhs[0]))
        mexErrMsgIdAndTxt( "FISHTRACKER:GetCurrentTracks:inputNotDouble",
                "Nfish must be  double.");
    else if(!mxIsStruct(prhs[1]))
        mexErrMsgIdAndTxt( "FISHTRACKER:GetCurrentTracks:inputNotStruct",
                "Tracks must be a structure.");
    else if(!mxIsStruct(prhs[2]))
        mexErrMsgIdAndTxt( "FISHTRACKER:GetCurrentTracks:inputNotStruct",
                "Savefields must be structure.");
    else if(!mxIsStruct(prhs[3]))
        mexErrMsgIdAndTxt( "FISHTRACKER:GetCurrentTracks:inputNotStruc",
                "SavefieldsSub must be structure."); 
    else if(!mxIsDouble(prhs[4]))
        mexErrMsgIdAndTxt( "FISHTRACKER:GetCurrentTracks:inputNotDouble",
                "SavefieldsIf must be double."); 


    
    /* get input arguments */
    nfish = mxGetScalar(prhs[0]);
    saveFieldsIn = prhs[2];
    saveFieldsOut = prhs[3];

    nfields = mxGetNumberOfFields(saveFieldsIn);
    nStructElems = mxGetNumberOfElements(prhs[1]); /*tracks */

    if ((mxGetNumberOfElements(prhs[4])!=nfields) || (mxGetNumberOfFields(saveFieldsOut)!=nfields))
        mexErrMsgIdAndTxt( "FISHTRACKER:GetCurrentTracks:inputSizeMismatch",
                "Arrays must be of the same length."); 

    saveFieldSegIf = mxGetData(prhs[4]);

    /* create output struct*/
    fnamesin = mxCalloc(nfields, sizeof(*fnamesin));
    fnamesout = mxCalloc(nfields, sizeof(*fnamesout));
    for (ifield=0; ifield< nfields; ifield++){
      fnamesin[ifield] = mxGetFieldNameByNumber(saveFieldsIn,ifield);
      fnamesout[ifield] = mxGetFieldNameByNumber(saveFieldsOut,ifield);
    }
    plhs[0] = mxCreateStructMatrix(1, nfish, nfields, fnamesout);

    
    /* copy the data */
    for(ifield=0; ifield<nfields; ifield++) {
	
      for (jstruct=0; jstruct<nStructElems; jstruct++) {

	if (saveFieldSegIf[ifield]) {
	  tmp = mxGetField(mxGetField(prhs[1],jstruct, "segment" ),0, fnamesin[ifield]);
	}
	else {
	  tmp = mxGetField(prhs[1],jstruct,fnamesin[ifield]);
	}
	
	mxSetField(plhs[0], jstruct, fnamesout[ifield], mxDuplicateArray(tmp));
      }
    }
    mxFree((void *)fnamesin);   
    mxFree((void *)fnamesout);   
    return;
}
