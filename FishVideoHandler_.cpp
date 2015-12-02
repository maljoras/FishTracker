/**
 * adapted from:
 * @brief mex interface for VideoCapture_
 * @author Kota Yamaguchi
 * @date 2012
 */

#include "mexopencv.hpp"
#include "VideoHandler.h"
using namespace std;
using namespace cv;


// Persistent objects

/// Last object id to allocate
int last_id = 0;
/// Object container
map<int,Ptr<VideoHandler> > obj_;


/// Field names for VideoHandler::Segments.
#define NSEGMENTFIELD 18
const char *segments_fields[NSEGMENTFIELD] = { "BoundingBox","Centroid","Area","Orientation","Size","MinorAxisLength","MajorAxisLength","Image","FilledImage","RotImage","RotFilledImage","FishFeatureC","FishFeatureCRemap","CenterLine","Thickness","mback","bendingStdValue","MSERregions"};

/**
 * Main entry called from Matlab
 * @param nlhs number of left-hand-side arguments
 * @param plhs pointers to mxArrays in the left-hand-side
 * @param nrhs number of right-hand-side arguments
 * @param prhs pointers to mxArrays in the right-hand-side
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    if (nrhs<1 || nlhs>4)
        mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
    
    // Determine argument format between constructor or (id,method,...)
    vector<MxArray> rhs(prhs,prhs+nrhs);
    int id = 0;
    string method;

   
    if ((nrhs==4) && (rhs[0].isChar() && (rhs[0].toString()=="camera"))){
        // CApture Constructor is called. Create a new object from argument
	if (!obj_.empty())
	    mexErrMsgIdAndTxt("mexopencv:error","Only one camera instance supported");
	obj_[++last_id] =  new VideoHandler(rhs[1].toInt(),rhs[2].toString(),rhs[3].toBool());
        plhs[0] = MxArray(last_id);
        return;

    } else if ((nrhs==2) && rhs[0].isChar()) {
        // Constructor is called. Create a new object from argument
	obj_[++last_id] =  new VideoHandler(rhs[0].toString(),rhs[1].toBool());
        plhs[0] = MxArray(last_id);
        return;
    }
    else if (rhs[0].isNumeric() && rhs[0].numel()==1 && nrhs>1) {

	id = rhs[0].toInt();
	if (obj_.find(id)==obj_.end())
	    return; // no object availalbe

        method = rhs[1].toString();
    }
    else
        mexErrMsgIdAndTxt("mexopencv:error","Invalid arguments");


    // Big operation switch
    VideoHandler& obj = *(obj_[id]);
    if (method == "delete") {
        if (nrhs!=2 || nlhs!=0)
            mexErrMsgIdAndTxt("mexopencv:error","Output not assigned");
	obj_[id].release();
        obj_.erase(id);
    }
    else if (method == "step") {
        if (nrhs!=2 || nlhs>4)
            mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments. Outputs: [seg,bwimg,frame,oframe]");

	obj.step();

	if (nlhs>0) {
	    if (obj.Segments.size()) {
		mxArray * p = mxCreateStructMatrix(obj.Segments.size(),1,NSEGMENTFIELD, segments_fields);
		if (!p)
		    mexErrMsgIdAndTxt("mexopencv:error", "Allocation error");
		

		for (int i=0; i<obj.Segments.size(); i++) {
		    mxSetField(p,i,"Image",MxArray(obj.Segments[i].Image,mxLOGICAL_CLASS));
		    mxSetField(p,i,"FilledImage",MxArray(obj.Segments[i].FilledImage));
		    mxSetField(p,i,"RotImage",MxArray(obj.Segments[i].RotImage,mxLOGICAL_CLASS));
		    mxSetField(p,i,"RotFilledImage",MxArray(obj.Segments[i].RotFilledImage));
		    mxSetField(p,i,"FishFeatureC",MxArray(obj.Segments[i].FishFeature));
		    mxSetField(p,i,"FishFeatureCRemap",MxArray(obj.Segments[i].FishFeatureRemap));
		    mxSetField(p,i,"BoundingBox",MxArray(obj.Segments[i].Bbox));
		    mxSetField(p,i,"Orientation",MxArray(-obj.Segments[i].Orientation + 90.));
		    mxSetField(p,i,"Size",MxArray(obj.Segments[i].Size));
		    mxSetField(p,i,"MajorAxisLength",MxArray(obj.Segments[i].MajorAxisLength));
		    mxSetField(p,i,"MinorAxisLength",MxArray(obj.Segments[i].MinorAxisLength));
		    mxSetField(p,i,"Centroid",MxArray(obj.Segments[i].Centroid));

		    mxSetField(p,i,"CenterLine",MxArray(obj.Segments[i].CenterLine));
		    mxSetField(p,i,"Thickness",MxArray(obj.Segments[i].Thickness));
		    mxSetField(p,i,"Area",MxArray(obj.Segments[i].Area));
		    mxSetField(p,i,"mback",MxArray(obj.Segments[i].mback));
		    mxSetField(p,i,"bendingStdValue",MxArray(obj.Segments[i].bendingStdValue));
		    //mxSetField(p,i,"MSERregions",MxArray(NULL));
		}
		plhs[0] = p;
	    }
	    else {
		plhs[0] = mxCreateDoubleMatrix( 0, 0, mxREAL );
	    }
	}

	if (nlhs>1) 
	    plhs[1] = MxArray(obj.BWImg,mxLOGICAL_CLASS);
	if (nlhs>2) 
	    plhs[2] = MxArray(obj.Frame);
	if (nlhs>3) 
	    plhs[3] = MxArray(obj.OFrame);


    }
    else if (method == "resetBkg") {
        if (nrhs!=2 || nlhs!=0)
            mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
	
	obj.resetBkg();
    }
    else if (method == "getFrame") {
        if (nrhs!=2 || nlhs!=1)
            mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
	cv::Mat frame;
	obj.getOFrame(&frame);
	plhs[0] = MxArray(frame);
    }
    else if (method == "setScale") {
        if (nrhs!=5 || nlhs!=0)
            mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
	
	vector<float> v(3);
	for (int i=0; i<3;i++) 
	    v[i] = rhs[2+i].toFloat();
	obj.setScale(v);
    }
    else if (method == "getScale") {
        if (nrhs!=2 || nlhs!=1)
            mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
	plhs[0] = MxArray(obj.getScale());
    }
    else if (method == "start") {
        if (nrhs!=2 || nlhs!=0)
            mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
	if (obj.start()!=0)
	    mexErrMsgIdAndTxt("mexopencv:error"," Could not start the VideoHandler");
	
    }
    else if (method == "stop") {
        if (nrhs!=2 || nlhs!=0)
            mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
	obj.stop();
    }

    // else if (method == "capget") {
    //     if (nrhs!=3)
    //         mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
    //     plhs[0] = MxArray(obj.pVideoCapture->get(CapProp[rhs[2].toString()]));
    // }

    // else if (method == "capset") {
    //     if (nrhs!=4)
    //         mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
    //     plhs[0] = MxArray(obj.pVideoCapture->set(CapProp[rhs[2].toString()],rhs[3].toDouble()));
    // }
    else if (method == "get") { // calls directly yhe properties ov VideoHandler
        if (nrhs!=3)
            mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
	plhs[0] = MxArray(obj.get(rhs[2].toString()));
    }
    else if (method == "set") { // calls directly yhe properties ov VideoHandler
        if (nrhs!=4)
            mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
	if (obj.set(rhs[2].toString(),rhs[3].toDouble())!=0) {
	    mexErrMsgIdAndTxt("mexopencv:error","Got an error returned");
	}
    }
    else
        mexErrMsgIdAndTxt("mexopencv:error","Unrecognized operation");
}
