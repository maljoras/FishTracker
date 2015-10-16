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
map<int,VideoHandler *> obj_;

/** Capture Property map for option processing
 */
const ConstMap<std::string,int> CapProp = ConstMap<std::string,int>
    ("PosMsec",       cv::CAP_PROP_POS_MSEC)       //!< Current position of the video file in milliseconds or video capture timestamp.
    ("PosFrames",     cv::CAP_PROP_POS_FRAMES)     //!< 0-based index of the frame to be decoded/captured next.
    ("AVIRatio",      cv::CAP_PROP_POS_AVI_RATIO)  //!< Relative position of the video file: 0 - start of the film, 1 - end of the film.
    ("FrameWidth",    cv::CAP_PROP_FRAME_WIDTH)    //!< Width of the frames in the video stream.
    ("FrameHeight",   cv::CAP_PROP_FRAME_HEIGHT)   //!< Height of the frames in the video stream.
    ("FPS",           cv::CAP_PROP_FPS)            //!< Frame rate.
    ("FourCC",        cv::CAP_PROP_FOURCC)         //!< 4-character code of codec.
    ("FrameCount",    cv::CAP_PROP_FRAME_COUNT)    //!< Number of frames in the video file.
    ("Format",        cv::CAP_PROP_FORMAT)         //!< Format of the Mat objects returned by retrieve() .
    ("Mode",          cv::CAP_PROP_MODE)           //!< Backend-specific value indicating the current capture mode.
    ("Brightness",    cv::CAP_PROP_BRIGHTNESS)     //!< Brightness of the image (only for cameras).
    ("Contrast",      cv::CAP_PROP_CONTRAST)       //!< Contrast of the image (only for cameras).
    ("Saturation",    cv::CAP_PROP_SATURATION)     //!< Saturation of the image (only for cameras).
    ("Hue",           cv::CAP_PROP_HUE)            //!< Hue of the image (only for cameras).
    ("Gain",          cv::CAP_PROP_GAIN)           //!< Gain of the image (only for cameras).
    ("Exposure",      cv::CAP_PROP_EXPOSURE)       //!< Exposure (only for cameras).
    ("ConvertRGB",    cv::CAP_PROP_CONVERT_RGB)    //!< Boolean flags indicating whether images should be converted to RGB.
    //("WhiteBalance",cv::CAP_PROP_WHITE_BALANCE)  //!< Currently not supported
    ("Rectification", cv::CAP_PROP_RECTIFICATION)  //!< Rectification flag for stereo cameras (note: only supported by DC1394 v 2.x backend currently)
;

/// Field names for VideoHandler::Segments.
#define NSEGMENTFIELD 12
const char *segments_fields[NSEGMENTFIELD] = { "Bbox","Center","Orientation","Size","Image","FilledImage","FilledColImage","RotImage","RotFilledImage","RotFilledColImage","CenterLine","Thickness"};

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


    if (nrhs==1) {
        // Constructor is called. Create a new object from argument
        obj_[++last_id] =  new VideoHandler(rhs[0].toString());
        plhs[0] = MxArray(last_id);
        return;
    }
    else if (rhs[0].isNumeric() && rhs[0].numel()==1 && nrhs>1) {
        id = rhs[0].toInt();
        method = rhs[1].toString();
    }
    else
        mexErrMsgIdAndTxt("mexopencv:error","Invalid arguments");



    // Big operation switch
    VideoHandler& obj = *(obj_[id]);
    if (method == "delete") {
        if (nrhs!=2 || nlhs!=0)
            mexErrMsgIdAndTxt("mexopencv:error","Output not assigned");
        obj_.erase(id);
    }
    else if (method == "step") {
        if (nrhs!=2 || nlhs>4)
            mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments. Outputs: [seg,cont,bwimg,frame]");

        obj.step();

	if (nlhs>0) {
	    int n =0;
	    for (int i=0;i<obj.Segments.size();i++)
		if ((obj.Segments[i].Size.width>0) && (obj.Segments[i].Size.height>0))
		    n++;
	    mxArray * p = mxCreateStructMatrix(n,1,NSEGMENTFIELD, segments_fields);
	    if (!p)
		mexErrMsgIdAndTxt("mexopencv:error", "Allocation error");

	    int j = 0;
	    for (int i=0; i<obj.Segments.size(); i++) {
		if ((obj.Segments[i].Size.width>0)  && (obj.Segments[i].Size.height>0)) {
		    mxSetField(p,j,"Bbox",MxArray(obj.Segments[i].Bbox));
		    mxSetField(p,j,"Orientation",MxArray(obj.Segments[i].Orientation));
		    mxSetField(p,j,"Size",MxArray(obj.Segments[i].Size));
		    mxSetField(p,j,"Center",MxArray(obj.Segments[i].Center));
		    mxSetField(p,j,"Image",MxArray(obj.Segments[i].Image,mxLOGICAL_CLASS));
		    mxSetField(p,j,"FilledImage",MxArray(obj.Segments[i].FilledImage));
		    mxSetField(p,j,"FilledColImage",MxArray(obj.Segments[i].FilledColImage));
		    mxSetField(p,j,"RotImage",MxArray(obj.Segments[i].RotImage,mxLOGICAL_CLASS));
		    mxSetField(p,j,"RotFilledImage",MxArray(obj.Segments[i].RotFilledImage));
		    mxSetField(p,j,"RotFilledColImage",MxArray(obj.Segments[i].RotFilledColImage));
		    mxSetField(p,j,"CenterLine",MxArray(obj.Segments[i].CenterLine));
		    mxSetField(p,j,"Thickness",MxArray(obj.Segments[i].Thickness));
		    j++;
		}
	    }
	    plhs[0] = p;
	}
	if (nlhs>1) 
	    plhs[1] = MxArray(obj.Contours);
	if (nlhs>2) 
	    plhs[2] = MxArray(obj.BWImg,mxLOGICAL_CLASS);
	if (nlhs>3) 
	    plhs[3] = MxArray(obj.Frame);
	if (nlhs>4) 
	    plhs[4] = MxArray(obj.OFrame);
    }
    else if (method == "getframe") {
        if (nrhs!=2 || nlhs!=1)
            mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
	
	plhs[0] = MxArray(obj.Frame);
    }
    else if (method == "setMinContourSize") {
	if (nrhs!=3 || nlhs!=0)
	    mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
	obj.minContourSize = rhs[2].toInt();
    }
    else if (method == "getMinContourSize") {
	if (nrhs!=2 || nlhs!=1)
	    mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
	plhs[0] =MxArray(obj.minContourSize);
    }
    else if (method == "setScale") {
        if (nrhs<3 || nlhs!=0)
            mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
	if (!rhs[2].isNumeric()) {
	    if (rhs[2].toString()=="on")
		obj.scaled = true;
	    else 
		obj.scaled = false;
	} else if (nrhs=6) {
	    vector<float> v(3);
	    for (int i=0; i<3;i++) 
		v[i] = rhs[2+i].toFloat();
	    vector<float> delta(1);
	    delta[0] = rhs[5].toFloat();
	    obj.setScale(v,delta);
	} else
	    mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
    }
    else if (method == "isScaled") {
        if (nrhs!=2 || nlhs!=1)
            mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
	
	plhs[0] = MxArray(obj.scaled);
    }
    else if (method == "isPlotting") {
        if (nrhs!=2 || nlhs!=1)
            mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
	
	plhs[0] = MxArray(obj.plotif);
    }
    else if (method == "setTimePos") {
        if (nrhs!=3 || nlhs!=0)
            mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
	obj.setTime(rhs[2].toDouble());
    }
    else if (method == "getTimePos") {
        if (nrhs!=2 || nlhs!=1)
            mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
	plhs[0] = MxArray(obj.pVideoCapture->get(cv::CAP_PROP_POS_MSEC));
    }

    else if (method == "setPlotting") {
        if (nrhs!=3 || nlhs!=0)
            mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
	if (!rhs[2].isNumeric()) {
	    if (rhs[2].toString()=="on")
		obj.plotif = true;
	    else 
		obj.plotif = false;
	} else
	    obj.plotif = rhs[2].toBool();
    }
    else if ((method == "get") ||  (method == "capget")) {
        if (nrhs!=3)
            mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
        plhs[0] = MxArray(obj.pVideoCapture->get(CapProp[rhs[2].toString()]));
    }

    else if ((method == "set")|| (method == "capset")) {
        if (nrhs!=4)
            mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
        plhs[0] = MxArray(obj.pVideoCapture->set(CapProp[rhs[2].toString()],rhs[3].toDouble()));
    }
    else if (method == "bsget") {
        if (nrhs!=3 || nlhs>1)
            mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
        string prop(rhs[2].toString());
        if (prop == "DetectShadows")
            plhs[0] = MxArray(obj.pBackgroundSubtractor->getDetectShadows());
        else if (prop == "Dist2Threshold")
            plhs[0] = MxArray(obj.pBackgroundSubtractor->getDist2Threshold());
        else if (prop == "History")
            plhs[0] = MxArray(obj.pBackgroundSubtractor->getHistory());
        else if (prop == "kNNSamples")
            plhs[0] = MxArray(obj.pBackgroundSubtractor->getkNNSamples());
        else if (prop == "NSamples")
            plhs[0] = MxArray(obj.pBackgroundSubtractor->getNSamples());
        else if (prop == "ShadowThreshold")
            plhs[0] = MxArray(obj.pBackgroundSubtractor->getShadowThreshold());
        else if (prop == "ShadowValue")
            plhs[0] = MxArray(obj.pBackgroundSubtractor->getShadowValue());
        else
            mexErrMsgIdAndTxt("mexopencv:error","Unrecognized property");
    }
    else if (method == "bsset") {
        if (nrhs!=4 || nlhs!=0)
            mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
        string prop(rhs[2].toString());
        if (prop == "DetectShadows")
            obj.pBackgroundSubtractor->setDetectShadows(rhs[3].toBool());
        else if (prop == "Dist2Threshold")
            obj.pBackgroundSubtractor->setDist2Threshold(rhs[3].toDouble());
        else if (prop == "History")
            obj.pBackgroundSubtractor->setHistory(rhs[3].toInt());
        else if (prop == "kNNSamples")
            obj.pBackgroundSubtractor->setkNNSamples(rhs[3].toInt());
        else if (prop == "NSamples")
            obj.pBackgroundSubtractor->setNSamples(rhs[3].toInt());
        else if (prop == "ShadowThreshold")
            obj.pBackgroundSubtractor->setShadowThreshold(rhs[3].toDouble());
        else if (prop == "ShadowValue")
            obj.pBackgroundSubtractor->setShadowValue(rhs[3].toDouble());
        else
            mexErrMsgIdAndTxt("mexopencv:error","Unrecognized property");
    }
    else
        mexErrMsgIdAndTxt("mexopencv:error","Unrecognized operation");
}
