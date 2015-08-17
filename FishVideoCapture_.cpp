/**
 * @file VideoCapture_.cpp
 * @brief mex interface for VideoCapture_
 * @author Kota Yamaguchi
 * @date 2012
 */
#include "mexopencv.hpp"
using namespace std;
using namespace cv;

// Persistent objects

/// Last object id to allocate
int last_id = 0;
/// Object container
map<int,VideoCapture> obj_;

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
    if (nrhs<1 || nlhs>2)
        mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
    
    // Determine argument format between constructor or (id,method,...)
    vector<MxArray> rhs(prhs,prhs+nrhs);
    int id = 0;
    string method;
    if (nrhs==1) {
        // Constructor is called. Create a new object from argument
        obj_[++last_id] = (rhs[0].isChar()) ? 
            VideoCapture(rhs[0].toString()) : VideoCapture(rhs[0].toInt());
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
    VideoCapture& obj = obj_[id];
    if (method == "delete") {
        if (nrhs!=2 || nlhs!=0)
            mexErrMsgIdAndTxt("mexopencv:error","Output not assigned");
        obj_.erase(id);
    }
    else if (method == "open") {
        if (nrhs!=3 || nlhs!=1)
            mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
        bool b = (rhs[2].isChar()) ?
            obj.open(rhs[2].toString()) : obj.open(rhs[2].toInt());
        plhs[0] = MxArray(b);
    }
    else if (method == "isOpened") {
        if (nrhs!=2|| nlhs!=1)
            mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
        plhs[0] = MxArray(obj.isOpened());
    }
    else if (method == "release") {
        if (nrhs!=2 || nlhs!=0)
            mexErrMsgIdAndTxt("mexopencv:error","Output not assigned");
        obj.release();
    }
    else if (method == "grab") {
        if (nrhs!=2|| nlhs!=1)
            mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
        plhs[0] = MxArray(obj.grab());
    }
    else if (method == "retrieve") {
        if (nrhs!=2|| nlhs!=1)
            mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
        Mat frame;
        if (obj.retrieve(frame)) {
            if (frame.type()==CV_8UC3)
                cvtColor(frame,frame,cv::COLOR_BGR2RGB);
            plhs[0] = MxArray(frame);
        }
        else
            plhs[0] = MxArray(mxCreateNumericMatrix(0,0,mxUINT8_CLASS,mxREAL));
    }
    else if (method == "read") {
        if (nrhs!=2)
            mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
	Mat oframe;
        if (obj.read(oframe)) {
            if (oframe.type()==CV_8UC3)
                cvtColor(oframe,oframe,cv::COLOR_BGR2RGB);
            plhs[0] = MxArray(oframe);
	    if (nlhs==2) 
		plhs[1] = MxArray(oframe);
        }
        else
	{
            plhs[0] = MxArray(mxCreateNumericMatrix(0,0,mxUINT8_CLASS,mxREAL));
	    if (nlhs==2) 
		plhs[1] = MxArray(mxCreateNumericMatrix(0,0,mxUINT8_CLASS,mxREAL));
	}
    }
    else if (method == "readSingle") {
        if (nrhs!=2)
            mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
        Mat oframe;
        Mat frame;
        if (obj.read(oframe)) {
	    if (oframe.type()==CV_8UC3) {
	      cvtColor(oframe,oframe,cv::COLOR_BGR2RGB);
	      oframe.convertTo(frame, CV_32FC3, 1/255.0);
	      plhs[0] = MxArray(frame);
	      if (nlhs==2) 
		  plhs[1] = MxArray(oframe);

	    } else {
 	      mexErrMsgIdAndTxt("mexopencv:error","Expect CV_8UC3 color movie");
	    }
        }
        else
	{
            plhs[0] = MxArray(mxCreateNumericMatrix(0,0,mxSINGLE_CLASS,mxREAL));
	    if (nlhs==2) 
		plhs[1] = MxArray(mxCreateNumericMatrix(0,0,mxUINT8_CLASS,mxREAL));
	}
    }
    else if (method == "readGraySingle") {
        if (nrhs!=2)
            mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
        Mat frame;
        Mat oframe;
        if (obj.read(oframe)) {
	    if (oframe.type()==CV_8UC3) {
	      cvtColor(oframe,frame,cv::COLOR_BGR2GRAY);
	      frame.convertTo(frame, CV_32FC1, 1/255.0);
	      plhs[0] = MxArray(frame);
	      if (nlhs==2) 
		{
		  cvtColor(oframe,oframe,cv::COLOR_BGR2RGB);
		  plhs[1] = MxArray(oframe);
		}

	    } else {
 	      mexErrMsgIdAndTxt("mexopencv:error","Expect CV_8UC3 color movie");
	    }
        }
        else
	{
            plhs[0] = MxArray(mxCreateNumericMatrix(0,0,mxSINGLE_CLASS,mxREAL));
	    if (nlhs==2) 
		plhs[1] = MxArray(mxCreateNumericMatrix(0,0,mxUINT8_CLASS,mxREAL));
	}

    }
    else if (method == "readInvertedGraySingle") {
        if (nrhs!=2)
            mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
        Mat frame;
	Mat oframe;
        if (obj.read(oframe)) {
	    if (oframe.type()==CV_8UC3) {
	      cvtColor(oframe,frame,cv::COLOR_BGR2GRAY);
	      frame.convertTo(frame, CV_32FC1, 1/255.0);
	      plhs[0] = MxArray(1.0 - frame);
	      if (nlhs==2) 
		{
		  cvtColor(oframe,oframe,cv::COLOR_BGR2RGB);
		  plhs[1] = MxArray(oframe);
		}
	    } else {
 	      mexErrMsgIdAndTxt("mexopencv:error","Expect CV_8UC3 color movie");
	    }
        }
        else
	{
            plhs[0] = MxArray(mxCreateNumericMatrix(0,0,mxSINGLE_CLASS,mxREAL));
	    if (nlhs==2) 
		plhs[1] = MxArray(mxCreateNumericMatrix(0,0,mxUINT8_CLASS,mxREAL));
	}

    }

    else if (method == "readGray") {
        if (nrhs!=2)
            mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
        Mat frame;
	Mat oframe;
        if (obj.read(oframe)) {
	    if (oframe.type()==CV_8UC3) {
	      cvtColor(oframe,frame,cv::COLOR_BGR2GRAY);
	      plhs[0] = MxArray(frame);	
	      if (nlhs==2) 
		{
		  cvtColor(oframe,oframe,cv::COLOR_BGR2RGB);
		  plhs[1] = MxArray(oframe);
		}
	    } else {
 	      mexErrMsgIdAndTxt("mexopencv:error","Expect CV_8UC3 color movie");
	    }
        }
        else
	{
            plhs[0] = MxArray(mxCreateNumericMatrix(0,0,mxUINT8_CLASS,mxREAL));
	    if (nlhs==2) 
		plhs[1] = MxArray(mxCreateNumericMatrix(0,0,mxUINT8_CLASS,mxREAL));
	}
    }
    else if (method == "readInvertedGray") {
        if (nrhs!=2)
            mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
        Mat frame;
	Mat oframe;
        if (obj.read(oframe)) {
	    if (oframe.type()==CV_8UC3) {
	      cvtColor(oframe,frame,cv::COLOR_BGR2GRAY);
	      plhs[0] = MxArray(255-frame);
	      if (nlhs==2) 
		{
		  cvtColor(oframe,oframe,cv::COLOR_BGR2RGB);
		  plhs[1] = MxArray(oframe);
		}

	    } else {
 	      mexErrMsgIdAndTxt("mexopencv:error","Expect CV_8UC3 color movie");
	    }
        }
        else
	{       
	    plhs[0] = MxArray(mxCreateNumericMatrix(0,0,mxUINT8_CLASS,mxREAL));
	    if (nlhs==2) 
		plhs[1] = MxArray(mxCreateNumericMatrix(0,0,mxUINT8_CLASS,mxREAL));
	}

    }

    else if ((method == "readScaledS") || (method == "readScaledU")) {
        if (nrhs!=4)
            mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
        if ((mxGetNumberOfElements(rhs[2]) !=3) || (!mxIsClass(rhs[2],"double")))
             mexErrMsgIdAndTxt("mexopencv:error","Provide RGB single vector scale");
        if ((mxGetNumberOfElements(rhs[3]) !=3)|| (!mxIsClass(rhs[3],"double")))
             mexErrMsgIdAndTxt("mexopencv:error","Provide RGB single vector delta");


        Mat oframe;
	Mat channel[3];
	Mat output ;
	double  *scale, *delta; //!!!!!! NEEDS DOUBLE
	scale = mxGetPr(rhs[2]);
	delta = mxGetPr(rhs[3]);
	float cr= (float) scale[0]/255.;
	float cg= (float) scale[1]/255.;
	float cb= (float) scale[2]/255.;
	float z = (float) (delta[0] + delta[1] + delta[2]);

        if (obj.read(oframe)) {
	    if (oframe.type()==CV_8UC3) {

	      split(oframe, channel);
	      for (int ii=0; ii<3; ii++) {
		channel[ii].convertTo(channel[ii], CV_32FC1);
	      }
	      output =  cb*channel[0] + cg*channel[1] + cr*channel[2] + z;

	      // substract mean
	      output = output - cv::mean(output);

	      if  (method == "readScaledU") // rescale back to U
		{ 
		  output += 0.5; // between 0..1
		  output.convertTo(output, CV_8UC1, 255);	      
		}
	      plhs[0] = MxArray(output);
	      if (nlhs==2) 
		{
		  cvtColor(oframe,oframe,cv::COLOR_BGR2RGB);
		  plhs[1] = MxArray(oframe);
		}

	    } else {
 	      mexErrMsgIdAndTxt("mexopencv:error","Expect CV_8UC3 color movie");
	    }
        }
        else
	{	      
	    if  (method == "readScaledU") // rescale back to U
		plhs[0] = MxArray(mxCreateNumericMatrix(0,0,mxUINT8_CLASS,mxREAL));
	    else
		plhs[0] = MxArray(mxCreateNumericMatrix(0,0,mxSINGLE_CLASS,mxREAL));
	 
	    if (nlhs==2) 
		plhs[1] = MxArray(mxCreateNumericMatrix(0,0,mxUINT8_CLASS,mxREAL));
	}
    }

    else if (method == "get") {
        if (nrhs!=3)
            mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
        plhs[0] = MxArray(obj.get(CapProp[rhs[2].toString()]));
    }
    else if (method == "set") {
        if (nrhs!=4)
            mexErrMsgIdAndTxt("mexopencv:error","Wrong number of arguments");
        plhs[0] = MxArray(obj.set(CapProp[rhs[2].toString()],rhs[3].toDouble()));
    }
    else
        mexErrMsgIdAndTxt("mexopencv:error","Unrecognized operation");
}
