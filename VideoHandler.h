
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "opencv2/video.hpp"

//#include <opencv2/opencv.hpp>
#include <iostream>
#include <string>   // for strings

#include <stdio.h>
#include <stdlib.h>

// Include Boost headers for system time and threading
#include "boost/thread.hpp"

using namespace std;
using namespace cv;

class Segment
{
public:

    /** Constructor. */
    Segment();

    /** Destructor. */
    virtual ~Segment();

    cv::Mat Image;
    cv::Mat FilledImage;
    cv::Mat RotImage;
    cv::Mat RotFilledImage;
    cv::Mat FishFeature;
    cv::Mat FishFeatureRemap;
    cv::Mat CenterLine;    
    Scalar mback;
    
    cv::Rect Bbox;
    float Orientation;
    cv::Point2f Centroid;
    double MajorAxisLength;
    double MinorAxisLength;
    double Area;
    double bendingStdValue;
    cv::Size2f Size;

    vector<double> Thickness;
};



/**
 * VideoHandler Class. Loads video frames from files and detects contours in a threaded meanner
 */
class VideoHandler
{
public:

    /**
     * initializes video
     */
    VideoHandler(const string fname);

    /** Destructor. */
    virtual ~VideoHandler();


    /**
     * computes the frame, foreground detection, and the contours in threaded manner
     */
    void step();

    /**
     * sets the scaling of RGB.
     */
    void setScale(vector<float> scale);

    /**
     * gets the scaling of RGB.
     */
    vector<float> getScale();

    /**
     * sets the time for video to start
     */
    void setTime(double msec);

    /**
     * sets properties
     */
    int set(const string prop, double value);
    
    /**
     * gets the properties
     */
    double get(const string prop);
    
    
    //properties
    cv::Ptr<cv::BackgroundSubtractorKNN> pBackgroundSubtractor;
    cv::Ptr<cv::VideoCapture>  pVideoCapture;

    cv::Mat OFrame;
    cv::Mat Frame;
    cv::Mat BWImg;

    vector<Segment> Segments;
    int ngoodmsk;
    vector<bool> goodmsk;

   
private:
    void stopAndInitialize();
    void startThread();
    int _readNextFrameThread();
    void getSegment(Segment * segm, vector<cv::Point>inContour, cv::Mat inBwImg, cv::Mat inFrame,cv::Mat inOFrame);
    void plotCurrentFrame();
    void makeGoodMsk();

    cv::Size featureSize; 

protected:

    cv::Mat m_NextFrame;
    cv::Mat m_NextOFrame;
    cv::Mat m_NextBWImg;
//    vector<vector<cv::Point> > m_NextContours;
    vector<vector<cv::Point> > Contours;
    //   vector<Segment> m_NextSegments;

    vector<float> m_Scale;
    float m_Delta;

    boost::mutex m_FrameMutex;

    bool scaled;
    bool plotif;
    int nprobe;

    bool colorfeature;
    int minWidth; // pixels for patch
    int minExtent;
    int maxExtent;
    int minArea;
    int maxArea;

    int featurewidth;
    int featureheight;
};

