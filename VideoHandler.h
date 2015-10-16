
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
    cv::Mat FilledColImage;
    cv::Mat RotImage;
    cv::Mat RotFilledImage;
    cv::Mat RotFilledColImage;

    cv::Rect Bbox;
    float Orientation;
    cv::Point2f Center;
    cv::Size2f Size;
    vector<cv::Point> CenterLine;
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
    void setScale(vector<float> scale, vector<float> delta);

    /**
     * sets the time for video to start
     */
    void setTime(double msec);


    //properties
    cv::Ptr<cv::BackgroundSubtractorKNN> pBackgroundSubtractor;
    cv::Ptr<cv::VideoCapture>  pVideoCapture;



    cv::Mat OFrame;
    cv::Mat Frame;
    cv::Mat BWImg;
    vector<vector<cv::Point> > Contours;
    vector<Segment> Segments;

    bool scaled;
    bool plotif;
    int minContourSize ;
    int nprobe;

private:
    int _readNextFrameThread();
    void getSegment(Segment * segm, vector<cv::Point>inContour, cv::Mat inBwImg, cv::Mat inFrame,cv::Mat inOFrame);
    void plotCurrentFrame();

protected:

    cv::Mat m_NextFrame;
    cv::Mat m_NextOFrame;
    cv::Mat m_NextBWImg;
    vector<vector<cv::Point> > m_NextContours;
    vector<Segment> m_NextSegments;

    vector<float> m_Scale;
    float m_Delta;

    boost::mutex m_FrameMutex;
};

