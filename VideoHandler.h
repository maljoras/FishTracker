
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "opencv2/video.hpp"

//#include <opencv2/opencv.hpp>
#include <iostream>
#include <string>   // for strings

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <glibmm/threads.h>
#include <glibmm/timer.h>
#include <glibmm/init.h>

#include "SaveVideoClass.h"
#include "FlyCapture2.h"

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
 * BackgroundThres Class. Similar to cv::BackgroundSubtractor
 */

class BackgroundThresholder
{
public:
    BackgroundThresholder();
    virtual ~BackgroundThresholder();

    void apply(cv::Mat frame, cv::Mat * bwimg);
    void clear(void);
    void setInverted(bool invertedif);
    int getNSkip(void);
    void setNSkip(int value);
    void setHistory(int value);
    int getHistory(void);
    double getThresScale(void);
    void setThresScale(double);
    

protected:
    int m_history;
    int m_nskip;
    int m_threstype;
    double m_adjustThresScale;    

private:
    float m_thres;
    cv::Mat m_meanImage;
    cv::Mat m_backImage;
    int m_istep;
    double m_adjustThresScaleCorrected;    
};


/**
 * VideoHandler Class. Loads video frames from files and detects contours in a threaded meanner
 */
class VideoHandler
{
public:

    /**
     * initializes video/camera
     */
    VideoHandler(const string fname,bool inKnnMethod);
    VideoHandler(int camIdx, const string fname,bool inKnnMethod);

    /** Destructor. */
    virtual ~VideoHandler();


    /**
     * computes the frame, foreground detection, and the contours in threaded manner. Returns -1 if falied 0 if success.
     */
    int step(vector<Segment> * pSegs, double * pTimeStamp, cv::Mat * pFrame);

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

    /**
     * stops the grabbing  (and deletd the capture/savevideo object etc )
     */
    void stop();
    
    /**
     * starts the grabbing/ writing etc (and builds the objects)
     */
    int start();


    /**
     * resets the background subtractor
     */

    void resetBkg();


    /**
     * get the current original frame (NOT GARANTUEED TO BE THE LATEST)
     */

    void getOFrame(cv::Mat * pFrame);

    /**
     * get the current bwimg frame (NOT GARANTUEED TO BE THE LATEST)
     */

    void getBWImg(cv::Mat * pBWImg);

    
private:
    
    void getSegment(Segment * segm, vector<cv::Point>inContour, cv::Mat inBwImg, cv::Mat inFrame,cv::Mat inOFrame);
    void plotFrame(cv::Mat pFrame,const string windowName);
    bool testValid(Segment * pSeg);
    void initPars();
    
    cv::Size m_featureSize; 
    bool m_camera;
    bool m_stopped;
    int m_camIdx;
    
    cv::Mat m_NextFrame;
    cv::Mat m_NextOFrame;
    cv::Mat m_NextBWImg;
    double m_NextTimeStamp;
    cv::Mat m_OFrame;
    cv::Mat m_Frame;
    cv::Mat m_BWImg;
    vector<Segment> m_Segments;
    double m_TimeStamp;
    
    // thread stuff
    void readNextFrameThread();
    void segmentThread();

    void reinitThreads();
    void startThreads();
    void deleteThreads();
    void waitThreads();
    
    Glib::Threads::Thread * m_nextFrameThread;
    Glib::Threads::Thread * m_segmentThread;
    Glib::Threads::Mutex m_NextFrameMutex;
    Glib::Threads::Mutex m_SegmentMutex;

    Glib::Threads::Cond m_emptySegmentCond;
    Glib::Threads::Cond m_emptyNextFrameCond;
    Glib::Threads::Cond m_availableSegmentCond;
    Glib::Threads::Cond m_availableNextFrameCond;

    bool m_keepSegmentThreadAlive;
    bool m_keepNextFrameThreadAlive;
    bool m_availableNextFrame;
    bool m_availableSegment;
    bool m_nextFrameThreadFinished;
    bool m_segmentThreadFinished;


    bool m_threadsAlive;

protected:
    cv::Ptr<cv::VideoCapture>  pVideoCapture;
    cv::Ptr<VideoSaver> pVideoSaver;
    cv::Ptr<cv::BackgroundSubtractorKNN> pBackgroundSubtractor;
    cv::Ptr<BackgroundThresholder> pBackgroundThresholder;

//    vector<vector<cv::Point> > m_NextContours;

    //   vector<Segment> m_NextSegments;

    
    vector<float> Scale;
    float Delta;
    string fname;


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

    bool computeSegments;
    
    bool resizeif;
    float resizescale;

    bool knnMethod;
    bool inverted;

};

