 
#include <thread>
#include <mutex>
#include <chrono>
#include <condition_variable>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "opencv2/video.hpp"

//#include <opencv2/opencv.hpp>
#include <iostream>
#include <string>   // for strings

#include <stdio.h>
#include <stdlib.h>
//#include <unistd.h>

#ifdef FLYCAPTURE
#include "SaveVideoClass.h"
#else
#include "SaveVideoClassBase.h"
#endif

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
    cv::Mat FilledImageFixedSize;
    cv::Mat FilledImageFixedSizeRotated;
    cv::Mat RotImage;
    cv::Mat RotFilledImage;
    cv::Mat IdentityFeature;
    cv::Mat IdentityFeatureRemap;
    cv::Mat CenterLine;    
    Scalar mback;
    
    cv::Rect Bbox;
    float Orientation;
    cv::Point2f Centroid;
    double MajorAxisLength;
    double MinorAxisLength;
    double Area;
    double bendingStdValue;
    double reversed;
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

    void apply(cv::Mat frame, cv::Mat * bwimg, cv::Mat * dframe);
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
    VideoHandler(int camIdx, const string fname,const string codec, bool inKnnMethod);

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
    void findBodyContours(cv::Mat inBwImg, vector<vector<cv::Point> > * newcontours);
    void plotFrame(cv::Mat pFrame,const string windowName);
    bool testValid(Segment * pSeg);
    void initPars();
    
    cv::Size m_featureSize; 


    bool m_camera; /* this property will be always false if flycapture
		    * is not present. However, one might still use
		    * capturing using the onecv capture feature. */
#ifdef FLYCAPTURE
    int m_camIdx;
#endif
    
    bool m_stopped;
    
    cv::Mat m_NextFrame;
    cv::Mat m_NextOFrame;
    cv::Mat m_NextBWImg;
    cv::Mat m_NextDFrame;
    double m_NextTimeStamp;
    cv::Mat m_OFrame;
    cv::Mat m_Frame;
    cv::Mat m_BWImg;
    vector<Segment> m_Segments;
    vector<Segment> m_NextSegments;
    double m_TimeStamp;
    
    // thread stuff
    void readNextFrameThread();
    void segmentThread();

    void reinitThreads();
    void startThreads();
    void deleteThreads();
    void waitThreads();
    
    std::thread * m_nextFrameThread; 
    std::thread * m_segmentThread;
    std::mutex m_NextFrameMutex;
    std::mutex m_SegmentMutex;

    std::condition_variable m_emptySegmentCond;
    std::condition_variable m_emptyNextFrameCond;
    std::condition_variable m_availableSegmentCond;
    std::condition_variable m_availableNextFrameCond;

    bool m_keepSegmentThreadAlive;
    bool m_keepNextFrameThreadAlive;
    bool m_availableNextFrame;
    bool m_availableSegment;
    bool m_nextFrameThreadFinished;
    bool m_segmentThreadFinished;


    bool m_threadsAlive;
    string m_codec;
    
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
    bool difffeature;
    int minWidth; // pixels for patch
    int minExtent;
    int maxExtent;
    int minArea;
    int maxArea;

    Size fixedSizeImage;
    int featurewidth;
    int featureheight;

    bool computeSegments;
    
    bool resizeif;
    float resizescale;

    bool knnMethod;
    bool inverted;

};

