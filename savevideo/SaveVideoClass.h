#include "FlyCapture2.h"

#include "FrameRateCounter.h"
#include <opencv2/core/core.hpp>

#include <opencv2/highgui/highgui.hpp>
//#include <opencv2/opencv.hpp>
#include <iostream>
#include <string>   // for strings
#include <sstream>
#include <fstream>   

// Include Boost headers for system time and threading
#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/thread.hpp"

using namespace std;
using namespace FlyCapture2;
using namespace boost::posix_time;
//using namespace cv;

//#include "Profiler.h"

/**
 * SaveVideo Class.  Saves and displays a video captured from a camaera.
 */ 
class VideoSaver  
{
public:
    /** Constructor. */
    VideoSaver();

    /** Destructor. */
    virtual ~VideoSaver();

    /**
     * Captures and writes the video to the file
     */ 
    int startCaptureAndWrite(const string fname, string codec);

    /**
     * returns the current frame
     */ 
    void getFrame(cv::Mat * pFrame);

    int close();

    int init(PGRGuid camIdx);

    /**
     * returns the current frame
     */ 
    double getWritingFPS();

    /**
     * Asks whether writing is finished
     */ 
    bool isFinished();

private:

   void _stopWriting();
   void _captureThread();
   void _captureAndWriteThread();

protected:

    FrameRateCounter m_FPSCounter;
    bool m_FirstImageGrabbed;
    bool m_KeepThreadAlive;
    bool m_KeepWritingAlive;
    bool m_WritingFinished;
    bool m_GrabbingFinished;
    boost::mutex m_FrameMutex;
    boost::mutex m_TimeMutex;

    FlyCapture2::Camera m_Camera;
    float m_FrameRateToUse;
    cv::Size m_FrameSize;

    cv::Mat m_Frame; 
    FlyCapture2::TimeStamp m_TimeStamp;

    std::fstream m_OutputFile;
    cv::VideoWriter m_Video;
};
