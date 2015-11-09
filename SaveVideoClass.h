#include "FlyCapture2.h"

#include "FrameRateCounter.h"
#include <opencv2/core/core.hpp>

#include <opencv2/highgui/highgui.hpp>
//#include <opencv2/opencv.hpp>
#include <iostream>
#include <string>   // for strings
#include <sstream>
#include <fstream>   


#include <glibmm/threads.h>
#include <glibmm/timer.h>
#include <glibmm/init.h>

// Include Boost headers for system time and threading

using namespace std;
using namespace FlyCapture2;
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
     * Captures without writing to video file
     */ 
    int startCapture();

    /**
     * makes the current frame available in "Frame"
     */ 
    int getFrame(cv::Mat * pFrame ,double * pTimeStamp, int *pFrameNumber);
    /**
     * returns the current framecounter
     */ 
    int getCurrentFrameNumber();
    
    int close();

    int init(PGRGuid camIdx);

    /**
     * returns the current frame rate of the video writing 
     */ 
    double getWritingFPS();

    /**
     * returns frame size 
     */
    cv::Size getFrameSize();

    /**
     * returns the theoretical FPS set (usually set by the camera)
     */
    double getFPS();

    /**
     * Asks whether writing is finished
     */ 
    bool isFinished();


private:

   void _stopWriting();
   void _captureThread();
   void _captureAndWriteThread();
   void waitForNewFrame();
   
protected:

    FrameRateCounter m_FPSCounter;
    Glib::Timer m_timer;  

    bool m_newFrameAvailable;
    bool m_KeepThreadAlive;
    bool m_KeepWritingAlive;
    bool m_WritingFinished;
    bool m_GrabbingFinished;
    bool m_writing;
    
    Glib::Threads::Mutex m_FrameMutex;
    Glib::Threads::Thread * m_captureThread;
    Glib::Threads::Thread * m_writingThread;
	
    FlyCapture2::Camera m_Camera;
    float m_FrameRateToUse;
    cv::Size m_FrameSize;
    
    
    cv::Mat m_Frame;
    FlyCapture2::TimeStamp m_TimeStamp;

    double m_LocalTimeStamp;

    int m_frameNumber;

    std::fstream m_OutputFile;
    cv::VideoWriter m_Video;
};
