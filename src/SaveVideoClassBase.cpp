

#include "SaveVideoClassBase.h"

#define ELAPSED(TIMER) static_cast<std::chrono::duration<double,std::milli>>(std::high_resolution_clock::now() - TIMER).count()
#define TIMER_ELAPSED ELAPSED(timer)
#define TIMER_START timer = std::high_resolution_clock::now();
#define TIMER_INIT auto timer = std::high_resolution_clock::now();

#define sleep(x) std::this_thread::sleep_for(std::chrono::milliseconds(x)) 
#define M_TIMER_ELAPSED ELAPSED(m_timer) 



VideoSaver::VideoSaver()
{
  m_KeepWritingAlive = false; // not yet in capture mode
  m_KeepThreadAlive = false; // not yet in capture mode
  m_WritingFinished = true;
  m_GrabbingFinished = true;
  m_writing = false;
  m_capturing = false;
  m_newFrameAvailable = false;
  m_isRGB = false;
  m_writingFrameNumber = 0;
  m_frameNumber = 0;

};

/****************************************************************************************/
VideoSaver::~VideoSaver()
{
  close();
  sleep(500);
}

int VideoSaver::close()
{
  if (isInit()) {
    if (isCapturing()) {
      stopWriting();
    }
    stopCamera();
  }
  return 0;

}

/****************************************************************************************/
int VideoSaver::stopCamera() {

  if (isInit()) {    
    m_Capture.release();
  }
  return 0;
}

/****************************************************************************************/
int VideoSaver::init(int camIdx)
{

  // Connect the camera
  m_Capture.open(camIdx);
  if (!m_Capture.isOpened()){
    std::cout << "Failed to connect to camera" << std::endl; 
    return -1;
  }
	
  // get frame rate 
  double fps = (double) m_Capture.get(cv::CAP_PROP_FPS);
  
  if (!fps) {
    // property not supported. manually check frame rate
    TIMER_INIT;
    int nframes = 50;
    cv::Mat frame;
    for (int i=0; i<nframes;i++) {
      if (!m_Capture.read(frame)) {
        std::cout << "Cannot establish frame rate of camera" << std::endl;
        m_Capture.release(); 
        return -1;
        break; 
      }
    }
    fps = ((double) nframes)/TIMER_ELAPSED;
  }
 
  // Set the frame rate.
  // Note that the actual recording frame rate may be slower,
  // depending on the bus speed and disk writing speed.
  m_FrameRateToUse = fps;
	
  std::cout << "Using frame rate " <<  m_FrameRateToUse << std::endl;

  //get the width and height
  double height = (double) m_Capture.get(cv::CAP_PROP_FRAME_HEIGHT);
  double width  = (double) m_Capture.get(cv::CAP_PROP_FRAME_WIDTH);
  m_FrameSize =  cv::Size(width,height);

  m_isRGB = false; // done in BGR in OpenCV
  return 0;
}

/****************************************************************************************/
int VideoSaver::stopWriting() 
{

  m_KeepWritingAlive = false;
  m_KeepThreadAlive = false;
  m_newFrameAvailableCond.notify_one();
  
  while ((!m_GrabbingFinished) || (!m_WritingFinished))
    sleep(250);
  
  if ((m_capturing) && (m_captureThread!=NULL)) {
    //std::cout << "Joining capture thread..." ;    
    m_captureThread->join();
    delete(m_captureThread);
    m_captureThread = NULL;
    m_capturing = false;
    std::cout << "done" << std::endl;    
  }

  if ((m_writing) && (m_writingThread != NULL)) {
    //std::cout << "Joining writing thread..." ;    
    m_writingThread->join();
    delete(m_writingThread);
    m_writingThread = NULL;
    m_writing = false;
    std::cout << "done" << std::endl;    
  }

  return 0;
}

/****************************************************************************************/
cv::Size VideoSaver::getFrameSize() {
  return m_FrameSize;
}
/****************************************************************************************/
int VideoSaver::getCurrentFrameNumber() {
  if (!m_GrabbingFinished) {
    return m_frameNumber;
  } else {
    std::cout << "Warning: grabbing finished!!" << std::endl;
    return -1;
  }
}

/****************************************************************************************/
int VideoSaver::getLostFrameNumber() {
  return m_frameNumber-m_writingFrameNumber;
}

/****************************************************************************************/
int VideoSaver::getCurrentFrame(cv::Mat * pFrame ,double * pTimeStamp, int *pFrameNumber) 
{  // this grabs a frame from outside. To not interfere with the writing thread it will just grab the current frame,
   // NO WAITING TIME. Will re-grab the same frame if called again. 
  if (!m_GrabbingFinished) {
       
    {
      std::unique_lock<std::mutex> lock(m_FrameMutex);
      
      if (m_Frame.size().width==0) {
	      *pFrame = cv::Mat::zeros(m_FrameSize,CV_8UC3);
      }
      else {
	      m_Frame.copyTo(*pFrame);
      }
      *pTimeStamp = m_TimeStamp;
      *pFrameNumber = m_frameNumber;
    }
    return 0;
  }
  else {
    std::cout << "WARNING getFrame  Failed!" << std::endl;
    return -1;
  }
}



/****************************************************************************************/
double VideoSaver::getFPS() 
{
  return (double) m_FrameRateToUse;
}

/****************************************************************************************/
bool VideoSaver::isFinished() 
{
  return m_WritingFinished && m_GrabbingFinished;
}


/****************************************************************************************/
bool VideoSaver::isInit() 
{
  return (m_Capture.isOpened());
}

/****************************************************************************************/
bool VideoSaver::isCapturing() 
{
  return (isInit() && (m_capturing));
}


/****************************************************************************************/
int VideoSaver::startCapture() {

  if (isFinished() && (isInit()) && (!isCapturing())) {
    // start thread to begin capture and populate Mat frame
    
    //start the grabbing thread
    m_KeepWritingAlive = false;  // not to be started
    m_WritingFinished = true;
    m_newFrameAvailable = false;
    std::cout <<  "Start video grabbing .." << std::endl;
    
    m_captureThread = new std::thread(&VideoSaver::captureThread,this);

    m_capturing = true;
    // wait for startup
    sleep(500);
    waitForNewFrame();

    return 0;

  } else {
    if (isInit()) {
        std::cout << "Warning: capture not yet finished !" << std::endl;    
      } else {
        std::cout << "Warning: camera not available!" << std::endl;    
      }
    return -1;    
  };      
}

void VideoSaver::waitForNewFrame() {
  
  std::unique_lock<std::mutex> lock(m_FrameMutex);

  while (!m_newFrameAvailable) {
    m_newFrameAvailableCond.wait(lock);
  }

}


/****************************************************************************************/
int VideoSaver::startCaptureAndWrite(const string inFname, string codec)
{
  // start the capture 
  if (startCapture()!=0) 
    return -1;

  // open file stream for the tpoints
  string fname = string(inFname);
  string txtfname;
  txtfname = fname + ".txt";

  m_OutputFile.open(txtfname.c_str(), std::ios::out );

  if (!m_OutputFile.is_open()) 
    {
      std::cout  << "Could not open the output text for write: " << txtfname << std::endl;
      return -1;
    }

  //start the video stream
  m_Video = cv::VideoWriter(fname,CV_FOURCC(codec[0],codec[1],codec[2],codec[3]),m_FrameRateToUse, m_FrameSize ,true);

  if (!m_Video.isOpened())
    {
      std::cout  << "Could not open the output video for write: " << fname << std::endl;
      return -1;
    }

  
  // start the writing thread
  std::cout <<  "Start video saving.." << std::endl;
  m_writing = true;
  m_writingThread = new std::thread(&VideoSaver::captureAndWriteThread,this);
  
  return 0;

}



/****************************************************************************************/
void VideoSaver::captureThread()
{

  m_GrabbingFinished = false;
  m_KeepThreadAlive = true;
  m_frameNumber = 0;
  m_newFrameAvailable = false;

  m_timer= std::chrono::highsystem_clock::now();

	  
  while (m_KeepThreadAlive) {

    double localtimestamp = -1.;
    

    cv::Mat frame;
    if (m_Capture.grab()) {
      localtimestamp = M_TIMER_ELAPSED;  
      m_Capture.retrieve(frame);
    } else {
      std::cout<< "Error: a grabbing error occured" << std::endl;
      break;
    }
   
    frame.convertTo(frame,CV_8UC3); // necessary ?


    // copy to frame variable and update times
    {  std::unique_lock<std::mutex> lock(m_FrameMutex); 

      m_Frame.release();
      m_Frame = frame.clone();
      
      
      m_TimeStamp =localtimestamp; 
      m_frameNumber++;
      m_newFrameAvailable = true;
      m_newFrameAvailableCond.notify_one();
    }

  } 
  m_newFrameAvailableCond.notify_one();

  std::cout << "Finished grabbing." << std::endl;    
  m_GrabbingFinished  = true;
}


/****************************************************************************************/
void VideoSaver::captureAndWriteThread()
{
  m_WritingFinished = false;

  // capture loop
  m_writingFrameNumber=0;
  m_KeepWritingAlive = true;


  int delayFound = 0;
  int grabbedFrameNumber;
  cv::Mat frame;
  double localTimeStamp;
  int frameNumber;
    
  while(m_KeepWritingAlive) {
    
    //auto start_t =  M_TIMER_ELAPSED;

    waitForNewFrame();

    {
      std::unique_lock<std::mutex> lock(m_FrameMutex);
      if (isRGB()) {
	      cv::cvtColor(m_Frame,frame ,CV_RGB2BGR);
      } else {
	      m_Frame.copyTo(frame);
      }
      localTimeStamp = m_TimeStamp;
      grabbedFrameNumber = m_frameNumber;
      frameNumber  = m_writingFrameNumber++;
      m_newFrameAvailable = false; // need to reset !!
    }
      
    m_Video.write(frame); // slow, thus out of the lock

    m_OutputFile << frameNumber 
      << "\t" << grabbedFrameNumber <<"\t" 
      <<  std::fixed << std::setprecision(8) 
      << localTimeStamp << std::endl;


    //stop_t = M_TIMER_ELAPSED;
    //double dur_ms = start_t - stop_t;	
    //delayFound = (int) (1000./m_FrameRateToUse - dur_ms);
    //if (delayFound>0) {
	  //  sleep(delayFound);
    //}
  }

  while (!m_GrabbingFinished)
    sleep(200);
    
  //close the files
  m_Video.release();
  m_OutputFile.close();


  std::cout << "Finished writing." << std::endl;    
  m_WritingFinished = true;
};

