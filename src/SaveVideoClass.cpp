

#include "SaveVideoClass.h"

using namespace FlyCapture2;


#define sleep(x) std::this_thread::sleep_for(std::chrono::milliseconds(x)) 
#define TIMER_ELAPSED ( std::clock() - m_timer ) / (double) CLOCKS_PER_SEC

VideoSaver::VideoSaver()
{
  m_KeepWritingAlive = false; // not yet in capture mode
  m_KeepThreadAlive = false; // not yet in capture mode
  m_WritingFinished = true;
  m_GrabbingFinished = true;
  m_writing = false;
  m_newFrameAvailable = false;
  

};

/****************************************************************************************/
VideoSaver::~VideoSaver()
{
  close();
}

int VideoSaver::close()
{
  _stopWriting();

  if (m_Camera.IsConnected()) {    
    Error error = m_Camera.Disconnect();
    if ( error != PGRERROR_OK ) {
      error.PrintErrorTrace();
      return -1;
    }
    sleep(1000);
  }
  return 0;
}

/****************************************************************************************/
int VideoSaver::init(PGRGuid camIdx)
{
  Error error;
  CameraInfo camInfo;
  
  // Connect the camera
  error = m_Camera.Connect( &camIdx );
  if ( error != PGRERROR_OK )
    {
      std::cout << "Failed to connect to camera" << std::endl;     
      return -1;
    }
    
  // Get the camera info and print it out
  error = m_Camera.GetCameraInfo( &camInfo );
  if ( error != PGRERROR_OK )
    {
      std::cout << "Failed to get camera info from camera" << std::endl;     
      return -1;
    }
  std::cout << camInfo.vendorName << " "
	    << camInfo.modelName << " " 
	    << camInfo.serialNumber << std::endl;
	
	
  //-----------------
  // get frame rate
  // Check if the camera supports the FRAME_RATE property
  PropertyInfo propInfo;
  propInfo.type = FRAME_RATE;
  error = m_Camera.GetPropertyInfo( &propInfo );
  if (error != PGRERROR_OK)
    {
      error.PrintErrorTrace();
      return -1;
    }

  m_FrameRateToUse = 15.0f;
  if ( propInfo.present == true )
    {
      // Get the frame rate
      Property prop;
      prop.type = FRAME_RATE;
      error = m_Camera.GetProperty( &prop );
      if (error != PGRERROR_OK)
        {
	  error.PrintErrorTrace();
	  return -1;
        }
      else
	{
	  // Set the frame rate.
	  // Note that the actual recording frame rate may be slower,
	  // depending on the bus speed and disk writing speed.
	  m_FrameRateToUse = prop.absValue;
	}
    }
  printf("Using frame rate of %3.1f\n", m_FrameRateToUse);

  //get the width and height
  Format7ImageSettings settings;
  unsigned int packetSize;
  float percentage;
  error = m_Camera.GetFormat7Configuration( &settings,&packetSize,&percentage );
  if ( error != PGRERROR_OK ) {
    error.PrintErrorTrace();
    return -1;
  }
  m_FrameSize =  cv::Size(settings.width,settings.height);

  settings.pixelFormat = PIXEL_FORMAT_RAW8;
  bool valid;
  Format7PacketInfo pinfo;
  error = m_Camera.ValidateFormat7Settings( &settings,&valid,&pinfo);
  if ( error != PGRERROR_OK ) {
    error.PrintErrorTrace();
    return -1;
  }

  if (!valid) {
    std::cout  << "Could not validate Format 7."  << std::endl;
    return -1;
  }

  error = m_Camera.SetFormat7Configuration( &settings,pinfo.recommendedBytesPerPacket);
  if ( error != PGRERROR_OK ) {
    error.PrintErrorTrace();
    return -1;
  }


  // set time stamping on
  EmbeddedImageInfo info;

  // Get configuration    
  error = m_Camera.GetEmbeddedImageInfo( &info );
  if ( error != PGRERROR_OK ) 
    {
      error.PrintErrorTrace();
      return -1;
    }

  info.timestamp.onOff = true;

  // Set configuration
  error = m_Camera.SetEmbeddedImageInfo( &info );
  if ( error != PGRERROR_OK ) 
    {
      error.PrintErrorTrace();
      return -1;
    }

  return 0;
}

/****************************************************************************************/
void VideoSaver::_stopWriting() 
{

  m_KeepWritingAlive = false;
  m_KeepThreadAlive = false;

  while ((!m_GrabbingFinished) || (!m_WritingFinished))
    sleep(100);

  m_captureThread->join();
  delete(m_captureThread);
  m_captureThread = NULL;

  if (m_writing) {
    m_writingThread->join();
    delete(m_writingThread);
    m_writingThread = NULL;
  }
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
int VideoSaver::getFrame(cv::Mat * pFrame ,double * pTimeStamp, int *pFrameNumber) 
{
  if (!m_GrabbingFinished) {
    waitForNewFrame();
    
    {
      std::unique_lock<std::mutex> lock(m_FrameMutex);
      
      if (m_Frame.size().width==0) {
	      *pFrame = cv::Mat::zeros(m_FrameSize,CV_8UC3);
      }
      else {
	      m_Frame.copyTo(*pFrame);
      }
      *pTimeStamp = m_LocalTimeStamp;
      *pFrameNumber = m_frameNumber;
      m_newFrameAvailable = false;
    }
    return 0;
  }
  else {
    std::cout << "WARNING getFrame  Failed!" << std::endl;
    return -1;
  }
}

#ifdef COMPILE_AS_EXECUTABLE
/****************************************************************************************/
double VideoSaver::getWritingFPS() 
{
  return m_FPSCounter.GetFrameRate();
}
#endif

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
int VideoSaver::startCapture() {

  if (isFinished()) {
    // start thread to begin capture and populate Mat frame
    Error error = m_Camera.StartCapture();
    if ( error == PGRERROR_ISOCH_BANDWIDTH_EXCEEDED )
    {
      std::cout << "Bandwidth exceeded" << std::endl;     
      return -1;
    }
    else if ( error != PGRERROR_OK )
    {
      std::cout << "Failed to start image capture" << std::endl;     
      return -1;
    } 
    
    
    //start the grabbing thread
    m_KeepWritingAlive = false;  // not to be started
    m_WritingFinished = true;
    m_newFrameAvailable = false;
    std::cout <<  "Start video grabbing .." << std::endl;

    m_captureThread = new std::thread(&VideoSaver::_captureThread,this);

    // wait for startup
    sleep(500);
    waitForNewFrame();
    return 0;

  } else {
    std::cout << "Warning: capture not finished !" << std::endl;
    return -1;
  } 
  
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
  m_writingThread = new std::thread(&VideoSaver::_captureAndWriteThread,this);
  
  return 0;

}

/****************************************************************************************/
void VideoSaver::_captureThread()
{

  m_GrabbingFinished = false;
  m_KeepThreadAlive = true;
  m_frameNumber = 0;
  m_newFrameAvailable = false;
  Image rgbImage;
  Image rawImage;
  Image rawImage2;

  m_timer= std::clock();
  
  while (m_KeepThreadAlive) {

    Error error = m_Camera.RetrieveBuffer( &rawImage );
    if ( error != PGRERROR_OK ) {
	     error.PrintErrorTrace();
	  }
      
    //get the time stamp
    const double localtimestamp = TIMER_ELAPSED;

    // convert to bgr
    rawImage2.DeepCopy(&rawImage); // not sure if really needed since we convert below...
    rawImage2.Convert(PIXEL_FORMAT_RGB, &rgbImage );

    // convert to Mat
    unsigned int rowBytes = (unsigned int) ((double) rgbImage.GetReceivedDataSize())/(double)rgbImage.GetRows();       

    // copy to frame variable and update times
    {
	     std::unique_lock<std::mutex> lock(m_FrameMutex);
	     m_Frame.release();
	     m_Frame = cv::Mat(rgbImage.GetRows(), rgbImage.GetCols(), CV_8UC3, rgbImage.GetData(),rowBytes);
	     // could this happen ?
	     if (m_Frame.size().width==0) 
	       m_Frame = cv::Mat::zeros(m_FrameSize,CV_8UC3);

	     //rawFrame.copyTo(m_Frame); 
	     m_TimeStamp = rawImage.GetTimeStamp();
	     m_frameNumber++;
	     m_LocalTimeStamp = localtimestamp;
	     m_newFrameAvailable = true;
       m_newFrameAvailableCond.notify_one();
    }

  } 
  rawImage.ReleaseBuffer();
  rawImage2.ReleaseBuffer();
  rgbImage.ReleaseBuffer();
  m_newFrameAvailableCond.notify_one();
    
  // stop the camera
  Error error = m_Camera.StopCapture();
  if ( error != PGRERROR_OK )
    error.PrintErrorTrace();
  
  m_GrabbingFinished  = true;
}


/****************************************************************************************/
void VideoSaver::_captureAndWriteThread()
{
  m_WritingFinished = false;
  Error error;

  // capture loop
  int frameNumber=0;
#ifdef COMPILE_AS_EXECUTABLE
  m_FPSCounter.Reset();
#endif
  m_KeepWritingAlive = true;


  int delayFound = 0;
  int grabbedFrameNumber;
  cv::Mat frame;
  FlyCapture2::TimeStamp timeStamp;
  double localTimeStamp;
  
  while(m_KeepWritingAlive) {
    
    //unsigned int micsec = ctime.microSeconds;
    const double currentTime =  TIMER_ELAPSED;

    {
	    std::unique_lock<std::mutex> lock(m_FrameMutex);
	    cv::cvtColor(m_Frame,frame,CV_RGB2BGR);	
	    localTimeStamp = m_LocalTimeStamp;
	    timeStamp = m_TimeStamp;
	    grabbedFrameNumber = m_frameNumber;
    }
      
    m_Video.write(frame); // slow, thus out of the lock


    if (timeStamp.microSeconds==0)  { // seems not to work
	     m_OutputFile << frameNumber << "\t" << grabbedFrameNumber <<"\t" <<  std::fixed << std::setprecision(5) << localTimeStamp << std::endl;
    }
    else { 
	     m_OutputFile << frameNumber << "\t" << grabbedFrameNumber <<"\t" <<  std::fixed << std::setprecision(5) << timeStamp.seconds << "\t" << timeStamp.microSeconds << std::endl;
    }

    frameNumber++; // this is the writing number 

#ifdef COMPILE_AS_EXECUTABLE
    m_FPSCounter.NewFrame();
#endif
    const double thisTime = TIMER_ELAPSED;
    const double seconds = thisTime - currentTime;	
    delayFound = static_cast<int> (1000./m_FrameRateToUse - seconds*1000.);
    if (delayFound>0) {
	     sleep(delayFound*1000);
    }
  }

  while (!m_GrabbingFinished)
    sleep(1000);
    
  //close the files
  m_Video.release();
  m_OutputFile.close();


  std::cout << "Finished writing" << std::endl;    
  m_WritingFinished = true;
};

