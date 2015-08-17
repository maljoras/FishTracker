

#include "SaveVideoClass.h"

using namespace std;
using namespace FlyCapture2;
using namespace boost::posix_time;
//using namespace cv;


VideoSaver::VideoSaver()
{
  m_KeepWritingAlive = false; // not yet in capture mode
  m_WritingFinished = false;
};


VideoSaver::~VideoSaver()
{ 
}

int VideoSaver::close()
{
  if (m_Camera.IsConnected()) {
    this->_stopWriting();
    
    while (!this->isFinished())
      usleep(1000);
    
    Error error = m_Camera.Disconnect();
    if ( error != PGRERROR_OK ) {
      error.PrintErrorTrace();
      return -1;
    }
  }
  return 0;
}

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
    cout  << "Could not validate Format 7."  << endl;
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


void VideoSaver::_stopWriting() 
{
  if (m_KeepWritingAlive) {
    m_KeepThreadAlive = false;
    m_KeepWritingAlive = false;
  }
} 

void VideoSaver::_captureThread()
{

  m_GrabbingFinished = false;
  m_FirstImageGrabbed = false;
  m_KeepThreadAlive = true;

  while (m_KeepThreadAlive) 
    {
      Image rawImage;
      Error error;
      error = m_Camera.RetrieveBuffer( &rawImage );
      if ( error != PGRERROR_OK )
	{
	  cout << "capture error\r";
	}
      
      //get the time stamp
      m_TimeMutex.lock();    
      m_TimeStamp = rawImage.GetTimeStamp();
      m_TimeMutex.unlock();    

      // convert to bgr
      Image rgbImage;
      rawImage.Convert(PIXEL_FORMAT_BGR, &rgbImage );

      // convert to Mat
      cv::Mat rawFrame;
      unsigned int rowBytes = (double) rgbImage.GetReceivedDataSize()/(double)rgbImage.GetRows();       
      rawFrame = cv::Mat(rgbImage.GetRows(), rgbImage.GetCols(), CV_8UC3, rgbImage.GetData(),rowBytes);

      // copy to frame variable
      m_FrameMutex.lock();    
      rawFrame.copyTo(m_Frame);
      m_FrameMutex.unlock();    

      m_FirstImageGrabbed  = true;
    } 
  m_GrabbingFinished  = true;
}

void VideoSaver::getFrame(cv::Mat * pFrame) 
{ 
  m_FrameMutex.lock();
  m_Frame.copyTo(*pFrame);
  m_FrameMutex.unlock();
}

double VideoSaver::getWritingFPS() 
{
  return m_FPSCounter.GetFrameRate();
}

bool VideoSaver::isFinished() 
{
  return m_WritingFinished;
}


int VideoSaver::startCaptureAndWrite(const string inFname, string codec)
{


  Error error;

  // open file stream for the tpoints
  string fname = string(inFname);
  string txtfname;
  txtfname = fname + ".txt";

  m_OutputFile.open(txtfname.c_str(), std::ios::out );

  if (!m_OutputFile.is_open()) 
    {
      cout  << "Could not open the output text for write: " << txtfname << endl;
      return -1;
    }

  //start the video stream
  m_Video = cv::VideoWriter(fname,CV_FOURCC(codec[0],codec[1],codec[2],codec[3]),m_FrameRateToUse, m_FrameSize ,true);

  if (!m_Video.isOpened())
    {
      cout  << "Could not open the output video for write: " << fname << endl;
      return -1;
    }

  // start thread to begin capture and populate Mat frame
  error = m_Camera.StartCapture();
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
  cout <<  "Start video grabbing .." << endl;
  boost::thread captureThread(&VideoSaver::_captureThread,this);

  // wait for startup
  while (!m_FirstImageGrabbed)
    usleep(1000);

  // start the writing thread
  cout <<  "Start video saving.." << endl;
  boost::thread captureAndWriteThread(&VideoSaver::_captureAndWriteThread,this);

  return 0;

}

void VideoSaver::_captureAndWriteThread()
{
  m_WritingFinished = false;
  Error error;

  // capture loop
  int frameNumber=0;
  m_FPSCounter.Reset();
  m_KeepWritingAlive = true;

  time_duration td;
  ptime thisTime,currentFrameTimestamp;
  int delayFound = 0;

  while(m_KeepWritingAlive)
    {
      //unsigned int micsec = ctime.microSeconds;
      currentFrameTimestamp =  microsec_clock::local_time();

      m_FrameMutex.lock();
      m_Video.write(m_Frame);
      m_FrameMutex.unlock();


      m_TimeMutex.lock();
      if (m_TimeStamp.microSeconds==0) // seems not to work
	m_OutputFile << frameNumber << "\t" <<  currentFrameTimestamp << endl;
      else 
	m_OutputFile << frameNumber << "\t" <<  m_TimeStamp.seconds << "\t" << m_TimeStamp.microSeconds << endl;

      m_TimeMutex.unlock();

      m_FPSCounter.NewFrame();
      frameNumber++;

      thisTime = microsec_clock::local_time();
      td = thisTime - currentFrameTimestamp;	
      delayFound = 1000./m_FrameRateToUse - td.total_milliseconds();
      if (delayFound>0)
	usleep(delayFound*1000);

    }

  while (!m_GrabbingFinished)
    usleep(1000);
    
  //close the files
  m_Video.release();
  m_OutputFile.close();

  // stop the camera
  error = m_Camera.StopCapture();
  if ( error != PGRERROR_OK )
    error.PrintErrorTrace();

  cout << "Finished writing" << endl;    
  m_WritingFinished = true;
};

