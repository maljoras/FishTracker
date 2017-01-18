#include "SaveVideoClass.h"

#define sleep(x) std::this_thread::sleep_for(std::chrono::milliseconds(x)) 
#define M_TIMER_ELAPSED ((double)(( std::clock() - m_timer ) / (double) CLOCKS_PER_SEC))

VideoSaverFlyCapture::VideoSaverFlyCapture()
{ // base class constructor are automatically called
  m_isFlyCapture = false;
};

/****************************************************************************************/
VideoSaverFlyCapture::~VideoSaverFlyCapture()
{
  if (isInit()) {
    VideoSaver::close();
  }
}

int VideoSaverFlyCapture::stopCamera() {

  if (isFlyCapture()) {
    if (m_Camera.IsConnected()) {    
      FlyCapture2::Error error = m_Camera.Disconnect();
      if ( error != FlyCapture2::PGRERROR_OK ) {
	error.PrintErrorTrace();
	return -1;
      }
    }
    return 0;
  } else {
    return VideoSaver::stopCamera();
  }
}


/****************************************************************************************/

int VideoSaverFlyCapture::init(int camIdx){
  
  m_isFlyCapture = false;

  // check whether FlyCap camera exists
  FlyCapture2::PGRGuid guid;
  FlyCapture2::BusManager busMgr;
  FlyCapture2::Error error;

  error = busMgr.GetCameraFromIndex(camIdx, &guid );
  if (error == FlyCapture2::PGRERROR_OK) {
    return init(guid);
  }
  else { // try OpenCV
    cout << "Cannot find FlyCapture camera. Use OpenCV instead." << endl;
    return VideoSaver::init(camIdx);
  }  
}

/****************************************************************************************/
int VideoSaverFlyCapture::init(FlyCapture2::PGRGuid camIdx)
{
  FlyCapture2::Error error;
  FlyCapture2::CameraInfo camInfo;
  
  // Connect the camera
  error = m_Camera.Connect( &camIdx );
  if ( error != FlyCapture2::PGRERROR_OK )
    {
      std::cout << "Failed to connect to camera" << std::endl;     
      return -1;
    }
    
  // Get the camera info and print it out
  error = m_Camera.GetCameraInfo( &camInfo );
  if ( error != FlyCapture2::PGRERROR_OK )
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
  FlyCapture2::PropertyInfo propInfo;
  propInfo.type = FlyCapture2::FRAME_RATE;
  error = m_Camera.GetPropertyInfo( &propInfo );
  if (error != FlyCapture2::PGRERROR_OK)
    {
      error.PrintErrorTrace();
      return -1;
    }

  m_FrameRateToUse = 15.0f;
  if ( propInfo.present == true )
    {
      // Get the frame rate
      FlyCapture2::Property prop;
      prop.type = FlyCapture2::FRAME_RATE;
      error = m_Camera.GetProperty( &prop );
      if (error != FlyCapture2::PGRERROR_OK)
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
  FlyCapture2::Format7ImageSettings settings;
  unsigned int packetSize;
  float percentage;
  error = m_Camera.GetFormat7Configuration( &settings,&packetSize,&percentage );
  if ( error != FlyCapture2::PGRERROR_OK ) {
    error.PrintErrorTrace();
    return -1;
  }
  m_FrameSize =  cv::Size(settings.width,settings.height);

  settings.pixelFormat = FlyCapture2::PIXEL_FORMAT_RAW8;
  bool valid;
  FlyCapture2::Format7PacketInfo pinfo;
  error = m_Camera.ValidateFormat7Settings( &settings,&valid,&pinfo);
  if ( error != FlyCapture2::PGRERROR_OK ) {
    error.PrintErrorTrace();
    return -1;
  }

  if (!valid) {
    std::cout  << "Could not validate Format 7."  << std::endl;
    return -1;
  }

  error = m_Camera.SetFormat7Configuration( &settings,pinfo.recommendedBytesPerPacket);
  if ( error != FlyCapture2::PGRERROR_OK ) {
    error.PrintErrorTrace();
    return -1;
  }


  // set time stamping on
  FlyCapture2::EmbeddedImageInfo info;

  // Get configuration    
  error = m_Camera.GetEmbeddedImageInfo( &info );
  if ( error != FlyCapture2::PGRERROR_OK ) 
    {
      error.PrintErrorTrace();
      return -1;
    }

  info.timestamp.onOff = true;

  // Set configuration
  error = m_Camera.SetEmbeddedImageInfo( &info );
  if ( error != FlyCapture2::PGRERROR_OK ) 
    {
      error.PrintErrorTrace();
      return -1;
    }


  m_isFlyCapture = true;
  m_isRGB = true;
  return 0;
}



/****************************************************************************************/
int VideoSaverFlyCapture::startCapture() {

  if (!isFlyCapture()) {
    return VideoSaver::startCapture();
    
  } else {

  
    if ((isFinished()) && (isInit()) && (!isCapturing())) {
      // start thread to begin capture and populate Mat frame
      FlyCapture2::Error error = m_Camera.StartCapture();
      if ( error == FlyCapture2::PGRERROR_ISOCH_BANDWIDTH_EXCEEDED )
      {
	std::cout << "Bandwidth exceeded" << std::endl;     
	return -1;
      }
      else if ( error != FlyCapture2::PGRERROR_OK )
      {
	std::cout << "Failed to start image capture" << std::endl;     
	return -1;
      } 
    
    
      //start the grabbing thread
      m_KeepWritingAlive = false;  // not to be started
      m_WritingFinished = true;
      m_newFrameAvailable = false;
      std::cout <<  "Start video grabbing .." << std::endl;

      m_captureThread = new std::thread(&VideoSaverFlyCapture::captureThread,this);

      m_capturing = true;

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
    } 
  }
}

/****************************************************************************************/
bool VideoSaverFlyCapture::isInit() 
{
  if (isFlyCapture())
    return (m_Camera.IsConnected());
  else
    return VideoSaver::isInit();
}



/****************************************************************************************/
void VideoSaverFlyCapture::captureThread()
{
  if (!isFlyCapture()) {
    VideoSaver::captureThread();
  } else {

  
    m_GrabbingFinished = false;
    m_KeepThreadAlive = true;
    m_frameNumber = 0;
    m_newFrameAvailable = false;
    FlyCapture2::Image rgbImage;
    FlyCapture2::Image rawImage;
    FlyCapture2::Image rawImage2;

    m_timer= std::clock();
  
    while (m_KeepThreadAlive) 
    {

      FlyCapture2::Error error = m_Camera.RetrieveBuffer( &rawImage );
      if ( error != FlyCapture2::PGRERROR_OK )
      {
	error.PrintErrorTrace();
      }
      
      //get the time stamp
      const double localtimestamp = M_TIMER_ELAPSED;

      // convert to bgr
      rawImage2.DeepCopy(&rawImage); // not sure if really needed since we convert below...
      rawImage2.Convert(FlyCapture2::PIXEL_FORMAT_RGB, &rgbImage );

      // convert to Mat
      unsigned int rowBytes = (double) rgbImage.GetReceivedDataSize()/(double)rgbImage.GetRows();       

      // copy to frame variable and update times
      { std::unique_lock<std::mutex> lock(VideoSaver::m_FrameMutex); 

	m_Frame.release();
	m_Frame = cv::Mat(rgbImage.GetRows(), rgbImage.GetCols(), CV_8UC3, rgbImage.GetData(),rowBytes);
	// could this happen ?
	if (m_Frame.size().width==0) 
	  m_Frame = cv::Mat::zeros(m_FrameSize,CV_8UC3);

	m_TimeStamp =  localtimestamp;
	m_frameNumber++;
	m_newFrameAvailable = true;
	VideoSaver::m_newFrameAvailableCond.notify_one();
      }


    } 
    rawImage.ReleaseBuffer();
    rawImage2.ReleaseBuffer();
    rgbImage.ReleaseBuffer();

    m_newFrameAvailableCond.notify_one();
    std::cout << "Finished grabbing." << std::endl;    
    m_GrabbingFinished  = true;

  }
}






