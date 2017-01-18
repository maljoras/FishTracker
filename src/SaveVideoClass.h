#include "FlyCapture2.h"
#include "SaveVideoClassBase.h"

using namespace std;


/**
 * SaveVideo Class.  Saves and displays a video captured from a camaera.
 */ 
class VideoSaverFlyCapture: public VideoSaver  
{
public:
  /** Constructor. */
  VideoSaverFlyCapture();

  /** Destructor. */
  virtual ~VideoSaverFlyCapture();

  /**
   */
  bool isFlyCapture() {
    return m_isFlyCapture;
  };

	
  /**
   * Captures and writes the video to the file
   */ 
  virtual int startCapture();

  virtual int init(int camIdx);
  
  int init(FlyCapture2::PGRGuid camIdx);

  virtual bool isInit();
  
  using VideoSaver::close;
  
private:

  FlyCapture2::Camera m_Camera;
  bool m_isFlyCapture;
  
protected:

  virtual int stopCamera();
  virtual void captureThread();	

};
