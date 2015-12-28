

#include "VideoHandler.h"
//#define DEBUG
//#define PLOTSEGMENTS

#define MAXCONTOUR 100
#define MAXVALIDCONTOUR 25
#define MAXOTSU 25
#define MINBACKTHRESSTEP 20

using namespace std;
using namespace cv;


Segment::Segment() {}
Segment::~Segment() {};


VideoHandler::VideoHandler(const string inFname,bool inKnnMethod)
{
  m_camera = false;
  fname = inFname;
  m_stopped = true;

  knnMethod = inKnnMethod;

  m_nextFrameThreadFinished = true;
  m_segmentThreadFinished = true;
  m_availableNextFrame = false;
  m_availableSegment = false;
  m_threadsAlive = false;

  Glib::init();
  initPars();
};

/****************************************************************************************/
VideoHandler::VideoHandler(int camIdx, const string inFname, bool inKnnMethod)
{
  m_camera = true;
  m_stopped = true;
  m_camIdx = camIdx;
  fname = inFname;

  knnMethod = inKnnMethod;

  m_nextFrameThreadFinished = true;
  m_segmentThreadFinished = true;
  m_availableNextFrame = false;
  m_availableSegment = false;
  m_threadsAlive = false;
  
  Glib::init();
  initPars();
}


/****************************************************************************************/
VideoHandler::~VideoHandler()
{
  if (!m_stopped)
    stop();
};


/****************************************************************************************/
void VideoHandler::initPars() {

  scaled = false;
  plotif = false;
  colorfeature = false;
  difffeature = false;
  computeSegments = true;

  resizeif = false;
  resizescale = 1.;

  inverted = false;
  
  nprobe = 7;
  Delta = 0;

  fixedSizeImage = cv::Size(0,0);
	  
  minWidth = 2 ;
  minExtent = 2;
  minArea = 2;
  maxArea = 10000;
  maxExtent = 10000;
  featureheight = 100;
  featurewidth =  20;
  m_featureSize  = Size(featureheight,featurewidth); // transposed image
  vector<float> tmp(3);
  Scale = tmp;
  for (int i=0; i<3;i++)
    Scale[i] = 0.33333;

}

/****************************************************************************************/
bool VideoHandler::testValid(Segment * pSeg) {

  double  extent = pSeg->MajorAxisLength + pSeg->MinorAxisLength;
  if ((extent<minExtent) || (extent>maxExtent) || (pSeg->Area<minArea)|| (pSeg->Area>maxArea) || (pSeg->MinorAxisLength<minWidth)) 
    return false;
  else 
    return true;
}

/****************************************************************************************/
void VideoHandler::plotFrame(cv::Mat frame,const string windowName) {

  if (frame.size().width>0) {
    Mat smallFrame;
    Size size(400,400);
    if (frame.size().width>0) {
      resize(frame,smallFrame,size);
      namedWindow( windowName, WINDOW_AUTOSIZE );
      imshow( windowName, smallFrame);
      waitKey(1);
    }
  }
}

/****************************************************************************************/
void VideoHandler::startThreads() {
  m_nextFrameThread = Glib::Threads::Thread::create(sigc::mem_fun(*this, &VideoHandler::readNextFrameThread));
  Glib::usleep(200);
  m_segmentThread = Glib::Threads::Thread::create(sigc::mem_fun(*this, &VideoHandler::segmentThread));
  Glib::usleep(200);
  m_threadsAlive = true;
}

/****************************************************************************************/
void VideoHandler::deleteThreads() {

  if (m_threadsAlive) {
    m_keepNextFrameThreadAlive = false;
    m_keepSegmentThreadAlive = false;

    m_emptyNextFrameCond.signal();
    m_emptySegmentCond.signal();
    //m_availableSegmentCond.signal();
    //m_availableNextFrameCond.signal();
    
    while ((!m_nextFrameThreadFinished) || (!m_segmentThreadFinished)) {
      Glib::usleep(200);
    }
    
    if ((m_segmentThread!=NULL)) {
      m_segmentThread->join();
    }

    if ((m_nextFrameThread!=NULL)) {
      m_nextFrameThread->join();
    }

    m_threadsAlive = false;
    
  }
}

/****************************************************************************************/
void VideoHandler::reinitThreads() {

  if (!m_stopped) {

    if (!m_camera) {
      waitThreads(); // Nextframe will be second last, segment last frame. 
      
      Glib::Threads::Mutex::Lock lock(m_NextFrameMutex);
      int iframe;
      iframe = pVideoCapture->get(cv::CAP_PROP_POS_FRAMES);
 
      if (iframe>0) {
	iframe = iframe-2; // two step  back;
	iframe = iframe<0?0:iframe;
	pVideoCapture->set(cv::CAP_PROP_POS_FRAMES,iframe);
      }

    }
    // signal that a new frame needs to be collected
    { Glib::Threads::Mutex::Lock lock(m_NextFrameMutex);
      m_availableNextFrame = false;
      m_emptyNextFrameCond.signal();
    }
    {// wait frame to finish
      Glib::Threads::Mutex::Lock lock(m_NextFrameMutex);
      while ((!m_availableNextFrame) && (!m_nextFrameThreadFinished) && (!m_segmentThreadFinished))
	m_availableNextFrameCond.wait(m_NextFrameMutex);
    }
    // signal that a new frame can be handled
    { Glib::Threads::Mutex::Lock lock(m_SegmentMutex);
      m_availableSegment = false;
      m_emptySegmentCond.signal();
    }
  }
}


/****************************************************************************************/
void VideoHandler::waitThreads() {

  { // first wait segment to finish
    Glib::Threads::Mutex::Lock lock(m_SegmentMutex);
    while ((!m_availableSegment) && (!m_nextFrameThreadFinished) && (!m_segmentThreadFinished))
      m_availableSegmentCond.wait(m_SegmentMutex);
  }

  {// wait next frame to finish
    Glib::Threads::Mutex::Lock lock(m_NextFrameMutex);
    while ((!m_availableNextFrame) && (!m_nextFrameThreadFinished) && (!m_segmentThreadFinished))
      m_availableNextFrameCond.wait(m_NextFrameMutex);
  }

  
}


/****************************************************************************************/
int VideoHandler::start() {

  if (m_stopped) {
    if ( m_camera) {
      FlyCapture2::BusManager busMgr;
      FlyCapture2::Error error;
      FlyCapture2::PGRGuid guid;
      
      cout << m_camIdx << endl;
      error = busMgr.GetCameraFromIndex( m_camIdx, &guid );
      if (error != PGRERROR_OK)
      {
	error.PrintErrorTrace();
	return -1;
      }
      cout << "got camera from index" << endl;
      pVideoSaver = new VideoSaver();
      cout << "created VideoSaver" << endl;
      if (pVideoSaver->init(guid)!=0){
	cout << "Warning: Error in initializing the camera\n" ;
	return -1;
      }
      if (fname.empty()) {
	cout << "start capture thread" << endl;
	if (pVideoSaver->startCapture()!=0) {
	  cout << "Warning: Error in starting the capture thread the camera \n";
	  return -1;
	}
      }
      else {
	cout << "start capture and write threads" << endl;
	if (pVideoSaver->startCaptureAndWrite(fname,string("X264"))!=0) {
	  cout << "Warning: Error in starting the writing thread the camera \n";
	  return -1;
	}
      }
    }
      else  {
      pVideoCapture = new VideoCapture(fname);
      //pVideoCapture->set(cv::CAP_PROP_CONVERT_RGB,true); // output RGB
    }

    if (knnMethod)
      pBackgroundSubtractor =  cv::createBackgroundSubtractorKNN(250,400,false);
    else
      pBackgroundThresholder =  new BackgroundThresholder();

    m_stopped=false;
    startThreads();
  }
  return 0;
}

/****************************************************************************************/
void VideoHandler::stop() {

  destroyAllWindows();

  if (!m_stopped) {

    deleteThreads();

    if (!m_camera)
      pVideoCapture.release();
    else {
      pVideoSaver.release();
    }

    if (knnMethod)
      pBackgroundSubtractor.release();
    else
      pBackgroundThresholder.release();
  } else {
    cout << "Warning: VideoHandler already stopped" << endl;
  }
  m_stopped = true;

    
}
/****************************************************************************************/
void VideoHandler::getOFrame(cv::Mat * pFrame) {
  {
    Glib::Threads::Mutex::Lock lock(m_SegmentMutex);
    if ((!m_camera) && (colorfeature))
      cvtColor(m_OFrame,*pFrame,cv::COLOR_BGR2RGB);
    else
      m_OFrame.copyTo(*pFrame); // might be not the latest one. Actually, likley the next one
  }
}
/****************************************************************************************/
void VideoHandler::getBWImg(cv::Mat * pBWImg) {
  {
    Glib::Threads::Mutex::Lock lock(m_SegmentMutex);
    m_BWImg.copyTo(*pBWImg); // might be not the latest one. Actually, likley the next one
  }
}

/****************************************************************************************/
int VideoHandler::step(vector<Segment> * pSegs,  double * pTimeStamp, cv::Mat * pFrame){

  
  if (m_stopped)  {
    if (start()!=0) {
      cout << "WARNING: step could not be excecuted" << endl;
      vector<Segment> tmp(0);
      *pSegs = tmp;
      *pFrame = Mat::zeros(0,0,CV_8UC1);
      return -1;
    }
  }

  cv::Mat bwimg;
  
  {
    Glib::Threads::Mutex::Lock lock(m_SegmentMutex);
    while ((!m_availableSegment) && (!m_nextFrameThreadFinished) && (!m_segmentThreadFinished))
      m_availableSegmentCond.wait(m_SegmentMutex);

    if ((m_nextFrameThreadFinished)|| (m_segmentThreadFinished)) {
      cout << "WARNING: threads are finished or not started. step not excecuted" << endl;
      vector<Segment> tmp(0);
      *pSegs = tmp;
      *pFrame = Mat::zeros(0,0,CV_8UC1);
      return -1;
    }
      
    // % always use the processed frame
    m_Frame.copyTo(*pFrame);


    vector <Segment> tmp(0);
    m_Segments = m_NextSegments;
    m_NextSegments = tmp; // clear memory

    *pSegs = m_Segments; // this should evoke a copy... ?!?!
    *pTimeStamp = m_TimeStamp;

    if (plotif)
      m_BWImg.copyTo(bwimg);

    // signal stuff
    m_availableSegment = false;
    m_emptySegmentCond.signal();
  }

  if (plotif) {
    bwimg.setTo(255,bwimg>0);
    plotFrame(bwimg,"bwimg"); 
    plotFrame(*pFrame,"frame");
  }
  return 0;
};

/****************************************************************************************/
void VideoHandler::readNextFrameThread()
{

  m_keepNextFrameThreadAlive = true;  
  m_availableNextFrame = false;
  m_nextFrameThreadFinished = false;
  
  while (m_keepNextFrameThreadAlive) {

    {
      Glib::Threads::Mutex::Lock lock(m_NextFrameMutex);
      // wait for segment thread  handling
      while ((m_availableNextFrame) && (m_keepNextFrameThreadAlive)) {
	m_emptyNextFrameCond.wait(m_NextFrameMutex);
      }
      if (!m_keepNextFrameThreadAlive)
	continue;

#ifdef DEBUG
      Glib::Timer timer;
      timer.start();
#endif

      Mat oframe;
      if (m_camera) {

	int frameNumber;
	if (pVideoSaver->getFrame(&oframe,&m_NextTimeStamp,&frameNumber)!=0) {
	  break;
	} // returns RGB
#ifdef DEBUG
	cout << oframe.size() << "  " << frameNumber << " " <<  m_NextTimeStamp << endl;
#endif
      
      } else {

	m_NextTimeStamp = pVideoCapture->get(cv::CAP_PROP_POS_MSEC)*1e-3; 
	if (!pVideoCapture->read(oframe)) { // should return RGB instead of BGR. But seems not to work.. NOW BGR
	  cout <<  "ERROR: File read error." << endl;
	  break;
	}
	//cvtColor(oframe,oframe,cv::COLOR_BGR2RGB);

#ifdef DEBUG
	cout << oframe.size() << " Time: " <<  m_NextTimeStamp << endl;
#endif

	if (oframe.type()!=CV_8UC3) {
	  cout <<  "ERROR: Expect CV_8UC3 color movie." << endl;
	  break;
	}
      }
    
      if ((resizeif) && (resizescale!=1)) {
	cv::resize(oframe,oframe,Size(0,0),resizescale,resizescale);
      }
    
      Mat frame;
    
      if (scaled) {
	cv::Mat channel[3];
	cv::Mat floatframe;
	vector<int> order(3);
	split(oframe, channel);
	for (int ii=0; ii<3; ii++) {
	  channel[ii].convertTo(channel[ii], CV_32FC1);
	  if (m_camera)
	    order[ii] = ii; // rgb
	  else
	    order[ii] = 2-ii; // bgr
	}
	floatframe =  ((float) Scale[0]/255.)*channel[order[0]] + ((float) Scale[1]/255.)*channel[order[1]] + ((float) Scale[2]/255.)*channel[order[2]] + (Delta); 
      
	// subtract mean
	//cv::Scalar globalmean = cv::mean(floatframe) ; // one channel anyway
	//floatframe -= (globalmean[0]-0.5);// 0.5 because should be between 0..1
	//floatframe.convertTo(frame, CV_8UC1,255.,(-globalmean[0]+0.5)*255); // better convert it back because of boundaries
	floatframe.convertTo(frame, CV_8UC1,255.);
      }
      else {
	if (m_camera)
	  cvtColor(oframe,frame,cv::COLOR_RGB2GRAY);
	else
	  cvtColor(oframe,frame,cv::COLOR_BGR2GRAY); // seems not to work for video capture...
      }
    
    
#ifdef DEBUG
      cout << "reading frame: " <<  timer.elapsed() << endl;  
      timer.start();
#endif
    
      // background computation (THIS IS THE SLOWEST PART!)
      if (knnMethod) {
	if (inverted) {
	  Mat iframe;
	  iframe = 255-frame;
	  pBackgroundSubtractor->apply(iframe, m_NextBWImg);
	}  else {
	  pBackgroundSubtractor->apply(frame, m_NextBWImg);
	}
      } else {
	pBackgroundThresholder->apply(frame, &m_NextBWImg, &m_NextDFrame);

      }
      
    
#ifdef DEBUG
      cout << "background sub: " <<  timer.elapsed() << endl;  
#endif

      m_NextFrame = frame;
      if (colorfeature)
	m_NextOFrame = oframe;
      else
	m_NextOFrame = frame;
    
      // thread stuff signal the end
      m_availableNextFrame = true;
      m_availableNextFrameCond.signal();
    }
  }
  m_nextFrameThreadFinished = true;
};
/****************************************************************************************/
void VideoHandler::segmentThread()
{

  m_keepSegmentThreadAlive = true;  
  m_availableSegment = false;
  m_segmentThreadFinished = false;
  
  while (m_keepSegmentThreadAlive) {
    
    {
      Glib::Threads::Mutex::Lock lock(m_SegmentMutex);
      
      while ((m_availableSegment) && (m_keepSegmentThreadAlive)){
	m_emptySegmentCond.wait(m_SegmentMutex);
      }
      if (!m_keepSegmentThreadAlive)
	continue;

      // start new Segment computation
      {
	Glib::Threads::Mutex::Lock lock(m_NextFrameMutex);
	
	while ((!m_availableNextFrame) && (m_keepSegmentThreadAlive)
	       && (!m_nextFrameThreadFinished)) {
	  m_availableNextFrameCond.wait(m_NextFrameMutex);
	}
	
	if ((!m_keepSegmentThreadAlive) || (m_nextFrameThreadFinished))
	  continue;

	if ((!colorfeature) && (difffeature) && (!knnMethod))
	  m_NextDFrame.copyTo(m_Frame);
	else
	  m_NextFrame.copyTo(m_Frame);
	
	if (colorfeature) 
	  m_NextOFrame.copyTo(m_OFrame);
	else if ((difffeature) && (!knnMethod)) // not strictly needed..
	  m_NextFrame.copyTo(m_OFrame);
	else
	  m_OFrame = m_Frame;

	  
	m_NextBWImg.copyTo(m_BWImg);
	
	m_TimeStamp = m_NextTimeStamp;
	
	m_availableNextFrame = false; 
	m_emptyNextFrameCond.signal();
      }
      
      
      if (computeSegments) {
	  
#ifdef DEBUG
	Glib::Timer timer;
	timer.start();
#endif
	
	// finding contours
	vector<vector<cv::Point> > contours;
	findFishContours(m_BWImg,&contours);
	
	
#ifdef DEBUG 
	cout << "contours: " <<  timer.elapsed() << endl;
	timer.start();
#endif
	
	Segment segm;
	vector <Segment> tmp(0);
	m_NextSegments.clear();
	m_NextSegments = tmp;

	int ii=0;
	int ssize;
	ssize = contours.size();
	if (ssize>MAXCONTOUR)
	  ssize = MAXCONTOUR;
	
	for( int i = 0; i< ssize; i++ ) {

	  getSegment(&segm,contours[i],m_BWImg,m_Frame,m_OFrame);
	      
	  if (testValid(&segm)) {
	    ii++;
	    m_NextSegments.push_back(segm);
	  }
	  if (ii>=MAXVALIDCONTOUR)
	    break;
	}
	
	
#ifdef DEBUG 
	cout << m_NextSegments.size() << endl;
	cout << "segments: " <<  timer.elapsed() << endl;  
#endif
      }
      
      // signal stuff
      m_availableSegment = true; 
      m_availableSegmentCond.signal();
    }
  }
  m_segmentThreadFinished = true;
}

/****************************************************************************************/
void VideoHandler::findFishContours(Mat inBwImg, vector<vector<cv::Point> > * newcontours) {

  // finding contours
  vector<Vec4i> hierarchy;
  vector<vector<cv::Point> > contours;

  findContours(inBwImg,contours,hierarchy,RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE, Point(0, 0));

  vector<bool>valid(contours.size());
  for (int i=0;i<contours.size();i++) {
    if ((contours[i].size()>10) && (contours[i].size()<maxArea)) // at least ten pixels
      valid[i] = true;
    else
      valid[i] = false;
  }

  Mat srel = getStructuringElement(MORPH_ELLIPSE,Size(3,3));

  // get the bounding boxes and erode 
  double mx =featureheight*featurewidth*3;
  for (int i=0;i<contours.size();i++) {
    if (!valid[i])
      continue;
    
    Rect bbox = boundingRect(contours[i]);
    double wh = bbox.width*bbox.height;
    if ((wh > mx) && (bbox.width+bbox.height<maxExtent) ){
      //try to split regions
      Mat roi;
      inBwImg(bbox).copyTo(roi); // make a true copy

      Mat tmp;
      roi.copyTo(tmp);
	
      // namedWindow("before",WINDOW_AUTOSIZE);
      // imshow("before",tmp>0);
      // waitKey(10);
      
      // MAYBE MORE EROSION ? 5x5 circular element,...
      dilate(roi,roi,srel,Point(-1,-1),1,BORDER_CONSTANT,0);
      erode(roi,roi,srel,Point(-1,-1),3,BORDER_CONSTANT,0);
      dilate(roi,roi,srel,Point(-1,-1),2,BORDER_CONSTANT,0);

      // namedWindow("after",WINDOW_AUTOSIZE);
      // imshow("after", roi>0);
      // waitKey(200);

      vector<vector<cv::Point> > localcontours;
      findContours(roi,localcontours,hierarchy,RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE, bbox.tl());
      int s=0;
      
      for (int j=0;j<localcontours.size();j++) {
	Rect bbox1 = boundingRect(localcontours[j]);
	if (bbox1.width*bbox1.height>minArea) {
	  (*newcontours).push_back(localcontours[j]);
	  s++;
	}
      }
    } else {
      if (bbox.width*bbox.height>minArea)  {
	(*newcontours).push_back(contours[i]);
      }
    }
  }
}
	


/****************************************************************************************/
void VideoHandler::getSegment(Segment * segm, vector<Point> inContour, Mat inBwImg, Mat inFrame, Mat inOFrame) {

  Mat srel = getStructuringElement(MORPH_ELLIPSE,Size(3,3));
  segm->Bbox = boundingRect(inContour);
  segm->Image = Mat(inBwImg(segm->Bbox)>0).clone();

  if (colorfeature)  {
    segm->FilledImage = inOFrame(segm->Bbox).clone();
  }
  else {
    segm->FilledImage = inFrame(segm->Bbox).clone();
  }

  // Mat Image2x;
  // getRectSubPix(inBwImg,Size(segm->Bbox.width*2,segm->Bbox.height*2),Point(segm->Bbox.x-segm->Bbox.width/2,segm->Bbox.y-segm->Bbox.height/2),Image2x,-1);
  // segm->Image2x = Image2x.clone()>0;

  segm->mback = mean(segm->FilledImage,segm->Image==0);
  //segm->FilledImage.setTo(segm->mback,msk);

  segm->Area  = (double) (cv::sum(segm->Image)[0])/255; // true seems to be 255...
  segm->MajorAxisLength = 0.;
  segm->MinorAxisLength = 0.;
  segm->Centroid = Point2f(segm->Bbox.x +segm->Bbox.width/2.,segm->Bbox.y + segm->Bbox.height/2.);

  
  // **** SEGMENT CALCULATION
  if ((inContour.size()>=5) && segm->Bbox.width>minWidth && segm->Bbox.height>minWidth){

    RotatedRect minEllipse = fitEllipse(inContour);
    segm->Orientation = minEllipse.angle;
 
    segm->Centroid = minEllipse.center; // better estimate
    segm->Size = minEllipse.size;

    if (segm->Size.width>segm->Size.height) {
      segm->MajorAxisLength = segm->Size.width;
      segm->MinorAxisLength = segm->Size.height;
    } else {
      segm->MajorAxisLength = segm->Size.height;
      segm->MinorAxisLength = segm->Size.width;
    }

    if (fixedSizeImage.width>0) {
      Mat tmpMat ;
      if (difffeature) {
	getRectSubPix(inFrame,fixedSizeImage,segm->Centroid,tmpMat,-1);
      } else {
	getRectSubPix(inOFrame,fixedSizeImage,segm->Centroid,tmpMat,-1);
      }
      segm->FilledImageFixedSize = tmpMat.clone(); // copy;
    }
    
    // get the rotation and the center points
    int nborder = 3;
    Mat img ;
    copyMakeBorder(segm->Image,img,nborder,nborder,nborder,nborder,BORDER_CONSTANT,0);

    //opening and closing
    dilate(img,img,srel,Point(-1,-1),2,BORDER_CONSTANT,0);
    erode(img,img,srel,Point(-1,-1),4,BORDER_CONSTANT,0);
    dilate(img,img,srel,Point(-1,-1),2,BORDER_CONSTANT,0);

    Mat F,L;
    distanceTransform(img,F,CV_DIST_L2, 3);
    Laplacian(F, L, -1, 5, 1, 0, BORDER_DEFAULT);//ksize 5

    // Threshold to get the center line points
    double minThres=0.66;
    double minVal,maxVal;
    minMaxIdx(L, &minVal, &maxVal, NULL, NULL, noArray());

    if ((minVal==0) && (minVal==0)) {
      // delete detection (see goodmsk)
      segm->MajorAxisLength = 0.;
      segm->MinorAxisLength = 0.;
      return;
    }
    
    threshold(L,L,minVal*minThres,255,THRESH_TOZERO_INV);

    vector<Point> locations;
    findNonZero(L!=0,locations);

    // //project points onto Major Axis
    Point2f localcenter(segm->Centroid.x-segm->Bbox.x,segm->Centroid.y-segm->Bbox.y);
    Point2f pvec(sin(segm->Orientation/180.*CV_PI),-cos(segm->Orientation/180.*CV_PI));
    vector<double> proj(locations.size());

    for (int i=0; i<locations.size();i++ ) {
      proj[i] = (locations[i].x-localcenter.x)*pvec.x + (locations[i].y-localcenter.y)*pvec.y;
    }

    Mat index;
    sortIdx(Mat(proj),index, CV_SORT_ASCENDING + CV_SORT_EVERY_COLUMN);

    // probe equidistant locations
    if (locations.size()<3) {
      segm->MajorAxisLength = 0.;
      segm->MinorAxisLength = 0.;
      return;
    }
      
    double pmax = proj[index.at<int>(locations.size()-2)];
    double pmin = proj[index.at<int>(1)];

    vector<Point> probe(nprobe);
    int j = 0;
    for (int i=0; i<nprobe;i++ ) {
      double ptarget = i*(pmax-pmin)/(nprobe-1) + pmin;
      while (proj[index.at<int>(j)]<ptarget && j<locations.size()-1) {
	j++;
      }
      probe[i] =  locations[index.at<int>(j)];
    } 

    // calulate thickness
    vector<double> thickness(nprobe);
    for (int i=0; i<nprobe; i++) {
      thickness[i] = F.at<float>(probe[i]);
    }
    
    if (thickness[0] + thickness[1]  < thickness[nprobe-1] + thickness[nprobe-2]) {
      // thinning. thus reverse points
      reverse(probe.begin(),probe.end());
      reverse(thickness.begin(),thickness.end());
      pvec= -pvec;
    }

    Mat centerLine(nprobe,2,CV_64FC1);
    for (int i=0; i<nprobe;i++ ) {
      probe[i].x = probe[i].x-nborder;
      probe[i].y = probe[i].y-nborder;
      centerLine.at<double>(i,0) = (double) probe[i].x + segm->Bbox.x;
      centerLine.at<double>(i,1) = (double) probe[i].y + segm->Bbox.y;
    }
    segm->CenterLine = centerLine;
    segm->Thickness = thickness;


    Point2f v(probe[0]-probe[2]);
    if (probe[0]==probe[2]) {  // could happen. Take the ellipse orientation instead
      v = pvec;
    } 
      
    v = v/norm(v);
    
    Point2f center(probe[0]); // center for rotations

#ifdef PLOTSEGMENTS
    if (segm->Size.height + segm->Size.width > 100 ){
      namedWindow("proj",WINDOW_AUTOSIZE);
      Mat timg = segm->FilledImage;
      Mat timg2;
      getRectSubPix(timg,Size(300,300),localcenter,timg2,-1);

      Point2f p1(150,150);
      Size sz= segm->Size;
      
      line(timg2,pvec*10+p1,p1,CV_RGB(255,255,145),1,8,0);
      line(timg2,v*10 + p1,  p1,CV_RGB(255,255,245),2,8,0);
      for (int i=0;i<nprobe;i++) {
	Point2f p2(probe[i]);
	circle(timg2,p2 + p1 - localcenter,5,CV_RGB(255-i*20,255-i*20,245-i*20),2,8,0);
      }
      circle(timg2,center + p1 - localcenter,5,CV_RGB(255,255,255),4,8,0);
      ellipse(timg2,p1,(sz)/2,segm->Orientation,0,360,CV_RGB(255,255,145),2,8,0);

      imshow("proj",timg2);
    }
#endif
    
    // get rotation matrix
    Size2f szin(segm->Image.size());
    Size2f szout(szin.height*fabs(v.y)+szin.width*fabs(v.x),szin.height*fabs(v.x)+szin.width*fabs(v.y));
    Point2f offs(szout.width-2*thickness[0],szout.height/2);
    offs = offs - center;
    
    Mat T = (Mat_<float>(2,3) << v.x, v.y, (1-v.x)*center.x-v.y*center.y + offs.x, -v.y,v.x,v.y*center.x+(1-v.x)*center.y + offs.y );

    // rotate the Images
    warpAffine(segm->Image,segm->RotImage,T,szout ,INTER_CUBIC,BORDER_CONSTANT,0);
    warpAffine(segm->FilledImage,segm->RotFilledImage,T,szout,INTER_CUBIC,BORDER_CONSTANT,segm->mback);

    Mat RotFilledMskImage;
    RotFilledMskImage = segm->RotFilledImage.clone();
    RotFilledMskImage.setTo(segm->mback,segm->RotImage==0);

    Point2f centerfish(szout.width - m_featureSize.width/2,szout.height/2);
    getRectSubPix(RotFilledMskImage,m_featureSize,centerfish,segm->FishFeature,-1);

    
    // get the output to matlab right
    transpose(segm->FishFeature,segm->FishFeature);
    flip(segm->FishFeature,segm->FishFeature,0);
    segm->FishFeature.convertTo(segm->FishFeature,CV_32FC1);

    // get bending stddev
    int cols = segm->RotFilledImage.cols;
    int rows = segm->RotFilledImage.rows;
    
    vector <Point2f> rotprobe(nprobe);
    for (int i=0;i<nprobe;i++) {
      rotprobe[i].x = T.at<float>(0,0)*probe[i].x + T.at<float>(0,1)*probe[i].y + T.at<float>(0,2);
      rotprobe[i].y = T.at<float>(1,0)*probe[i].x + T.at<float>(1,1)*probe[i].y + T.at<float>(1,2);
    }
    
    Mat comy = Mat::zeros(1,cols,CV_32FC1);
    for (int j=0;j<cols;j++) {
      float minm = cols;
      int mini = 0;
      for (int i=0;i<nprobe;i++) {
	float m = abs(rotprobe[i].x - (float) j);
	if (m<minm) {
	  mini = i;
	  minm = m;
	}
      }
      comy.at<float>(j) = rotprobe[mini].y;
    }
    
    Mat temp;
    Mat stdValue;
    meanStdDev(comy, temp, stdValue,comy>0);
    segm->bendingStdValue = stdValue.at<double>(0);

    
    // rotate the fixed size
    if (fixedSizeImage.width>0) {
      Mat tmpMat,rotTmpMat2x,rotTmpMat;

      Size fixedSize2x(fixedSizeImage.width*2.,fixedSizeImage.height*2.);
      Point2f center2(fixedSize2x.width/2.,fixedSize2x.height/2.); 
      Mat T2 = (Mat_<float>(2,3) << v.x, v.y, (1-v.x)*center2.x-v.y*center2.y, -v.y,v.x,v.y*center2.x+(1-v.x)*center2.y);

      if (difffeature) {
	getRectSubPix(inFrame,fixedSize2x,segm->Centroid,tmpMat,-1);
      }
      else {
	getRectSubPix(inOFrame,fixedSize2x,segm->Centroid,tmpMat,-1);
      }
      warpAffine(tmpMat,rotTmpMat2x,T2,fixedSize2x,INTER_CUBIC,BORDER_CONSTANT,segm->mback);
      getRectSubPix(rotTmpMat2x,fixedSizeImage,center2,rotTmpMat,-1);
      segm->FilledImageFixedSizeRotated = rotTmpMat.clone(); // copy;

    }

    
    // REMAP ??
    bool REMAP=false;
    
    if (REMAP) {
      
      int nconv = ceil(((float) cols)/nprobe/2);
      if (nconv==0) {
	return;
      }


      Mat kernel = getGaussianKernel(2*nconv+1, -1, CV_32FC1);
      filter2D(comy, comy, -1 , kernel.t(), Point(-1,-1), 0, BORDER_REPLICATE );

      Mat x = Mat::zeros(1,cols,CV_32FC1);
      for (int i=0;i<cols;i++) {
	x.at<float>(i) = (float) i;
      }
      Mat y = Mat::zeros(rows,1,CV_32FC1);
      for (int i=0;i<rows;i++) {
	y.at<float>(i) = (float) i;
      }

    
 
      // remap
      Mat Y;
      repeat(y,1,cols,Y);
      Mat X;
      repeat(x,rows,1,X);

      for (int i=0;i<cols;i++) {
	Y.col(i) += comy.at<float>(i) - (float) rows /2.;
      }
      //Y = max(min(Y,rows-1),0);

      // Mat map1,map2;
      // convertMaps(X,Y,map1,map2,CV_16SC2,false);
      Mat fishfeature;
      remap(segm->RotFilledImage,fishfeature,X,Y,INTER_CUBIC,BORDER_REPLICATE,segm->mback);

      //Point2f centerfish(szout.width - m_featureSize.width/2,szout.height/2);
      Mat ffr;
      getRectSubPix(fishfeature,m_featureSize,centerfish,ffr,-1);
      segm->FishFeatureRemap = ffr.clone(); //copy
      
      // get output right
      transpose(segm->FishFeatureRemap,segm->FishFeatureRemap);
      flip(segm->FishFeatureRemap,segm->FishFeatureRemap,0);
      segm->FishFeatureRemap.convertTo(segm->FishFeatureRemap,CV_32FC1);

    }
#ifdef PLOTSEGMENTS
    if (segm->Size.height + segm->Size.width > 100 ){
      namedWindow("fish",WINDOW_AUTOSIZE);
      imshow("fish",segm->FishFeature);


      
      Mat timg6;
      Point2f p8(fishfeature.size());
      getRectSubPix(fishfeature,Size(300,300),p8/2,timg6,-1);
      namedWindow("fishfeature",WINDOW_AUTOSIZE);
      imshow("fishfeature",timg6);


      Mat timg4(segm->RotFilledImage);
      for (int i=0;i<nprobe;i++) {
    	circle(timg4, rotprobe[i],4,CV_RGB(0,255,245),1,8,0);
      }

      for (int i=0;i<comy.cols;i++) {
    	circle(timg4, Point2f(i,comy.at<float>(i)),1,CV_RGB(0,255,245),1,8,0);
      }
 
      Mat timg3;
      Point2f p5 = Point2f(szout);
      getRectSubPix(timg4,Size(300,300),p5/2,timg3,-1);

      namedWindow("rotimg",WINDOW_AUTOSIZE);
      imshow("rotimg",timg3);
      waitKey(0); 
    }
#endif
  }
  
};

/****************************************************************************************/
int VideoHandler::set(const string prop, double value){

  if (prop=="scaled") {
    waitThreads();
    scaled = (bool) (value!=0);
    reinitThreads();
  }
  else  if (prop=="delta") {
    waitThreads();
    Delta = (float) value;
    reinitThreads();
  }
  else  if (prop=="resizeif") {
    waitThreads();
    resizeif = (bool) (value!=0);
    
    if (knnMethod)
      pBackgroundSubtractor->clear(); 
    else
      pBackgroundThresholder->clear(); 

    reinitThreads();
  }
  else if (prop=="resizescale") {
    waitThreads();
    if (knnMethod)
      pBackgroundSubtractor->clear(); 
    else
      pBackgroundThresholder->clear(); 

    resizescale = (float) (value);
    reinitThreads();
  }
  else if (prop=="featureheight") {
    waitThreads();
    featureheight = (int) value;
    m_featureSize.width = featureheight; // transposed
    reinitThreads();
  }
  else if (prop=="featurewidth") {
    waitThreads();
    featurewidth = (int) value;
    m_featureSize.height = featurewidth; // transposed
    reinitThreads();

  }
  else if (prop=="fixedSize") {
    waitThreads();
    fixedSizeImage = cv::Size((int) value,(int) value);
    reinitThreads();
  }
  else if (prop=="plotif") {
    plotif = (bool) (value!=0);
  }
  else if (prop=="colorfeature") {
    waitThreads();
    colorfeature = (bool) (value!=0);
    reinitThreads();
  }
  else if (prop=="difffeature") {
    waitThreads();
    difffeature = (bool) (value!=0);
    reinitThreads();
  }
  else if (prop=="nprobe") {
    waitThreads();
    nprobe = (int) value;
    reinitThreads();
  }
  else if (prop=="minWidth") {
    minWidth = (int) value;
  }
  else if (prop=="minExtent") {
    minExtent = (int) value;
  }
  else if (prop=="maxExtent") {
    maxExtent = (int) value;
  }
  else if (prop=="maxArea") {
    maxArea = (int) value;
  }
  else if (prop=="minArea") {
    minArea = (int) value;
  }
  else if (prop=="computeSegments") {
    computeSegments = (bool) value!=0;
  }
  else if ((prop=="inverted") && (!m_stopped)) {
    inverted = (bool) value!=0;
    if (!knnMethod) {
      waitThreads();
      pBackgroundThresholder->setInverted(inverted);
      reinitThreads();
    }
  }
  else if (prop=="timePos") {
    if ((!m_camera) && (!m_stopped)) {
      waitThreads();
      pVideoCapture->set(cv::CAP_PROP_POS_MSEC,value);
      reinitThreads();
    }  else {
      // do nothing
    }
  }
  else if (prop=="FPS")  {
    if ((!m_camera) && (!m_stopped))
      pVideoCapture->set(cv::CAP_PROP_FPS,value);
    // else do nothing.. might want to implement a change in FPS. but then writing has to stop... 
  }
  else if (prop == "PosFrames") {
    if ((!m_camera) && (!m_stopped)) 
      pVideoCapture->set(cv::CAP_PROP_POS_FRAMES,(int) value);
    // else do nothing.. 
  }
  else if ((prop=="DetectShadows") && (!m_stopped) && (knnMethod) ) {
    waitThreads();
    pBackgroundSubtractor->setDetectShadows(value!=0);
    reinitThreads();
  }
  else if ((prop == "Dist2Threshold")&& (!m_stopped) && (knnMethod)) {
    waitThreads();
    pBackgroundSubtractor->setDist2Threshold(value);
    reinitThreads();
  }
  else if ((prop == "History")&& (!m_stopped)) {
    waitThreads();
    if (knnMethod)
      pBackgroundSubtractor->setHistory((int) value);
    else
      pBackgroundThresholder->setHistory((int) value);
    reinitThreads();
  }
  else if ((prop == "kNNSamples")&& (!m_stopped) && (knnMethod)) {
    waitThreads();
    pBackgroundSubtractor->setkNNSamples((int)value);
    reinitThreads();
  }
  else if ((prop == "NSamples")&& (!m_stopped) && (knnMethod)) {
    waitThreads();
    pBackgroundSubtractor->setNSamples((int) value);
    reinitThreads();
  }
  else if ((prop == "ShadowThreshold")&& (!m_stopped) && (knnMethod)) {
    waitThreads();
    pBackgroundSubtractor->setShadowThreshold(value);
    reinitThreads();
  }
  else if ((prop == "ShadowValue")&& (!m_stopped)&& (knnMethod)) {
    waitThreads();
    pBackgroundSubtractor->setShadowValue(value);
    reinitThreads();
  }
  else if ((prop == "nskip")&& (!m_stopped)) {
    if (!knnMethod) {
      waitThreads();
      pBackgroundThresholder->setNSkip((int) value);
      reinitThreads();
    }
  }
  else if ((prop == "adjustThresScale")&& (!m_stopped)) {
    if (!knnMethod) {
      waitThreads();
      pBackgroundThresholder->setThresScale((double) value);
      reinitThreads();
    }
  }
  else {
    cout <<  "ERROR: Could not set property " << prop <<"! \n";
    return -1;
  }
  return 0;
}

/****************************************************************************************/
double VideoHandler::get(const string prop){

  if (prop=="scaled") {
    return (double) scaled;
  }
  else if (prop=="plotif") {
    return(double) plotif;
  }
  else if (prop=="fixedSize") {
    return(double) fixedSizeImage.width;
  }
  else if (prop=="camera") {
    return(double) m_camera;
  }
  else if (prop=="delta") {
    return(double) Delta;
  }
  else if (prop=="resizeif") {
    return(double) resizeif;
  }
  else if (prop=="knnMethod") {
    return(double) knnMethod;
  }
  else if (prop=="inverted") {
    return(double) inverted;
  }
  else if (prop=="resizescale") {
    return(double) resizescale;
  }
  else if (prop=="colorfeature") {
    return (double) colorfeature;
  }
  else if (prop=="difffeature") {
    return (double) difffeature;
  }
  else if (prop=="featureheight") {
    return (double) featureheight;
  }
  else if (prop=="featurewidth") {
    return (double) featurewidth;
  }
  else if (prop=="nprobe") {
    return (double) nprobe;
  }
  else if (prop=="computeSegments") {
    return (double) computeSegments;
  }
  else if (prop=="minWidth") {
    return (double) minWidth;
  }
  else if (prop=="minExtent") {
    return (double)  minExtent;
  }
  else if (prop=="maxExtent") {
    return (double) maxExtent;
  }
  else if (prop=="maxArea") {
    return (double) maxArea;
  }
  else if (prop=="timePos") {
    if ((!m_camera) && (!m_stopped))
      return ((double) pVideoCapture->get(cv::CAP_PROP_POS_MSEC))/1000.;
    else
      return (double) -1;
  }
  else if (prop=="minArea") {
    return (double) minArea;
  }   else  if ((prop == "DetectShadows") && (!m_stopped) ) {
    if (knnMethod)
      return (double) pBackgroundSubtractor->getDetectShadows();
    else
      return (double) -1;
  }
  else if ((prop == "Dist2Threshold")&& (!m_stopped)  ) {
    if (knnMethod)
      return (double) pBackgroundSubtractor->getDist2Threshold();
    else
      return (double) -1;
  }
  else if ((prop == "History")&& (!m_stopped) ) {
    if (knnMethod)
      return (double) pBackgroundSubtractor->getHistory();
    else
      return (double) pBackgroundThresholder->getHistory();
  }
  else if ((prop == "kNNSamples")&& (!m_stopped) ) {
    if (knnMethod) 
      return (double) pBackgroundSubtractor->getkNNSamples();
    else
      return -1;
  }
  else if ((prop == "NSamples")&& (!m_stopped) ) {
    if (knnMethod) 
      return (double) pBackgroundSubtractor->getNSamples();
    else
      return (double) -1;
  }
  else if ((prop == "ShadowThreshold")&& (!m_stopped) ){
    if (knnMethod)
      return (double) pBackgroundSubtractor->getShadowThreshold();
    else
      return (double) -1;
  }
  else if ((prop == "ShadowValue")&& (!m_stopped)) {
    if  (knnMethod)
      return (double) pBackgroundSubtractor->getShadowValue();
    else
      return (double) -1;
  }
  else if ((prop == "nskip")&& (!m_stopped)) {
    if  (knnMethod)
      return (double) -1;
    else
      return (double) pBackgroundThresholder->getNSkip();
  }
  else if ((prop == "adjustThresScale")&& (!m_stopped)) {
    if  (knnMethod)
      return (double) -1;
    else
      return (double) pBackgroundThresholder->getThresScale();
  }
  else if ((prop == "FrameWidth") && (!m_stopped)) {
    double width;
    if (!m_camera) 
      width =  (double) pVideoCapture->get(cv::CAP_PROP_FRAME_WIDTH);
    else  {
      Size frameSize = pVideoSaver->getFrameSize();
      width = (double) frameSize.width;
    }
    if (resizeif)
      width = round(width*resizescale);
    return width;
  }
  else if ((prop == "FrameHeight") && (!m_stopped)) {
    double height;
    if (!m_camera) 
      height = (double) pVideoCapture->get(cv::CAP_PROP_FRAME_HEIGHT);
    else  {
      Size frameSize = pVideoSaver->getFrameSize();
      height =(double) frameSize.height;
    }
    if (resizeif)
      height = round(height*resizescale);
    return height;
  }
  else if ((prop == "FPS")  && (!m_stopped)){
    if (!m_camera) 
      return (double) pVideoCapture->get(cv::CAP_PROP_FPS);
    else 
      return (double) pVideoSaver->getFPS();
  }
  else if ((prop == "FrameCount")  && (!m_stopped)){
    if (!m_camera) 
      return (double) pVideoCapture->get(cv::CAP_PROP_FRAME_COUNT);
    else
      return (double) 0;
  } else if ((prop == "PosFrames") && (!m_stopped)) {
    if (!m_camera) 
      return (double) pVideoCapture->get(cv::CAP_PROP_POS_FRAMES);
    else
      return (double) pVideoSaver->getCurrentFrameNumber();
  }
  else {
    cout <<  "ERROR: Could not read property " << prop <<". Maybe VideoHandler is not started! \n";
    return -1;
  }
  return 0;
}

// /** Capture Property map for option processing
//  */
// const ConstMap<std::string,int> CapProp = ConstMap<std::string,int>
//     ("PosMsec",       cv::CAP_PROP_POS_MSEC)       //!< Current position of the video file in milliseconds or video capture timestamp.
//     ("PosFrames",     cv::CAP_PROP_POS_FRAMES)     //!< 0-based index of the frame to be decoded/captured next.
//     ("AVIRatio",      cv::CAP_PROP_POS_AVI_RATIO)  //!< Relative position of the video file: 0 - start of the film, 1 - end of the film.
//     ("FrameWidth",    cv::CAP_PROP_FRAME_WIDTH)    //!< Width of the frames in the video stream.
//     ("FrameHeight",   cv::CAP_PROP_FRAME_HEIGHT)   //!< Height of the frames in the video stream.
//     ("FPS",           cv::CAP_PROP_FPS)            //!< Frame rate.
//     ("FourCC",        cv::CAP_PROP_FOURCC)         //!< 4-character code of codec.
//     ("FrameCount",    cv::CAP_PROP_FRAME_COUNT)    //!< Number of frames in the video file.
//     ("Format",        cv::CAP_PROP_FORMAT)         //!< Format of the Mat objects returned by retrieve() .
//     ("Mode",          cv::CAP_PROP_MODE)           //!< Backend-specific value indicating the current capture mode.
//     ("Brightness",    cv::CAP_PROP_BRIGHTNESS)     //!< Brightness of the image (only for m_cameras).
//     ("Contrast",      cv::CAP_PROP_CONTRAST)       //!< Contrast of the image (only for cameras).
//     ("Saturation",    cv::CAP_PROP_SATURATION)     //!< Saturation of the image (only for cameras).
//     ("Hue",           cv::CAP_PROP_HUE)            //!< Hue of the image (only for cameras).
//     ("Gain",          cv::CAP_PROP_GAIN)           //!< Gain of the image (only for cameras).
//     ("Exposure",      cv::CAP_PROP_EXPOSURE)       //!< Exposure (only for cameras).
//     ("ConvertRGB",    cv::CAP_PROP_CONVERT_RGB)    //!< Boolean flags indicating whether images should be converted to RGB.
//     //("WhiteBalance",cv::CAP_PROP_WHITE_BALANCE)  //!< Currently not supported
//     ("Rectification", cv::CAP_PROP_RECTIFICATION)  //!< Rectification flag for stereo cameras (note: only supported by DC1394 v 2.x backend currently)
// ;


/****************************************************************************************/
void VideoHandler::setScale(vector<float> scale) {
  waitThreads();
  for (int i=0;i<scale.size();i++) {
    Scale[i] = scale[i];
  }
  reinitThreads();
};

/****************************************************************************************/
vector< float> VideoHandler::getScale() {
  return Scale;
};

/****************************************************************************************/
void VideoHandler::resetBkg() {
  if (!m_stopped) {
    waitThreads();
    if (knnMethod) {
      pBackgroundSubtractor->clear(); // similar background anyway. do not reset 
    }
    else {
      pBackgroundThresholder->clear();
    }
    reinitThreads();
  }
};



//*******************************************************************************//
//* Thresholder 
//*******************************************************************************//

BackgroundThresholder::BackgroundThresholder() {
  m_history = 250;
  m_nskip = 5;
  m_adjustThresScale = 1;
  setInverted(false);
  m_istep = 0;
};

BackgroundThresholder::~BackgroundThresholder() {};

void BackgroundThresholder::clear() { // background image will be overwritten
  m_istep = 0;
}

void BackgroundThresholder::setHistory(int value) {
   m_history = value;
}
int BackgroundThresholder::getHistory() {
   return m_history;
}
void  BackgroundThresholder::setThresScale(double value) {
  m_adjustThresScale = value;
  if (m_threstype==THRESH_BINARY) // inverted 
    m_adjustThresScaleCorrected = 1 + 1 - m_adjustThresScale;
   else
    m_adjustThresScaleCorrected = m_adjustThresScale;
  return ;
}
double  BackgroundThresholder::getThresScale() {
  return m_adjustThresScale;
}

void BackgroundThresholder::setNSkip(int value) {
   m_nskip = value;
}
int BackgroundThresholder::getNSkip() {
   return m_nskip;
}


void BackgroundThresholder::setInverted(bool invertedif) {
  if (invertedif) {
    m_threstype = THRESH_BINARY;
  }
  else {
    m_threstype = THRESH_BINARY_INV;
  }
  setThresScale(m_adjustThresScale);
}

void BackgroundThresholder::apply(cv::Mat frame,cv::Mat *bwimg, cv::Mat *dframe) {

  if (frame.size().width==0)
    return;

  cv::Mat floatframe;
  frame.convertTo(floatframe,CV_32FC1);
  
  // subtract background

  if (m_istep>MINBACKTHRESSTEP)  { // has to be a rough estimate of the mean already
    *dframe = floatframe - m_meanImage;
    (*dframe).convertTo(*dframe,CV_8UC1);

    // apply threshold
    if (m_istep<MAXOTSU+MINBACKTHRESSTEP)  {
      m_thres = threshold(*dframe,*bwimg,m_thres,255,THRESH_OTSU + m_threstype);
    } else {
      threshold(*dframe,*bwimg,m_thres*m_adjustThresScaleCorrected,255, m_threstype);
    }
  } else {
    *bwimg = Mat::zeros(frame.size(),CV_8UC1);
    *dframe = Mat::zeros(frame.size(),CV_8UC1);
  }
  
  //update mean
  if (m_istep==0) {
    m_meanImage = floatframe - 127;
  }  else {
    int n;
    if (m_istep<m_history) {
      n = m_istep+1;
    } else {
      if ((m_istep % m_nskip)==0) {
	n = m_history/m_nskip;
      }  else {
	n = 0;
      }
    }

    if (n!=0) {
      // update 
      m_meanImage += 127;
      m_meanImage *= n-1;
      m_meanImage = (m_meanImage + floatframe)/n - 127;

    }
  }

  m_istep++;
}


/****************************************************************************************/
/****************************************************************************************/
// // // main
int main() {
  string vid;
  //vid = "/data/videos/longterm/longterm8.avi";
  //vid = "/home/malte/data/zebra/videos/longterm/longterm8.avi";
  vid = "/home/malte/Videos/Blongterm10.avi";

  VideoHandler vh(vid,true);
  //VideoHandler vh(0,"/home/malte/test_VeideoHandler.avi",false);
  cout << "init" << endl;
  vh.set("plotif",true);
  cout << "plotif" << endl;
  //namedWindow( "Patch", WINDOW_AUTOSIZE );
  vh.set("scaled",false);
  cout << "scaled" << endl;
  vh.set("colorfeature",true);
  cout << "colorfeature" << endl;

  vector<float> scale(3);
  scale[0] = -0.592655;
  scale[1] = -0.206795;
  scale[2] = 2.14107;
  vh.setScale(scale);
  
  cout << "set pars" << endl;
  vh.start();
  cout << "start loop " << endl;
  vh.set("inverted",false);
  vh.set("History",250);
  vh.set("History",250);
  vh.set("resizeif",false);
  vh.set("resizescale",0.5);;
  
  namedWindow("bwimg",WINDOW_AUTOSIZE);
  vector<Segment> segm;
  cv::Mat frame;
  double timeStamp;

  for (int i=0; i<100; i++) {
    if (vh.step(&segm,&timeStamp,&frame)!=0)
      break;
    cout << i << " (" << timeStamp <<") : Segment Size " <<  segm.size() << "\n";
    if (segm.size()>0)
      cout << "MajorAxisLength[0]: " << segm[0].MajorAxisLength << endl;
    waitKey(10);
  }
  vh.stop();
  vh.start();
  vh.stop();

  VideoHandler vh2(vid,true);

  for (int i=0; i<25; i++) {
    if (vh2.step(&segm,&timeStamp,&frame)!=0)
      break;
    waitKey(10);
    cout << i << ": Segment Size " << segm.size() << "\n";
  }
  vh2.stop();
  vh2.start();
  vh2.stop();
    
  return 0;
}


