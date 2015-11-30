

#include "VideoHandler.h"
//#define DEBUG
//#define PLOTSEGMENTS

#define MAXCONTOUR 250
#define MAXVALIDCONTOUR 50
#define MAXOTSU 50
#define MINBACKTHRESSTEP 20

using namespace std;
using namespace cv;


Segment::Segment() {}
Segment::~Segment() {};


VideoHandler::VideoHandler(const string fname,bool inKnnMethod)
{
  camera = false;
  m_fname = fname;
  stopped = true;
  m_threadExists = false;
  knnMethod = inKnnMethod;
    
  Glib::init();
  initPars();
};

/****************************************************************************************/
VideoHandler::VideoHandler(int camIdx, const string fname,bool inKnnMethod)
{
  camera = true;
  stopped = true;
  m_camIdx = camIdx;
  m_fname = fname;
  m_threadExists = false;
  knnMethod = inKnnMethod;

  Glib::init();
  initPars();
}


/****************************************************************************************/
VideoHandler::~VideoHandler()
{
  if (!stopped)
    stop();
};


/****************************************************************************************/
void VideoHandler::initPars() {

  scaled = false;
  plotif = false;
  colorfeature = false;
  computeSegments = true;

  resizeif = false;
  resizescale = 1.;

  inverted = false;
  
  nprobe = 7;
  m_Delta = 0;

  minWidth = 2 ;
  minExtent = 2;
  minArea = 2;
  maxArea = 10000;
  maxExtent = 10000;
  featureheight = 100;
  featurewidth =  20;
  featureSize  = Size(featureheight,featurewidth); // transposed image
  vector<float> tmp(3);
  m_Scale = tmp;
  for (int i=0; i<3;i++)
    m_Scale[i] = 0.33333;

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
void VideoHandler::plotFrame(cv::Mat frame) {

  Mat smallFrame;
  Size size(400,400);
  if (frame.size().width>0) {
    resize(frame,smallFrame,size);
    namedWindow( "Frame", WINDOW_AUTOSIZE );
    imshow( "Frame", smallFrame);
    waitKey(1);
  }

}

/****************************************************************************************/
void VideoHandler::startThread() {

  ///deleteThread();
  m_readThread = Glib::Threads::Thread::create(sigc::mem_fun(*this, &VideoHandler::_readNextFrameThread));
}
/****************************************************************************************/
void VideoHandler::waitThread() {
  Glib::Threads::Mutex::Lock lock(m_FrameMutex);
}

/****************************************************************************************/
void VideoHandler::deleteThread() {

  waitThread();
  if (m_threadExists) {
    //m_readThread->join();
      m_threadExists = false;
  }
}


/****************************************************************************************/
void VideoHandler::initialize() {

  if (!stopped) {

    waitThread();
    
    if (!camera) {
      int iframe;
      iframe = pVideoCapture->get(cv::CAP_PROP_POS_FRAMES);
 
      //cout << "frame: " << iframe << endl;
      if (iframe>0) {
	iframe = iframe-1; // step one back;
	pVideoCapture->set(cv::CAP_PROP_POS_FRAMES,iframe);
      }
    }

    // start a new thread
    startThread();
    Glib::usleep(200); // need to wait a bit to acquire the lock in the thread
  }
}

/****************************************************************************************/
int VideoHandler::start() {

  if (stopped) {
    if ( camera) {
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
      if (m_fname.empty()) {
	cout << "start capture thread" << endl;
	if (pVideoSaver->startCapture()!=0) {
	  cout << "Warning: Error in starting the capture thread the camera \n";
	  return -1;
	}
      }
      else {
	cout << "start capture and write threads" << endl;
	if (pVideoSaver->startCaptureAndWrite(m_fname,string("X264"))!=0) {
	  cout << "Warning: Error in starting the writing thread the camera \n";
	  return -1;
	}
      }
    }
      else  {
      pVideoCapture = new VideoCapture(m_fname);
      pVideoCapture->set(cv::CAP_PROP_CONVERT_RGB,true); // output RGB
    }

    if (knnMethod)
      pBackgroundSubtractor =  cv::createBackgroundSubtractorKNN(250,400,false);
    else
      pBackgroundThresholder =  new BackgroundThresholder();
    stopped=false;
    startThread();
    Glib::usleep(200); // need to wait a bit to acquire the lock in the thread
    waitThread();
  }
  return 0;
}

/****************************************************************************************/
void VideoHandler::stop() {

  deleteThread();
  destroyAllWindows();

  if (!stopped) {

    if (!camera)
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
  stopped = true;
}

/****************************************************************************************/
void VideoHandler::getOFrame(cv::Mat * pFrame) {
  *pFrame = m_OFrame;
}


/****************************************************************************************/
void VideoHandler::step(vector<Segment> * pSeg, cv::Mat * pOFrame, cv::Mat * pFrame,cv::Mat * pBWImg){

  if (stopped) {
    if (start()!=0) {
      cout << "WARNING: step not excecuted" << endl;
      return;
    }
  }
  
  waitThread();
  //output
  *pOFrame = m_NextOFrame.clone();
  *pFrame = m_NextFrame.clone();
  *pBWImg = m_NextBWImg.clone();

  m_OFrame = m_NextOFrame.clone();

  
  // start a new thread to load in the background
  startThread();
  
#ifdef DEBUG
  Glib::Timer timer;
  timer.start();
#endif


  // finding contours
  vector<Vec4i> hierarchy;
  vector<vector<cv::Point> > contours;
  if (computeSegments) {
    findContours(*pBWImg,contours,hierarchy,RETR_EXTERNAL, CHAIN_APPROX_SIMPLE, Point(0, 0));
  }


#ifdef DEBUG 
  cout << "contours: " <<  timer.elapsed() << endl;
  timer.start();
#endif


  if (computeSegments) {
    pSeg->clear();
    Segment segm;
    int ii=0;
    int ssize;
    ssize = contours.size();
    if (ssize>MAXCONTOUR)
      ssize = MAXCONTOUR;
    
    for( int i = 0; i< ssize; i++ ) {
      getSegment(&segm,contours[i],*pBWImg,*pFrame,*pOFrame);
      if (testValid(&segm)) {
	ii++;
	pSeg->push_back(segm);
      }
      if (ii>=MAXVALIDCONTOUR)
	break;
    }
  }

 
#ifdef DEBUG 
  cout << "segments: " <<  timer.elapsed() << endl;  
#endif
  
  if (plotif)
    plotFrame(*pFrame);

};

/****************************************************************************************/
void VideoHandler::_readNextFrameThread()
{
  Glib::Threads::Mutex::Lock lock(m_FrameMutex);
  if (stopped) {
    cout <<  "Error: VideoHandler stopped. Cannot read next Frame." << endl;
    return;
  }

#ifdef DEBUG
  Glib::Timer timer;
  timer.start();
  m_threadExists = true;
#endif
  
  Mat oframe,frame;
  if (camera) {

    double timeStamp;
    int frameNumber;
    if (pVideoSaver->getFrame(&oframe,&timeStamp,&frameNumber)!=0) {
      return;
    } // returns RGB
#ifdef DEBUG
    cout << oframe.size() << endl;
    cout << frameNumber << " " <<  timeStamp << endl;
#endif

  } else {
    if (!pVideoCapture->read(oframe)) { // should return RGB instead of BGR. But seems not to work..
      cout <<  "File read error";
      return;
    }
    cvtColor(oframe,oframe,cv::COLOR_BGR2RGB);
    
    if (oframe.type()!=CV_8UC3) {
      cout <<  "ERROR: Expect CV_8UC3 color movie";
      return;
    }
  }

  if ((resizeif) && (resizescale!=1)) {
    cv::resize(oframe,oframe,Size(0,0),resizescale,resizescale);
  }
  
  
  if (scaled) {
    cv::Mat channel[3];
    cv::Mat floatframe;
    split(oframe, channel);
    for (int ii=0; ii<3; ii++) {
      channel[ii].convertTo(channel[ii], CV_32FC1);
    }
    floatframe =  ((float) m_Scale[0]/255.)*channel[0] + ((float) m_Scale[1]/255.)*channel[1] + ((float) m_Scale[2]/255.)*channel[2];// + (m_Delta); //might not need delta because mean subtracted below anyway
    
    // subtract mean
    cv::Scalar globalmean = cv::mean(floatframe) ; // one channel anyway
    //floatframe -= (globalmean[0]-0.5);// 0.5 because should be between 0..1
    floatframe.convertTo(frame, CV_8UC1,255.,(-globalmean[0]+0.5)*255); // better convert it back because of boundaries
  }
  else
    cvtColor(oframe,frame,cv::COLOR_RGB2GRAY);


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
    pBackgroundThresholder->apply(frame, &m_NextBWImg);
  }
      

#ifdef DEBUG
  cout << "background sub: " <<  timer.elapsed() << endl;  
#endif
  // // finding contours
  // vector<Vec4i> hierarchy;
  // m_NextContours.clear();
  // findContours(m_NextBWImg,m_NextContours,hierarchy,cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE, Point(0, 0));
  // // get the fish patches and segments
  // vector<Segment> segms(m_NextContours.size());
  // for( int i = 0; i< m_NextContours.size(); i++ )
  //   getSegment(&(segms[i]),m_NextContours[i],m_NextBWImg,frame,oframe);
  
  // m_NextSegments = segms;
  m_NextFrame = frame;
  if (colorfeature)
    m_NextOFrame = oframe;
  else
    m_NextOFrame = frame;

};

/****************************************************************************************/
void VideoHandler::getSegment(Segment * segm, vector<Point> inContour, Mat inBwImg, Mat inFrame, Mat inOFrame) {


  Mat mContours = Mat(inContour);
  segm->Bbox = boundingRect(mContours);
  segm->Image = Mat(inBwImg(segm->Bbox)>0);
  if (colorfeature) 
    segm->FilledImage = Mat(inOFrame(segm->Bbox));
  else
    segm->FilledImage = Mat(inFrame(segm->Bbox));

  // getRectSubPix(inBwImg,Size(segm->Bbox.width*2,segm->Bbox.height*2),Point(segm->Bbox.x-segm->Bbox.width/2,segm->Bbox.y-segm->Bbox.height/2),segm->Image2x,-1);
  // segm->Image2x = segm->Image2x>0;
  // getRectSubPix(tmpFrame,Size(segm->Bbox.width*2,segm->Bbox.height*2),Point(segm->Bbox.x-segm->Bbox.width/2,segm->Bbox.y-segm->Bbox.height/2),segm->FilledImage2x,-1);

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
    
    // get the rotation and the center points
    int nborder = 3;
    Mat img ;
    copyMakeBorder(segm->Image,img,nborder,nborder,nborder,nborder,BORDER_CONSTANT,0);

    //opening and closing
    dilate(img,img,Mat(),Point(-1,-1),1,BORDER_CONSTANT,0);
    erode(img,img,Mat(),Point(-1,-1),2,BORDER_CONSTANT,0);
    dilate(img,img,Mat(),Point(-1,-1),1,BORDER_CONSTANT,0);

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

    Point2f centerfish(szout.width - featureSize.width/2,szout.height/2);
    getRectSubPix(RotFilledMskImage,featureSize,centerfish,segm->FishFeature,-1);

    // get the output to matlab right
    transpose(segm->FishFeature,segm->FishFeature);
    flip(segm->FishFeature,segm->FishFeature,0);
    segm->FishFeature.convertTo(segm->FishFeature,CV_32FC1);

    

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

    // REMAP ??
    bool REMAP=true;
    
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

      //Point2f centerfish(szout.width - featureSize.width/2,szout.height/2);
      getRectSubPix(fishfeature,featureSize,centerfish,segm->FishFeatureRemap,-1);

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
    waitThread();
    scaled = (bool) (value!=0);
    initialize();
  }
  else  if (prop=="delta") {
    waitThread();
    m_Delta = (float) value;
    initialize();
  }
  else  if (prop=="resizeif") {
    waitThread();
    resizeif = (bool) (value!=0);
    initialize();
  }
  else if (prop=="resizescale") {
    waitThread();
    resizescale = (float) (value);
    initialize();
  }
  else if (prop=="featureheight") {
    waitThread();
    featureheight = (int) value;
    featureSize.width = featureheight; // transposed
    initialize();
  }
  else if (prop=="featurewidth") {
    waitThread();
    featureheight = (int) value;
    featureSize.height = featurewidth; // transposed
    initialize();

  }
  else if (prop=="plotif") {
    plotif = (bool) (value!=0);
  }
  else if (prop=="colorfeature") {
    waitThread();
    colorfeature = (bool) (value!=0);
    initialize();
  }
  else if (prop=="nprobe") {
    waitThread();
    nprobe = (int) value;
    initialize();
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
  else if ((prop=="inverted") && (!stopped)) {
    inverted = (bool) value!=0;
    if (!knnMethod) {
      waitThread();
      pBackgroundThresholder->setInverted(inverted);
      initialize();
    }
  }
  else if (prop=="timePos") {
    if ((!camera) && (!stopped)) {
      waitThread();
      pVideoCapture->set(cv::CAP_PROP_POS_MSEC,value);
      //pBackgroundSubtractor->clear(); // similar background anyway. do not reset
      initialize();
    }  else {
      // do nothing
    }
  }
  else if (prop=="FPS")  {
    if ((!camera) && (!stopped))
      pVideoCapture->set(cv::CAP_PROP_FPS,value);
    // else do nothing.. might want to implement a change in FPS. but then writing has to stop... 
  }
  else if (prop == "PosFrames") {
    if ((!camera) && (!stopped)) 
      pVideoCapture->set(cv::CAP_PROP_POS_FRAMES,(int) value);
    // else do nothing.. 
  }
  else if ((prop=="DetectShadows") && (!stopped) && (knnMethod) ) {
    waitThread();
    pBackgroundSubtractor->setDetectShadows(value!=0);
    initialize();
  }
  else if ((prop == "Dist2Threshold")&& (!stopped) && (knnMethod)) {
    waitThread();
    pBackgroundSubtractor->setDist2Threshold(value);
    initialize();
  }
  else if ((prop == "History")&& (!stopped)) {
    waitThread();
    if (knnMethod)
      pBackgroundSubtractor->setHistory((int) value);
    else
      pBackgroundThresholder->setHistory((int) value);
    initialize();
  }
  else if ((prop == "kNNSamples")&& (!stopped) && (knnMethod)) {
    waitThread();
    pBackgroundSubtractor->setkNNSamples((int)value);
    initialize();
  }
  else if ((prop == "NSamples")&& (!stopped) && (knnMethod)) {
    waitThread();
    pBackgroundSubtractor->setNSamples((int) value);
    initialize();
  }
  else if ((prop == "ShadowThreshold")&& (!stopped) && (knnMethod)) {
    waitThread();
    pBackgroundSubtractor->setShadowThreshold(value);
    initialize();
  }
  else if ((prop == "ShadowValue")&& (!stopped)&& (knnMethod)) {
    waitThread();
    pBackgroundSubtractor->setShadowValue(value);
    initialize();
  }
  else if ((prop == "nskip")&& (!stopped)) {
    if (!knnMethod) {
      waitThread();
      pBackgroundThresholder->setNSkip((int) value);
      initialize();
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
  else if (prop=="camera") {
    return(double) camera;
  }
  else if (prop=="delta") {
    return(double) m_Delta;
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
  else if (prop=="featureheight") {
    return (double) featureheight;
  }
  else if (prop=="featurewidth") {
    return (double) featureheight;
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
    if ((!camera) && (!stopped))
      return (double) pVideoCapture->get(cv::CAP_PROP_POS_MSEC);
    else
      return (double) -1;
  }
  else if (prop=="minArea") {
    return (double) minArea;
  }   else  if ((prop == "DetectShadows") && (!stopped) ) {
    if (knnMethod)
      return (double) pBackgroundSubtractor->getDetectShadows();
    else
      return (double) -1;
  }
  else if ((prop == "Dist2Threshold")&& (!stopped)  ) {
    if (knnMethod)
      return (double) pBackgroundSubtractor->getDist2Threshold();
    else
      return (double) -1;
  }
  else if ((prop == "History")&& (!stopped) ) {
    if (knnMethod)
      return (double) pBackgroundSubtractor->getHistory();
    else
      return (double) pBackgroundThresholder->getHistory();
  }
  else if ((prop == "kNNSamples")&& (!stopped) ) {
    if (knnMethod) 
      return (double) pBackgroundSubtractor->getkNNSamples();
    else
      return -1;
  }
  else if ((prop == "NSamples")&& (!stopped) ) {
    if (knnMethod) 
      return (double) pBackgroundSubtractor->getNSamples();
    else
      return (double) -1;
  }
  else if ((prop == "ShadowThreshold")&& (!stopped) ){
    if (knnMethod)
      return (double) pBackgroundSubtractor->getShadowThreshold();
    else
      return (double) -1;
  }
  else if ((prop == "ShadowValue")&& (!stopped)) {
    if  (knnMethod)
      return (double) pBackgroundSubtractor->getShadowValue();
    else
      return (double) -1;
  }
  else if ((prop == "nskip")&& (!stopped)) {
    if  (knnMethod)
      return (double) -1;
    else
      return (double) pBackgroundThresholder->getNSkip();
  }
  else if ((prop == "FrameWidth") && (!stopped)) {
    double width;
    if (!camera) 
      width =  (double) pVideoCapture->get(cv::CAP_PROP_FRAME_WIDTH);
    else  {
      Size frameSize = pVideoSaver->getFrameSize();
      width = (double) frameSize.width;
    }
    if (resizeif)
      width = round(width*resizescale);
    return width;
  }
  else if ((prop == "FrameHeight") && (!stopped)) {
    double height;
    if (!camera) 
      height = (double) pVideoCapture->get(cv::CAP_PROP_FRAME_HEIGHT);
    else  {
      Size frameSize = pVideoSaver->getFrameSize();
      height =(double) frameSize.height;
    }
    if (resizeif)
      height = round(height*resizescale);
    return height;
  }
  else if ((prop == "FPS")  && (!stopped)){
    if (!camera) 
      return (double) pVideoCapture->get(cv::CAP_PROP_FPS);
    else 
      return (double) pVideoSaver->getFPS();
  }
  else if ((prop == "FrameCount")  && (!stopped)){
    if (!camera) 
      return (double) pVideoCapture->get(cv::CAP_PROP_FRAME_COUNT);
    else
      return (double) 0;
  } else if ((prop == "PosFrames") && (!stopped)) {
    if (!camera) 
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
//     ("Brightness",    cv::CAP_PROP_BRIGHTNESS)     //!< Brightness of the image (only for cameras).
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
  waitThread();
  for (int i=0;i<scale.size();i++) {
    m_Scale[i] = scale[i];
  }
  initialize();
};

/****************************************************************************************/
vector< float> VideoHandler::getScale() {
  return m_Scale;
};

/****************************************************************************************/
void VideoHandler::resetBkg() {
  if (!stopped) {
    waitThread();
    if (knnMethod) {
      pBackgroundSubtractor->clear(); // similar background anyway. do not reset 
    }
    else {
      pBackgroundThresholder->clear();
    }
    initialize();
  }
};



//*******************************************************************************//
//* Thresholder 
//*******************************************************************************//

BackgroundThresholder::BackgroundThresholder() {
  m_history = 250;
  m_nskip = 5;
  m_threstype = THRESH_BINARY_INV;
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
void BackgroundThresholder::setNSkip(int value) {
   m_nskip = value;
}
int BackgroundThresholder::getNSkip() {
   return m_nskip;
}

void BackgroundThresholder::setInverted(bool invertedif) {
  if (invertedif)
    m_threstype = THRESH_BINARY;
  else
    m_threstype = THRESH_BINARY_INV;
}

void BackgroundThresholder::apply(cv::Mat frame,cv::Mat * bwimg) {

  if (frame.size().width==0)
    return;

  cv::Mat floatframe;
  frame.convertTo(floatframe,CV_32FC1);
  
  // subtract background

  if (m_istep>MINBACKTHRESSTEP)  { // has to be a rough estimate of the mean already
    cv::Mat dframe;  
    dframe = floatframe - m_meanImage;
    dframe.convertTo(dframe,CV_8UC1);

    // apply threshold
    if (m_istep<MAXOTSU)  {
      m_thres = threshold(dframe,*bwimg,m_thres,255,THRESH_OTSU + m_threstype);
    } else {
      threshold(dframe,*bwimg,m_thres,255, m_threstype);
    }
  } else
    *bwimg = Mat::zeros(frame.size(),CV_8UC1);
  
  
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
  vid = "/data/videos/longterm/longterm8.avi";
  //VideoHandler vh(vid,false);
  VideoHandler vh(0,"/home/malte/test_VeideoHandler.avi",false);
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

  vector<Segment> segm;
  cv::Mat frame, oframe, bwimg;
  namedWindow("bwimg",WINDOW_AUTOSIZE);
  for (int i=0; i<100; i++) {
    vh.step(&segm,&oframe,&frame,&bwimg);

    imshow("bwimg",bwimg);
    waitKey(100);
    cout << i << ": Segment Size" <<  segm.size() << "\n";
  }
  vh.stop();
  vh.start();
  vh.stop();

  VideoHandler vh2(vid,true);

  for (int i=0; i<25; i++) {
    vh2.step(&segm,&oframe,&frame,&bwimg);
    waitKey(100);
    cout << i << ": Segment Size" << segm.size() << "\n";
  }
  vh2.stop();
  vh2.start();
  vh2.stop();
    
  return 0;
}


