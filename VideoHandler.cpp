

#include "VideoHandler.h"
#include <iostream>
#include <unistd.h>

using namespace std;
using namespace cv;
using namespace boost::posix_time;

Segment::Segment() {}
Segment::~Segment() {};


VideoHandler::VideoHandler(const string fname)
{

  scaled = false;
  plotif = false;
  colorfeature = false;

  nprobe = 7;
  m_Delta = 0;

  minWidth = 2 ;
  minExtent = 2;
  minArea = 2;
  maxArea = 10000;
  maxExtent = 10000;
  ngoodmsk = 0;
  featureheight = 100;
  featurewidth =  20;
  featureSize  = Size(featureheight,featurewidth); // transposed image
  vector<float> tmp(3);
  m_Scale = tmp;
  for (int i; i<3;i++)
    m_Scale[i] = 0.33333;

  pVideoCapture = new VideoCapture(fname);
  pBackgroundSubtractor =  cv::createBackgroundSubtractorKNN(250,400,false);

  _readNextFrameThread(); // directly fill to avoid time stuff
};


VideoHandler::~VideoHandler()
{
  m_FrameMutex.lock(); // wait for lock
  usleep(500);
  m_FrameMutex.unlock();
};

void VideoHandler::setTime(double msec) {

  m_FrameMutex.lock(); // wait for lock
  pVideoCapture->set(cv::CAP_PROP_POS_MSEC,msec);
  //pBackgroundSubtractor->clear(); // similar background anyway. do not reset
  m_FrameMutex.unlock();
  stopAndInitialize();
};



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
      // delete detection (see goodmsk)
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


 //    if (segm->Size.height + segm->Size.width > 100 ){
//       namedWindow("proj",WINDOW_AUTOSIZE);
//       Mat timg = segm->FilledImage;
//       Mat timg2;
//       getRectSubPix(timg,Size(300,300),localcenter,timg2,-1);

//       Point2f p1(150,150);
//       Size sz= segm->Size;
      
//       line(timg2,pvec*10+p1,p1,CV_RGB(255,255,145),1,8,0);
//       line(timg2,v*10 + p1,  p1,CV_RGB(255,255,245),2,8,0);
//       for (int i=0;i<nprobe;i++) {
// 	Point2f p2(probe[i]);
// 	circle(timg2,p2 + p1 - localcenter,5,CV_RGB(255-i*20,255-i*20,245-i*20),2,8,0);
//       }
//       circle(timg2,center + p1 - localcenter,5,CV_RGB(255,255,255),4,8,0);
//       ellipse(timg2,p1,(sz)/2,segm->Orientation,0,360,CV_RGB(255,255,145),2,8,0);

//       imshow("proj",timg2);


// //      waitKey(0); 
//     }

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
    meanStdDev(comy, temp, stdValue, comy>0);
    segm->bendingStdValue = stdValue.at<float>(0);

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

      //Point2f centerfish(szout.width - featureSize.width/2,szout.height/2);
      getRectSubPix(fishfeature,featureSize,centerfish,segm->FishFeatureRemap,-1);
      segm->FishFeatureRemap.convertTo(segm->FishFeatureRemap,CV_32FC1);
    }
    
    // if (segm->Size.height + segm->Size.width > 100 ){
    //   namedWindow("fish",WINDOW_AUTOSIZE);
    //   imshow("fish",segm->FishFeature);


      
    //   Mat timg6;
    //   Point2f p8(fishfeature.size());
    //   getRectSubPix(fishfeature,Size(300,300),p8/2,timg6,-1);
    //   namedWindow("fishfeature",WINDOW_AUTOSIZE);
    //   imshow("fishfeature",timg6);


    //   Mat timg4(segm->RotFilledImage);
    //   for (int i=0;i<nprobe;i++) {
    // 	circle(timg4, rotprobe[i],4,CV_RGB(0,255,245),1,8,0);
    //   }

    //   // for (int i=0;i<comy.cols;i++) {
    //   // 	circle(timg4, Point2f(i,comy.at<float>(i)+rotimg.cols/2.),1,CV_RGB(255,255,245),1,8,0);
    //   // }
    //   for (int i=0;i<comy.cols;i++) {
    // 	circle(timg4, Point2f(i,comy.at<float>(i)),1,CV_RGB(0,255,245),1,8,0);
    //   }
 
    //   Mat timg3;
    //   Point2f p5 = Point2f(szout);
    //   getRectSubPix(timg4,Size(300,300),p5/2,timg3,-1);

    //   namedWindow("rotimg",WINDOW_AUTOSIZE);
    //   imshow("rotimg",timg3);
    //   waitKey(0); 
    // }

  }
  
};


void VideoHandler::makeGoodMsk() {

  vector<bool> msk(Segments.size());
  ngoodmsk = 0;

  for (int i=0;i<Segments.size();i++){ 
    double  extent = Segments[i].MajorAxisLength + Segments[i].MinorAxisLength;
    if ((extent<minExtent) || (extent>maxExtent) || (Segments[i].Area<minArea)|| (Segments[i].Area>maxArea) || (Segments[i].MinorAxisLength<minWidth)) 
      msk[i] = false;
    else {
      msk[i] = true;
      ngoodmsk++;
    }
  }
  goodmsk = msk;
}


void VideoHandler::plotCurrentFrame() {

  namedWindow( "Frame", WINDOW_AUTOSIZE );
  imshow( "Frame", Frame);
  //waitKey(100);
}

void VideoHandler::startThread() {


  boost::thread readNextFrameThread(&VideoHandler::_readNextFrameThread,this);

}

void VideoHandler::stopAndInitialize() {

  m_FrameMutex.lock(); // wait for lock

  int iframe;
  iframe = pVideoCapture->get(cv::CAP_PROP_POS_FRAMES);
  if (iframe>0)
    iframe = iframe-1; // step one back;
  pVideoCapture->set(cv::CAP_PROP_POS_FRAMES,iframe);

  m_FrameMutex.unlock();
  // start a new thread
  startThread();
}


void VideoHandler::step(){

 

  m_FrameMutex.lock(); // wait for thread 
  OFrame = m_NextOFrame;
  Frame = m_NextFrame;
  BWImg = m_NextBWImg;
  m_FrameMutex.unlock();

  // start a new thread to load in the background
  startThread();
  
  // ptime thisTime,thatTime;
  // time_duration td;
  // thatTime =  microsec_clock::local_time();

  // // finding contours
  vector<Vec4i> hierarchy;
  findContours(BWImg,Contours,hierarchy,RETR_EXTERNAL, CHAIN_APPROX_SIMPLE, Point(0, 0));
    
  // get the fish patches and segments
  vector<Segment> segms(Contours.size());
  for( int i = 0; i< Contours.size(); i++ )
    getSegment(&(segms[i]),Contours[i],BWImg,Frame,OFrame);

  Segments = segms;
  
  // make mask
  makeGoodMsk();

  // thisTime = microsec_clock::local_time();
  // td = thisTime - thatTime;	
  // cout << td.total_milliseconds() << "\n";

  
  if (plotif)
    plotCurrentFrame();

};

void VideoHandler::setScale(vector<float> scale) {
  m_FrameMutex.lock(); 
  m_Scale = scale;
  m_FrameMutex.unlock(); 
  stopAndInitialize();
};
vector< float> VideoHandler::getScale() {
  return m_Scale;
};

int VideoHandler::set(const string prop, double value){

  if (prop=="scaled") {
    m_FrameMutex.lock();
    scaled = (bool) (value>0);
    m_FrameMutex.unlock();
    stopAndInitialize();
  }
  else  if (prop=="delta") {
    m_FrameMutex.lock();
    m_Delta = (float) value;
    m_FrameMutex.unlock();
    stopAndInitialize();
  }
  else if (prop=="featureheight") {
    m_FrameMutex.lock();
    featureheight = (int) value;
    featureSize.width = featureheight; // transposed
    m_FrameMutex.unlock();
    stopAndInitialize();
  }
  else if (prop=="featurewidth") {
    m_FrameMutex.lock();
    featureheight = (int) value;
    featureSize.height = featurewidth; // transposed
    m_FrameMutex.unlock();
    stopAndInitialize();
  }
  else if (prop=="plotif") {
    plotif = (bool) (value>0);
  }
  else if (prop=="colorfeature") {
    m_FrameMutex.lock();
    colorfeature = (bool) (value>0);
    m_FrameMutex.unlock();
    stopAndInitialize();
  }
  else if (prop=="nprobe") {
    m_FrameMutex.lock();
    nprobe = (int) value;
    m_FrameMutex.unlock();
    stopAndInitialize();
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
  }  else {
    cout <<  "ERROR: Could not set property " << prop <<"! \n";
    return -1;
  }
  return 0;
}

double VideoHandler::get(const string prop){

  if (prop=="scaled") {
    return (double) scaled;
  }
  else if (prop=="plotif") {
    return(double) plotif;
  }
  else if (prop=="delta") {
    return(double) m_Delta;
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
  else if (prop=="minArea") {
    return (double) minArea;
  }  else {
    cout <<  "ERROR: Could not read property " << prop <<"! \n";
    return -1;
  }
  return 0;
}

int VideoHandler::_readNextFrameThread()
{
  m_FrameMutex.lock();
  Mat oframe,frame;
  if (pVideoCapture->read(oframe)) {

    if (oframe.type()!=CV_8UC3) {
      m_FrameMutex.unlock();
      cout <<  "ERROR: Expect CV_8UC3 color movie";
      return -1;
    }

    if (scaled) {
	Mat channel[3];

	split(oframe, channel);
	for (int ii=0; ii<3; ii++) {
	  channel[ii].convertTo(channel[ii], CV_32FC1,1/255.);
	}
	frame =  m_Scale[2]*channel[0] + m_Scale[1]*channel[1] + m_Scale[0]*channel[2] + m_Delta;

	// substract mean
	frame = frame - cv::mean(frame);
	frame += 0.5; // between 0..1
	frame.convertTo(frame, CV_8UC1, 255);
    }
    else
      cvtColor(oframe,frame,cv::COLOR_BGR2GRAY);

    // backgound computation (THIS IS THE SLOWEST PART!)
    pBackgroundSubtractor->apply(frame, m_NextBWImg);
      
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
    

  } else {

    cout <<  "File read error";
    m_FrameMutex.unlock();
    return -1;
  }

  m_FrameMutex.unlock();

  return 0;
};


// // // main
int main() {

  VideoHandler vh("/home/malte/data/zebra/videos/test_wtih_bu_video1.avi");
  //namedWindow( "Patch", WINDOW_AUTOSIZE );
  vh.set("scaled",false);
  vh.set("colorfeature",true);
  for (int i=0; i<30; i++) {
    vh.step();
    waitKey(100);
    cout << i << ": Segment Size" <<  vh.Segments.size() << "\n";
  }

  return 0;
}
