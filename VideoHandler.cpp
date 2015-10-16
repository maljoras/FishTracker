

#include "VideoHandler.h"
#include <iostream>
#include <unistd.h>


using namespace std;
using namespace cv;

Segment::Segment() {}
Segment::~Segment() {};


VideoHandler::VideoHandler(const string fname)
{

  scaled = false;
  plotif = false;
  minContourSize = 5;
  nprobe = 5;
  m_Delta = 0;

  for (int i; i<3;i++)
    m_Scale.push_back(0.33333);

  pVideoCapture = new VideoCapture(fname);
  pBackgroundSubtractor =  cv::createBackgroundSubtractorKNN(250,400,false);
  _readNextFrameThread(); // read first frame;
};


VideoHandler::~VideoHandler()
{
  m_FrameMutex.lock(); // wait for lock
  m_FrameMutex.unlock();
};

void VideoHandler::setTime(double msec) {

  m_FrameMutex.lock(); // wait for lock
  pVideoCapture->set(cv::CAP_PROP_POS_MSEC,msec);
  //pBackgroundSubtractor->clear(); // similar background anyway
  m_FrameMutex.unlock();
  
  _readNextFrameThread(); // read first frame;

};


void VideoHandler::setScale(vector<float> scale, vector<float> delta) {
  m_Scale = scale;
  for (int i=0;i<scale.size();i++)
    m_Scale[i] = m_Scale[i]/255.;
  m_Delta = 0.;
  for (int i=0;i<delta.size();i++)
    m_Delta += delta[i];
};


void VideoHandler::getSegment(Segment * segm, vector<Point> inContour, Mat inBwImg, Mat inFrame, Mat inOFrame) {


  Mat mContours = Mat(inContour);
  segm->Bbox = boundingRect(mContours);
  segm->Image = Mat(inBwImg(segm->Bbox)>0);
  segm->FilledImage = Mat(inFrame(segm->Bbox));
  segm->FilledColImage = Mat(inOFrame(segm->Bbox));

  if ((inContour.size()>=5) && segm->Bbox.width>5 && segm->Bbox.height>5){
    RotatedRect minEllipse = fitEllipse(inContour);
    segm->Orientation = minEllipse.angle;
    segm->Center = minEllipse.center;
    segm->Size = minEllipse.size;

    // get the rotation and the center points
    Mat img;
    int nborder = 3;
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
	return;
      }

    threshold(L,L,minVal*minThres,255,THRESH_TOZERO_INV);


    vector<Point> locations;
    bool trans = false;
    if (img.size().height<img.size().width) {
      // transpose for a better mapping
      findNonZero(L.t()!=0,locations);
      trans = true;
    } else {
      findNonZero(L!=0,locations);
      trans = false;
    }

    vector<Point> probe(nprobe);
    for (int i=0; i<nprobe;i++ ) {
      probe[i] = locations[i*(locations.size()-1)/(nprobe-1)];
    }

    if (trans) // swap
      for (int i=0; i<nprobe;i++ ) {
	int y = probe[i].y;
	probe[i].y = probe[i].x;
	probe[i].x = y;
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
    }
    for (int i=0; i<nprobe;i++ ) {
      probe[i].x = probe[i].x-nborder;
      probe[i].y = probe[i].y-nborder;
    }
    
    segm->CenterLine = probe;
    segm->Thickness = thickness;

    Point2f v = Point2f(probe[0]-probe[1]);
    double d = norm(v);

    double a = v.x/d;
    double b = v.y/d;

    Size2f szin = Size2f(img.size());
    Size2f szout = Size2f(szin.height*fabs(b)+szin.width*fabs(a),szin.height*fabs(a)+szin.width*fabs(b));

    Point2f center = Point2f(probe[0]);
    Point2f offs = Point2f(szout.width-2*thickness[0],szout.height/2)-center;

    Mat T = (Mat_<float>(2,3) << a, b, (1-a)*center.x-b*center.y + offs.x, -b,a,b*center.x+(1-a)*center.y + offs.y );

    Mat rotImg;
    warpAffine(img,rotImg,T,szout,INTER_NEAREST,BORDER_CONSTANT,0);
    segm->RotImage = rotImg;

    // rotate the other Images
    Scalar mback = mean(segm->FilledImage,segm->Image==0);
    Scalar mbackcol = mean(segm->FilledColImage,segm->Image==0);
    warpAffine(segm->FilledImage,segm->RotFilledImage,T,szout,INTER_NEAREST,BORDER_CONSTANT,mback);
    warpAffine(segm->FilledColImage,segm->RotFilledColImage,T,szout,INTER_NEAREST,BORDER_CONSTANT,mbackcol);

  }

};


void VideoHandler::plotCurrentFrame() {

  namedWindow( "Frame", WINDOW_AUTOSIZE );
  imshow( "Frame", Frame);
  //waitKey(100);
}

void VideoHandler::step(){

  m_FrameMutex.lock(); // wait for lock

  // need to clone every thing ?
  OFrame = m_NextOFrame;
  Frame = m_NextFrame;
  BWImg = m_NextBWImg;
  Contours = m_NextContours;
  Segments = m_NextSegments;

  m_FrameMutex.unlock();

  // start a new thread
  boost::thread readNextFrameThread(&VideoHandler::_readNextFrameThread,this);
  //_readNextFrameThread();
  if (plotif)
    plotCurrentFrame();

};


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
	  channel[ii].convertTo(channel[ii], CV_32FC1);
	}
	frame =  m_Scale[2]*channel[0] + m_Scale[1]*channel[1] + m_Scale[0]*channel[2] + m_Delta;

	// substract mean
	frame = frame - cv::mean(frame);
	frame += 0.5; // between 0..1
	frame.convertTo(frame, CV_8UC1, 255);
    }
    else
      cvtColor(oframe,frame,cv::COLOR_BGR2GRAY);

    // backgound computation
    pBackgroundSubtractor->apply(frame, m_NextBWImg);

    // finding contours
    vector<Vec4i> hierarchy;
    m_NextContours.clear();
    findContours(m_NextBWImg,m_NextContours,hierarchy,cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE, Point(0, 0));


    // get the fish patches and segments
    vector<Segment> segms(m_NextContours.size());
    for( int i = 0; i< m_NextContours.size(); i++ )
      getSegment(&(segms[i]),m_NextContours[i],m_NextBWImg,frame,oframe);

    m_NextSegments = segms;
    m_NextFrame = frame;
    m_NextOFrame = oframe;

  } else {

    cout <<  "File read error";
    m_FrameMutex.unlock();
    return -1;
  }


  m_FrameMutex.unlock();
  return 0;
};


// // main
// int main() {

//   VideoHandler vh("/home/malte/Videos/5Zebrafish_nocover_22min.avi");
//   namedWindow( "Patch", WINDOW_AUTOSIZE );
//   for (int i=0; i<4; i++) {
//     vh.step();
//     cout << i << ": Segment Size" <<  vh.Segments.size() << "\n";
//   }

//   return 0;
// }
