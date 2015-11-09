#include "FlyCapture2.h"
#include "SaveVideoClass.h"

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <boost/algorithm/string.hpp>

#include <iostream>
#include <sstream>
#include <string>   // for strings
#include <ncurses.h>

#define MINDELAY 250  // waiting time in msec between video frames for the display

using namespace std;
using namespace FlyCapture2;

//using namespace cv;

//#include "Profiler.h"


static void help()
{
    cout
        << "------------------------------------------------------------------------------" << endl
        << "Capturing of videos from the Grashopper. Codecs: X264 MPEG DIVX               "       << endl
        << "if WAITMSEC=0 no display is shown. Otherwise to number of msec to wait. "       << endl
        << "Usage:"                                                                         << endl
        << "         ./savevideo inputvideoName [CODEC WAITMSEC] "                          << endl
        << "Example: ./savevideo test.avi X264 500"                                         << endl
        << "------------------------------------------------------------------------------" << endl
        << endl;
};


int selectCameras(VideoSaver ***pppVideoSavers, unsigned int * pNumSelected)
  {
    BusManager busMgr;
    Error error;
    unsigned int numCameras;
    PGRGuid guid;


    error = busMgr.GetNumOfCameras(&numCameras);
    if (error != PGRERROR_OK)
      {
        error.PrintErrorTrace();
        return -1;
      }

    if (numCameras==0) 
      {
	cout << "Cannot find a camera..." << endl;
	return -1;
      }

    // select one or both cameras
    unsigned int numSelected=0;
    unsigned int selected[numCameras];

    initscr();	
    for ( unsigned int i = 0; i < numCameras; i++)
      {
        error = busMgr.GetCameraFromIndex( i, &guid );
        if (error != PGRERROR_OK)
        {
	  error.PrintErrorTrace();
	  return -1;
        }

	// Connect the camera
	Camera camera;
	error = camera.Connect( &guid);
	if ( error != PGRERROR_OK )
	  {
	    cout << "Failed to connect to camera" << endl;     
	    return -1;
	  }

	// Get the camera info and print it out
	CameraInfo camInfo;
	stringstream ss;

	error = camera.GetCameraInfo( &camInfo );
	if ( error != PGRERROR_OK )
	  {
	    cout << "Failed to get camera info from camera" << endl;     
	    return -1;
	  }
        ss  << camInfo.vendorName << " " << camInfo.modelName << " " << camInfo.serialNumber << " " 
	    << camInfo.sensorResolution;

	string tmpstr = string(ss.str()); // copy

	if (numCameras>1) 
	  {
	    int key= 0;
	    clear();
	    printw("** Should this camera be used for recording?\n");
	    printw("%c: %s \n y or n :", char(65+i),tmpstr.c_str());
	    refresh();
	    key = getch();
	    if (key=='y') 
	      {
		selected[numSelected++] = i;
		cout << " selected." << endl;
	      }
	    else
	      {
		cout << " NOT selected." << endl;
	      }
	  } else  
	  {
	    printw("Selected camera %s for recording.\n",tmpstr.c_str());
	    selected[numSelected++] = i;
	  }


	camera.Disconnect();

      }
    endwin();	 // % end the curses mode

    if (!numSelected) 
      {
	cout << "need at least one camera selected.!" <<endl;
	return -1;
      }

    // construct the objects and init
    (*pppVideoSavers) =  new VideoSaver*[numSelected];    
    for ( unsigned int i = 0; i < numSelected; i++) 
      {
        error = busMgr.GetCameraFromIndex( selected[i], &guid );
        if (error != PGRERROR_OK)
        {
	  error.PrintErrorTrace();
	  return -1;
        }

	(*pppVideoSavers)[i] = new VideoSaver();

	if ((*pppVideoSavers)[i]->init(guid)!=0)
	  return -1;
      }

    (*pNumSelected) = numSelected;
    return 0;
  }


int main(int argc, const char *argv[]) 
  {

    //parse inputs
    int waitMsecs;
    if ((argc ==1) || (argc >4)) 
      {
	help(); 
	cout << "Not enough parameters" << endl;
	return -1;
      } else if (argc == 4) {
      if (!(istringstream(argv[3]) >> waitMsecs)) {
	help(); 
	cout << "Third parameter needs to be an int" << endl;
	return -1;  
      }
    } else {
      waitMsecs = MINDELAY;
    }

    std::string fourcc;
    if (argc==2) {
      fourcc = string("X264");
    } else {
      fourcc = argv[2];
    }

    boost::to_upper(fourcc);
    
    if (fourcc.length()!=4) {
      help(); 
      cout << "Second parameter needs to be an FOURCC code" << endl;
      return -1;  
    }

    const string fname = argv[1];

    unsigned int numCameras;

    // init 
    VideoSaver **ppVideoSavers;

    if (selectCameras(&ppVideoSavers,&numCameras)!=0) 
      return -1;

    // start capturing
    for ( unsigned int i = 0; i < numCameras; i++)
      {
	string local_fname= std::string(fname);
	if (numCameras>1) {
	  local_fname = char(65 + i) + local_fname;
	}
	if (ppVideoSavers[i]->startCaptureAndWrite(local_fname,fourcc)!=0)
	  return -1;
      }

    // wait and display loop
    char key = 0;
    int waitAmount  = waitMsecs > MINDELAY? waitMsecs: MINDELAY;
    cv::Mat frame;
    string imname;

    while(key != 'q') 
      {

	key = cv::waitKey(waitAmount); 
	
	if (waitMsecs>0) {
	  //usleep(delayFound*1000);
	  for ( unsigned int i = 0; i < numCameras; i++)
	    {
	      //plot
	      imname[0]= char(i+65);
	      cv::Mat frame;
	      double timeStamp;
	      int frameNumber;
	      ppVideoSavers[i]->getFrame(&frame,&timeStamp,&frameNumber);
		      
	      cv::imshow(imname.c_str(), frame);
	      cout << char(i+65) << ": Current Frame rate " <<  ppVideoSavers[i]->getWritingFPS() << "Hz  | " ;
	    };
	  cout << "\r" ;
	}
	
	for ( unsigned int i = 0; i < numCameras; i++)
	  {
	    if (ppVideoSavers[i]->isFinished())
	      key = 'q'; // quit
	  }    
	
      }

    
    for ( unsigned int i = 0; i < numCameras; i++)
      ppVideoSavers[i]->close();

    return 0;

  }
  

