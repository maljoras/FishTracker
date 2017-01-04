# xy.Tracker
Online tracking of groups of small animals in a planar environment in real time. See [bioRxiv](http://dx.doi.org/10.1101/071308) for details on the tracking methods. If you use xyTracker, please cite the  [bioRxiv](http://dx.doi.org/10.1101/071308) paper. 


xyTracker is a modular Matlab platform for real-time tracking of
animals for group learning experiments in planar environments. 
It tracks animals in real time and performs an online classification of the animal identity to ensure that animal identity 
is not lost for long tracking experiments. The name "xy" derives from Chinese "XiaoYu" (小鱼).

##Installation

The tracking system includes 3 versions (in order of increasing performance): A purely matlab based version, an OpenCV/matlab version, and an optimized MEX/C++/OpenCV version. All versions have the same matlab interface. The system automatically chooses the fastest available version depending what additional packages are installed on the system.  **However, the purely matlab version is depreciated since it is slow and not well maintained**. Please install the OpenCV/Mex version as described below.

Note further, that the xy.Tracker also uses (and includes a copy of) [networkComponents](http://www.mathworks.com/matlabcentral/fileexchange/42040-find-network-components) and parts of the project [munkres-cpp](https://github.com/kaajo/munkres-cpp). 

###Installation instructions

For the matlab-based version, only matlab is needed (no additional matlab toolbox license should be necessary; the vision toolbox is used for a different plotting interface, but not necessary). 

For fast tracking, install [OpenCV](http:///www.opencv.org) version >=3.0. and the excellent [mexopencv](https://github.com/kyamagu/mexopencv) project. **Please follow the installation instructions of the [mexopencv](https://github.com/kyamagu/mexopencv) toolbox to  install OpenCV and mexopencv**.  

Optionally, for grabbing from ptGray cameras the FlyCapture SDK is needed. One can download it from the [ptGrey website](http://www.ptgrey.com). 

After installing these prerequisites, one needs to compile the xyTracker source code. 

#### Windows
Call the provided matlab file for compilation. Thus, start matlab, and run the file
~~~~
>> installHelpersOnWindows;
~~~~
Make sure that the mexopencv directory is in the matlab path.  This will compile and install the source code if OpenCV and mexopencv are properly installed. If FlyCaptureSDK functionality should be compiled as well, the path settings in installHelpersOnWindows might have to be adjusted. 

#### Linux
For compilation, one needs to specify the Matlab path and the path to mexopencv and the path to FlycaptureSDK (if available) to compile with  
~~~~
$ make  MATLABDIR=/my/path/to/MATLAB/ FLYCAPINCLUDEDIR=/usr/include/flycapture MEXOPENCVDIR=/mypath/to/mexopencv
~~~~

To avoid library incompatibility with Matlab's packaged libraries one needs to preload the OpenCV and other libraries. With Linux, eg., put an alias into the .bashrc file (with possible adjustements of the paths):

~~~~
export PRELOAD_LIBS=/usr/lib64/libstdc++.so.6:/usr/lib64/libtiff.so.5:/usr/lib/libflycapture.so:/usr/lib64/libglibmm-2.4.so:/usr/lib64/libglib-2.0.so.0:/usr/lib64/libsigc-2.0.so.0:/usr/lib/libflycapture.so.2:`ls /usr/local/lib/libopencv*.so |xargs|tr ' ' ':'`
alias matlab='LD_PRELOAD=$PRELOAD_LIBS /opt/MATLAB/R2014b/bin/matlab '
~~~~ 

### Stimulation environment

Finally, for stimulus presenting functionality, the [PsychToolbox](http://psychtoolbox.org/) has to be installed and added to the matlab path. 

##Test
One can test the installation by running (from the directory where +xy directory can be seen)
~~~~
>> xy.Tracker.runSimpleTest()
~~~~

##Usage 
In MATLAB your need to add the path where the +xy package folder is located. Then, try (without arguments) to get some documentation about the parameters: 
~~~~
>> xy.Tracker;  
~~~~

To track a video file (with ui-dialog and automatic detection of the
number of bodies :
~~~~
>> T = xy.Tracker([]);  
>> T.track();   
>> T.plot(); % make some result plots  
>> T.save(); % save the T object and all results  
~~~~

To track video file 'myvideo.avi' having 3 animals, write (use nindiv=-1
for GUI selection, default is nindiv=[] for auto selection of the
number of animals) and known approximate length (e.g. 100 px) and width
(e.g. 30 px), do 
~~~~
>> T = xy.Tracker('myvideo.avi','nindiv',3,'bodywidth',30,'bodylength',100);  
~~~~

Set higher level of display and track first 20 seconds  
~~~~
>> T.setDisplay(3)  
>> T.track([0,20])  
~~~~

![Tracking screenshot](https://github.com/maljoras/xyTracker/blob/master/pics/track.png)

Or turn off the display for fastest tracking. 
~~~~
>> T.setDisplay(0)  
>> T.track()
~~~~

For further analysis the results of the tracking process are saved in the T.res field which can be obtained by
~~~~
>> res = T.getTrackingResults();
~~~~

res.pos is a [nFrames x [x,y] x nindiv] array of the position of the animal with unique ID.
res.tracks has a number of fields (controlled be the property T.saveFields and changed by the methods T.addSaveFields/T.removeSaveField).   
For example:

~~~~
>> res.tracks
   
   ans =  
                           id: [2101x5 double]  
                            t: [2101x5 double]  
                     centroid: [2101x5x2 double]  
                    classProb: [2101x5x5 double]  
                         bbox: [2101x5x4 double]  
               assignmentCost: [2101x5 double]  
                     velocity: [2101x5x2 double]  
    consecutiveInvisibleCount: [2101x5 double]  
          segment_Orientation: [2101x5 double]  
                   centerLine: [4-D double]  
                    thickness: [2101x5x7 double]   
      segment_MinorAxisLength: [2101x5 double]  
      segment_MajorAxisLength: [2101x5 double]  
             segment_reversed: [2101x5 double]  
                   identityId: [2101x5 double]  
~~~~

Each of the field as the dimensions [nFrames x nindiv x nDims] where
nDims are additional dimensions dependent on the field. For instance,
to plot the x-velocity for each tracked body after tracking:

~~~~
>> res = T.getTrackingResults([0,20]);  
>> v = T.getResField(res,'velocity',1); % 1 means delete "invisible" frames 
>> vx = squeeze(v(:,:,1)); 
>> plot(res.t,vx);  
~~~~

To plot the traces
~~~~
>> T.plotTrace()
~~~~

![Trace](https://github.com/maljoras/xyTracker/blob/master/pics/trace.jpg)


Or in 3D (the centroid positions versus time), one can do the following:
~~~~
>> res = T.getTrackingResults();
>> pos = T.interpolateInvisible(res,'pos',5); % interpolate lost detections and Boxcar smooth with n=5 frames
>> plot3(res.t,squeeze(p(:,1,:)),squeeze(p(:,2,:)));
>> xlabel('Time [s]'); ylabel('x-position [px]');zlabel('y-position [px]');
~~~~

![Trace](https://github.com/maljoras/xyTracker/blob/master/pics/trace3d.jpg)

Or a more fancy plot of the "centerline" the body-axis of the bodies for each frame. Moments of large amount of "bending" of the fish's body are plotted in red. The position of the head is marked by a circle. 

~~~~
>> T.plotCenterLine([10,20]) % 10 to 20 seconds for all identities   
~~~~
![Trace](https://github.com/maljoras/xyTracker/blob/master/pics/centerline.jpg)



Further examples can be found in the 'exps' and 'figs' directories.

## Video Grabbing

For Video grabbing and online tracking currently only
[PointGrey](https://www.ptgrey.com) cameras are supported via the
FlyCaptureSDK library. After installation and compilation (see above) animals can be tracked from a camera with
camera index CamIdx and simultaneously saved to raw video file
'myfile.avi' (camera index and ROI can be set with the flycap2 tool)

~~~~
>> T = xy.Tracker({0,'myfile.avi'});   
>> T.track(); % tracks indefintely. Close tracking window for stop or set end time  
>> T.save();  
>> clear ft; % to stop the background video recording  
~~~~

Tracking stops when the tracking display window is closed and results
are saved. Note that the video is continuously recorded in the
background until the ft object is cleared.

### Grabbing with other cameras

Although xyTracker is currently only tested with ptGrey cameras, it can also be used with other cameras through the OpenCV VideoCapture functionality for video grabbing. That is, try
for real-time tracking of 4 animals (with a camera with index 0) 

~~~~
>> T = xy.Tracker('0','nindiv',4);
>> T.setDisplay(1); % to test the tracking; note,however, that plotting will introduce a delay
>> T.track();
>> clear all; % will release the video hardware
~~~~

Note that the camera index is given as a string in matlab. Cameras indices 0-9 are supported. However, in contrast to the FlyCapture method (see above), the grabbed video is currently not encoded and saved in the background, when using this VideoCapture method.  


## Stimulation

Using the Matlab PsychToolbox it is possible to do online feedback experiments. To perform a new experiment one has to derive a class from xy.stimulus.Presenter and overload (at least) the "stepStimulus" method. See xy.stimulus.PresenterFlash for an example

To perform the experiment:  

~~~~
>> opts = [];   
>> opts.stimulus.presenter  = xy.stimulus.Presenter';  
>> opts.stmif = 1;  
>> opts.nindiv = 3;  
>> opts.detector.inverted = 1; % tracking in IR  
>> opts.stimulus.screen = 1; % X-window screen number  
>> opts.bodywidth = 30;  
>> opts.bodylength = 150;  % set approx. length and width if known to avoid auto-estimation
>>   
>> T = xy.Tracker({0,'savevideo.avi'},opts);  
>> T.setDisplay(0);  % to avoid displaying delays  
>> T.track(); % tracks until xy.stimulus.Presenter/isFinished returns true.  
>> T.save();  
>> clear T; % stop background avi recording  
~~~~

For debugging, stimulus presenter objects can be tested on
simulated animal tracks without the need for calling the
xy.Tracker.track() method:

~~~~
>> xy.helper.stimulusSimulator(xy.stimulus.Presenter')
~~~~


