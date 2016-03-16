# fish.Tracker
Online tracking of fish in real time. 


FishTracker is a modular Matlab platform form fast tracking of fish for group learning experiments. 
It tracks fish in realtime and performs a online classification of the fish indentity to ensure that fish identity 
is not lost for long tracking experiments. 

##Installation

For the OpenCV version one needs to install [OpenCV](http:///www.opencv.org) version >=3.0. For the Matlab/OpenCV  functionality one needs to install the excellent [mexopencv](https://github.com/kyamagu/mexopencv) project. Note the FishTracker also uses [networkComponents](http://www.mathworks.com/matlabcentral/fileexchange/42040-find-network-components) and parts of the project [munkres-cpp](https://github.com/kaajo/munkres-cpp). 

One need to add the Matlab directory in the provided Makefile and compile with  
~~~~
$ make
~~~~

##Usage 
In MATLAB your need to add the path where the +fish package folder is located. Then, try (without arguments) to get some documentation about the parameters: 
~~~~
>> fish.Tracker;  
~~~~

To track a video file (with uidialog and automatic detection of the number if fish) :
~~~~
>> ft = fish.Tracker([]);  
>> ft.track();   
>> ft.plot(); % make some result plots  
>> ft.save(); % save the ft object and all results  
~~~~

To track video file 'myvideo.avi' having 3 fish write (use 'nfish',-1 for GUI selection, default is 'nfish', [] for auto selection of the number of fish)
~~~~
>> ft = fish.Tracker('myvideo.avi','nfish',3);  
~~~~

Set higher level of display and track first 20 seconds  
~~~~
>> ft.setDisplay(3)  
>> ft.track([0,20])  
~~~~

![Tracking screenshot](https://github.com/maljoras/FishTracker/blob/master/pics/track.png)

Or turn off the display for fastest tracking. 
~~~~
>> ft.setDisplay(0)
>> ft.track()
~~~~

For further analsys the results of the tracking process are saved in the ft.res field which can be obtaioned by
~~~~
>> res = ft.getTrackingResults();
~~~~

res.pos is a [nFrames x [x,y] x nFish] array of the position of the fish with unique ID.
res.tracks has a number of fields (controlled be the property ft.saveFields and changed by the methods ft.addSaveFields/ft.removeSaveField).   
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
                       fishId: [2101x5 double]  
~~~~

Each of the field as the dimensions [nFrames x nFish x nDims] where nDims are additional dimensions dependent on the field. For instance, to plot the x-velocity for each from and fish after tracking:

~~~~
>> res = ft.getTrackingResults();  
>> vx = res.tracks.velocity(:,:,1);  
>> plot(res.t,vx);  
~~~~

To plot the traces
~~~~
>> ft.plotTrace()
~~~~

![Trace](https://github.com/maljoras/FishTracker/blob/master/pics/trace.jpg)


Or in 3d (the centroid positions versus time), one can do the following:
~~~~
>> res = ft.getTrackingResults();
>> pos = ft.interpolateInvisible('pos',5); % interpolate lost detections and smooth Gaussian with std=5 frames
>> plot3(res.t,squeeze(p(:,1,:)),squeeze(p(:,2,:)));
>> xlabel('Time [s]'); ylabel('x-position [px]');zlabel('y-position [px]');
~~~

![Trace](https://github.com/maljoras/FishTracker/blob/master/pics/trace3d.jpg)

Further examples can be found in the 'exps' directory.

## Video Grabbing

For Video grabbing and online tracking currently only [PointGrey](https://www.ptgrey.com) cameras are supported via the FlyCaptureSDK library. After installation and compilation (edit the Makefile to add the directories) fish can be tracked from a camera if camera index CamIdx and simultanously saved to raw video file 'myfile.avi' 

~~~~
>> ft = fish.Tracker({0,'myfile.avi'});   
>> ft.track(); % tracks indefintely. Close trackiung window for stop or set end time  
>> ft.save();  
>> clear ft; % to stop the background video recording  
~~~~

Tracking stops when the tracking display window is closed and results are saved. Note that the video is continously recorded in the background until the ft object is cleaerd. 

## Stimulation

Using the Matlab PsychToolbox it is possible to do online feedback experiments. To peform a new experiment one has to derive a class from fish.stimulus.Presenter and overload (at least) the "step" function. See FishStimulusPresenterCue for an example

To perform the experiment:  

~~~~
>> opts = [];   
>> opts.stimulus.presenter  = 'fish.stimulus.Presenter';  
>> opts.stmif = 1;  
>> opts.nfish = 3;  
>> opts.detector.inverted = 1; % tracking in IR  
>> opts.stimulus.screen = 1; % X-window screen number  
>> opts.fishwidth = 30;  
>> opts.fishlength = 150;  
>>   
>> ft = fish.Tracker({0,'savevideo.avi'},opts);  
>> ft.setDisplay(0);  % to avoid displaying delays  
>> ft.track(); % tracks until fish.stimulus.Presenter/isFinished returns true.  
>> ft.save();  
>> clear ft; % stop background avi recording  
~~~~
