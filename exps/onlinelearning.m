
LOAD = 1
COMPUTE =1

if LOAD && ~exist('ft','var')
  videoFile = ''; %'/data/videos/onlinelearning/test.avi';
  opts = [];
  opts.detector.inverted = 1;
  opts.nanimals = 4;
  opts.stmif = 1;
  opts.stimulus.screen = 1;

  opts.stimulus.presenter = 'xy.stimulus.PresenterOnlineLearningCue';

  ft = xy.Tracker({0,videoFile},opts);
end


if COMPUTE

  if ~exist('sbbox','var')
    %sbbox = ft.calibrateStimulusScreen();
    sbbox = [81,83,1570,1231];
  end
<<<<<<< HEAD
  ft.setOpts('stimulus.screenBoundingBox',sbbox,'display.displayEveryNFrame',100,'stimulus.stmSize',100,'stimulus.stmBkgType','texture','stimulus.stmLambda',0,'stimulus.stmTime',15,'stimulus.gapTime',1,'stimulus.midline',0.5,'stimulus.usePredIdentityId',false);
  
  ft.setDisplay(1);  
=======
  opts = [];
  opts.avgVelocity = 6;
  opts.maxVelocity = 80;
  
  
  opts.stimulus.screenBoundingBox = sbbox;
  opts.display.displayEveryNFrame = 100;
  opts.stimulus.stmSize = ft.fishlength;
  opts.stimulus.usePredIdentityId = false;

  opts.stimulus.defaultColor = [1,1,1];
  opts.stimulus.borderWidth = 0.05;

  opts.stimulus.gapTime = 0;
>>>>>>> development

  
  ft.setOpts(opts);

  ft.setDisplay(0);  
  
  ft.track();
  %ft.save();

  

end


