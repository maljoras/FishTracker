
LOAD = 1
COMPUTE =1

if LOAD && ~exist('ft','var')
  videoFile = ''; %'/data/videos/onlinelearning/test.avi';
  opts = [];
  opts.detector.inverted = 1;
  opts.nfish = 4;
  opts.stmif = 1;
  opts.stimulus.screen = 1;
  %opts.stimulus.screenBoundingBox = sbbox;
  
  opts.stimulus.presenter = 'fish.stimulus.PresenterOnlineLearningCue';



  ft = fish.Tracker({0,videoFile},opts);
end


if COMPUTE

  if ~exist('sbbox','var')
    %sbbox = ft.calibrateStimulusScreen();
    sbbox = [81,83,1570,1231];
  end
  ft.setOpts('stimulus.screenBoundingBox',sbbox,'display.displayEveryNFrame',100,'stimulus.stmSize',100,'stimulus.stmBkgType','texture','stimulus.stmLambda',0,'stimulus.stmTime',15,'stimulus.gapTime',1,'stimulus.midline',0.5,'stimulus.usePredFishId',false);
  
  ft.setDisplay(1);  

  
  ft.track();
  %ft.save();

  

end


