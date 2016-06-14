
LOAD = 1
COMPUTE =1

if LOAD && ~exist('ft','var')
  videoFile = ''; %'/data/videos/onlinelearning/test.avi';
  opts = [];
  opts.detector.inverted = 1;
  opts.nfish = 4;
  opts.stmif = 1;
  opts.stimulus.screen = 1;
  %opts.stimulus.presenter = 'fish.stimulus.PresenterOnlineLearningCue';
  opts.stimulus.presenter = 'fish.stimulus.Presenter';
  
  opts.fishwidth = 25;
  opts.fishlength =  90;
  
  
  ft = fish.Tracker({0,videoFile},opts);
end


if COMPUTE

  if ~exist('sbbox','var')
    %sbbox = ft.calibrateStimulusScreen();
    sbbox = [81,83,1570,1231];
  end
  opts = [];
  opts.avgVelocity = 6;
  opts.maxVelocity = 80;
  
  
  opts.stimulus.screenBoundingBox = sbbox;
  opts.display.displayEveryNFrame = 100;
  opts.stimulus.stmSize = ft.fishlength;
  opts.stimulus.usePredFishId = false;

  opts.stimulus.defaultColor = [1,1,1];
  opts.stimulus.borderWidth = 0.05;

  opts.stimulus.gapTime = 0;

  
  ft.setOpts(opts);

  ft.setDisplay(0);  
  
  ft.track();
  ft.save();

  

end


