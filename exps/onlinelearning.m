
LOAD = 1
COMPUTE =1

if LOAD && ~exist('xyT','var')
  videoFile = ''; %'/data/videos/onlinelearning/test.avi';
  opts = [];
  opts.detector.inverted = 1;
  opts.nbody = 4;
  opts.stmif = 1;
  opts.stimulus.screen = 1;

  opts.stimulus.presenter = 'xy.stimulus.PresenterOnlineLearningCue';

  xyT = xy.Tracker({0,videoFile},opts);
end


if COMPUTE

  if ~exist('sbbox','var')
    %sbbox = xyT.calibrateStimulusScreen();
    sbbox = [81,83,1570,1231];
  end

  opts = [];
  opts.avgVelocity = 6;
  opts.maxVelocity = 80;
  
  
  opts.stimulus.screenBoundingBox = sbbox;
  opts.display.displayEveryNFrame = 100;
  opts.stimulus.stmSize = xyT.bodylength;
  opts.stimulus.usePredIdentityId = false;

  opts.stimulus.defaultColor = [1,1,1];
  opts.stimulus.borderWidth = 0.05;

  opts.stimulus.gapTime = 0;

  
  xyT.setOpts(opts);

  xyT.setDisplay(0);  
  
  xyT.track();
  %xyT.save();

  

end


