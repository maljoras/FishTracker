
LOAD = 1
COMPUTE =1

if LOAD && ~exist('T','var')
  videoFile = ''; %'/data/videos/onlinelearning/test.avi';
  opts = [];
  opts.detector.inverted = 1;
  opts.nindiv = 4;
  opts.stmif = 1;
  opts.stimulus.screen = 1;

  opts.stimulus.presenter = 'xy.stimulus.PresenterOnlineLearningCue';

  T = xy.Tracker({0,videoFile},opts);
end


if COMPUTE

  if ~exist('sbbox','var')
    %sbbox = T.calibrateStimulusScreen();
    sbbox = [81,83,1570,1231];
  end

  opts = [];
  opts.avgVelocity = 6;
  opts.maxVelocity = 80;
  
  
  opts.stimulus.screenBoundingBox = sbbox;
  opts.display.displayEveryNFrame = 100;
  opts.stimulus.stmSize = T.bodylength;
  opts.stimulus.usePredIdentityId = false;

  opts.stimulus.defaultColor = [1,1,1];
  opts.stimulus.borderWidth = 0.05;

  opts.stimulus.gapTime = 0;

  
  T.setOpts(opts);

  T.setDisplay(0);  
  
  T.track();
  %T.save();

  

end


