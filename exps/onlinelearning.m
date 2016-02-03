
LOAD = 0;
COMPUTE =1;

if LOAD || ~exist('ft','var')
  videoFile = '/data/videos/onlinelearning/test.avi';
  opts = [];
  opts.detector.inverted = 1;
  opts.nfish = 3;
  opts.stmif = 1;
  opts.stimulus.screen = 1;
  opts.stimulus.presenter = 'FishStimulusPresenterOnlineLearning';
  opts.fishwidth = 30;
  opts.fishlength = 150;
  ft = FishTracker({0,''},opts);
end


if COMPUTE

  if ~exist('sbbox','var')
    %sbbox = ft.calibrateStimulusScreen();
    sbbox = [150,53,1580,1256];
  end
  ft.setOpts('stimulus.screenBoundingBox',sbbox,'display.displayEveryNFrame',100);
  ft.setDisplay(0);  
  
  
  ft.track();

  
  

end

