
LOAD = 1
COMPUTE =1

if LOAD && ~exist('ft','var')
  videoFile = '/data/videos/onlinelearning/dotsInTheMiddleTest1.avi';
  opts = [];
  opts.detector.inverted = 1;
  
  opts.nfish = 5;

  opts.stmif = 1;
  opts.display.videoHandler = true;
  opts.stimulus.screen = 1;
  opts.stimulus.presenter = 'fish.stimulus.PresenterOnlineLearningCue';
  opts.detector.adjustThresScale = 1.05;
  
  opts.fishwidth = 25;
  opts.fishlength =  90;
  
  
  ft = fish.Tracker({0,videoFile},opts);
end


if COMPUTE

  if ~exist('sbbox','var')
    %sbbox = ft.calibrateStimulusScreen();
    sbbox = [88,63,1584,1247];
  end
  opts = [];
  opts.avgVelocity = 5;
  opts.maxVelocity = 20;
    

  opts.display.displayEveryNFrame = 100;
  
  ostm = [];
  ostm.screenBoundingBox = sbbox;  
  ostm.usePredFishId = false;

  ostm.stmCol= parula(ft.nfish);
  ostm.stmLambda= 0.0;
  
  ostm.stmSize = ones(ft.nfish,1)*ft.fishlength.*(1+rand(ft.nfish,1));  

  ostm.stmType = 'fishtextures';
  ostm.stmSizeFactor = 1.2;
  ostm.stmShift = ft.fishlength; % in px of fish.tracker frame
  ostm.stmShiftOri = 0;
  ostm.fishVelThres = 2;
  
  ostm.funLRStm = @(fishIds) fishIds;
  ostm.lrif = 1;
  ostm.lrSwitchProb = 1;
  ostm.midline = 0;
  ostm.colBackground = [0,0,0];
  
  ostm.colBorder = [1,1,1];
  ostm.borderWidth = 0.05;

  ostm.nRound = 1000;

  ostm.stmTime = 10; % midline is 0: thus switches back and forth
  ostm.gapTime = 0;
  ostm.adaptationTime = 10;
  
  ostm.signalTime = 0;
  ostm.testingProb = 0;
  ostm.stmBkgType = 'none';;
  ostm.borderWidth = 0.05;

  opts.stimulus = ostm;
  
  ft.setOpts(opts);

  %ft.setDisplay(0);  
  ft.setDisplay('tracks',false);  
  ft.track();
  ft.save();

  

end


