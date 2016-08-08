
LOAD = 1
COMPUTE =0

if LOAD || ~exist('ft','var')
  videoFile = ''; %'/data/videos/onlinelearning/test.avi';
  opts = [];
  opts.detector.inverted = 1;
  opts.nanimals = 3;
  opts.stmif = 1;

  opts.stimulus.screen = 1;
  opts.stimulus.screenBoundingBox = [150,53,1580,1256];
  opts.stimulus.presenter = 'xy.stimulus.PresenterPlaneIndiv';

  
  opts.stimulus.stmBkgType = 'none'; %plane
  opts.stimulus.sideMarkerCol = [1,1,0;1,0,1];
  opts.stimulus.testingProb = 0.1;
  opts.stimulus.stmCol = [1,1,1]; % stimulus color (RGB [1,1,1] for white)    
  opts.stimulus.colBackgound = [0,0,0]; % background color (RGB [0,0,0] for black)

  
  opts.stimulus.stmTime = 5; 
  opts.stimulus.adaptationTime = 1;
  opts.stimulus.stmSize = 100;

  opts.stimulus.nRoundsPerGroup = 4; % switch betweem left & right stimulation
  opts.stimulus.nPauseStmGroupsPerExp = 2;  %   
  opts.stimulus.nanimalsStimPerExps = [1,5,3,4,5];

  
  %opts.fishwidth = 30;
  %opts.fishlength = 150;


  ft = xy.Tracker({0,videoFile},opts);
end


if COMPUTE
  %xy.helper.stimulusSimulator('xy.stimulus.PresenterPlaneIndiv','nanimals',5,'stimulus',opts.stimulus,'timeFactor',3)

  
  
  if ~exist('sbbox','var')
    %sbbox = ft.calibrateStimulusScreen();
    sbbox = [150,53,1580,1256];
  end
  ft.setOpts('stimulus.screenBoundingBox',sbbox,'display.displayEveryNFrame',100);
  ft.setDisplay(0);  

  
  ft.track();
  ft.save();

  

end


