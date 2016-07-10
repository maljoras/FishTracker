
LOAD = 1
COMPUTE =1

PLOT =0
TEST = 0;

path = '/data/videos/onlinelearning';
VIDID = 

if LOAD && ~exist('ft','var')


  opts = [];
  opts.nfish = 3;
  
  videoFile = [path filesep mfilename sprintf('F%d-%d.avi',opts.nfish,VIDID)];

  opts.detector.inverted = 1;

  opts.stmif = 1;
  opts.display.videoHandler = true;
  opts.stimulus.screen = 1;
  opts.stimulus.presenter = 'fish.stimulus.PresenterTrackTextureRegions';
  opts.detector.adjustThresScale = 1.05;
  
  opts.fishwidth = 25;
  opts.fishlength =  90;
  
  
  ft = fish.Tracker({0,videoFile},opts);
end


if COMPUTE

  if ~exist('sbbox','var')
    %sbbox = fish.Tracker.calibrateStimulusScreen();
    %sbbox =  [115,90,1550,1214];
    sbbox = [114,69,1557,1226];
  end
  opts = [];
  opts.avgVelocity = 5;
  opts.maxVelocity = 30;
    

  opts.display.displayEveryNFrame = 100;
  
  ostm = [];
  ostm.screenBoundingBox = sbbox;  
  ostm.usePredFishId = false;

  ostm.stmCol= parula(ft.nfish);
  ostm.stmSwitchInt= 3; % in sec
  ostm.stmSwitchIntCV= 0.5; 
  
  ostm.stmSize = ft.fishlength;  
  ostm.regSizeFactorScale = [0,0.5,1.5];
  ostm.stmSizeFactor = 1;  
  ostm.xRegions = [1/3,2/3];

  ostm.stmShift = ft.fishlength/4; % in px of fish.tracker frame
  ostm.stmShiftOri = 0;

  ostm.stmShiftOriSTD = -1; % random;
  ostm.stmShiftCV = 2;  
  ostm.stmSizeFactorCV = 0.5;
  
  ostm.stmVelThres = [2];
  ostm.stmBorderThres = 0.05;
  
  ostm.funLRStm = @(fishIds) fishIds;

  ostm.colBackground = [0,0,0];
  ostm.colBorder = [1,1,1];
  ostm.borderWidth = 0.0;


  ostm.adaptationTime = 30;
  ostm.tmax = 3600*2;
  
  
  
  opts.stimulus = ostm;
  ft.setOpts(opts);

  %ft.setDisplay(0);  
  ft.setDisplay('tracks',false);  

  if TEST
    fish.helper.stimulusSimulator(ft.stimulusPresenter);
  else
    ft.track();
    ft.save();
  end
  
  clear ft;
  

end




if PLOT
  dagresults = 1;
  ft.setDefaultResultType(dagresults);
  res = ft.getTrackingResults();
  pos = ft.interpolateInvisible(res,'pos');
  
  info = res.tracks.stmInfo;
  stmmsk = info(:,1,1)==ft.stimulusPresenter.ID_STIMULUS;


  som = nanmean(pos,3);
  ssom = nanstd(sqrt(sum(bsxfun(@minus,pos,som).^2,2)),[],3);
  
  return
  
  pos = res.pos;
  velocity = ft.deleteInvisible(res,'velocity');
  vel = sqrt(velocity(:,:,1).^2+velocity(:,:,2).^2);
  
  mpos1 = squeeze(nanmean(pos(stmmsk,:,:),1));
  mpos2 = squeeze(nanmean(pos(nonmsk,:,:),1));

  mvel1 = nanmean(vel(stmmsk,:),1);
  mvel2 = nanmean(vel(nonmsk,:),1);
  svel1 = fish.helper.stderr(vel(stmmsk,:),1);
  svel2 = fish.helper.stderr(vel(nonmsk,:),1);

  
  figure;
  subplot(2,1,1);
  plot(res.tabs,squeeze(pos(:,1,:)));
  title('x')
  hold on;
  plot(res.tabs,stmmsk*100-100,'linewidth',2);

  subplot(2,1,2);
  plot(res.tabs,squeeze(pos(:,2,:)));
  title('y')
  hold on;
  plot(res.tabs,stmmsk*100-100,'linewidth',2);

  
  figure;
  a = subplot(2,1,1);
  vel = ft.deleteInvisible(res,'velocity');
  vabs = sqrt(vel(:,:,1).^2 + vel(:,:,2).^2);
  vabs(vabs>100) = NaN;
  
  acc = [diff(vabs);zeros(1,size(vabs,2))];
  %maxacc = 20;
  %vabs(abs(acc)>maxacc) = NaN;
  
  plot(res.tabs,vabs);
  ylabel('Velocity [px/s]');
  xlabel('Time [s]');
    hold on;
  plot(res.tabs,stmmsk*10-10,'linewidth',2);

  
  b = subplot(2,1,2);
  plot(res.tabs(1:end-1),conv2(diff(vabs),ones(5,1)/5,'same'));
  ylabel('Accelleration [px/s^2]');
  xlabel('Time [s]');
  hold on;
  plot(res.tabs,stmmsk*1-1,'linewidth',2);

  linkaxes([a,b],'x')
  
end
