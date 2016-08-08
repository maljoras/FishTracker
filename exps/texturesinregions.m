
LOAD = 1
COMPUTE =1

PLOT =0
TEST = 0;

path = '/data/videos/onlinelearning/new';
VIDID = 2;

if LOAD && ~exist('xyT','var')


  opts = [];
  opts.nbody = 3;
  
  videoFile = [path filesep mfilename sprintf('F%d-%d.avi',opts.nbody,VIDID)];

  opts.detector.inverted = 1;

  opts.stmif = 1;
  opts.display.videoHandler = true;
  opts.stimulus.screen = 1;
  opts.stimulus.presenter = 'xy.stimulus.PresenterTrackTextureRegions';
  opts.detector.adjustThresScale = 1.05;
  
  opts.bodywidth = 35;
  opts.bodylength = 150;
  opts.detector.history = 5000;
  
  xyT = xy.Tracker({0,videoFile},opts);
end


if COMPUTE

  if ~exist('sbbox','var')
    %sbbox = xy.Tracker.calibrateStimulusScreen(0,1);
    %sbbox =  [115,90,1550,1214];
    sbbox = [147,71,1555,1237];
  end
  opts = [];
  opts.avgVelocity = 5;
  opts.maxVelocity = 30;
    

  opts.display.displayEveryNFrame = 100;
  
  ostm = [];
  ostm.screenBoundingBox = sbbox;  
  ostm.usePredIdentityId = false;

  ostm.stmCol= parula(xyT.nbody);
  ostm.stmOnInt= 2; % in sec
  ostm.stmOnIntCV= 0.1;   

  ostm.stmOffInt= 4; % in sec
  ostm.stmOffIntCV= 0.5; 
  
  ostm.stmSize = xyT.bodylength;  
  ostm.regSizeFactorScale = [0,1,2];
  ostm.stmSizeFactor = 1;  
  ostm.xRegions = [1/3,2/3];

  ostm.stmShift = xyT.bodylength/4; % in px of xy.tracker frame
  ostm.stmShiftOri = 0;

  ostm.stmShiftOriSTD = -1; % random;
  ostm.stmShiftCV = 2;  
  ostm.stmSizeFactorCV = 0.25;
  
  ostm.stmVelThres = [2];
  ostm.stmBorderThres = 0.05;
  
  ostm.funLRStm = @(identityIds) identityIds;

  ostm.colBackground = [0,0,0];
  ostm.colBorder = [1,1,1];
  ostm.borderWidth = 0.00;


  ostm.adaptationTime = 30;
  ostm.tmax = 3600*2;
  
  
  
  opts.stimulus = ostm;
  xyT.setOpts(opts);

  %xyT.setDisplay(0);  
  xyT.setDisplay('tracks',false);  

  if TEST
    xy.helper.stimulusSimulator(xyT.stimulusPresenter);
  else
    xyT.track();
    xyT.save();
  end
  
  clear xyT;
  

end




if PLOT
  dagresults = 1;
  xyT.setDefaultResultType(dagresults);
  res = xyT.getTrackingResults();
  pos = xyT.interpolateInvisible(res,'pos');
  
  info = res.tracks.stmInfo;
  stmmsk = info(:,1,1)==xyT.stimulusPresenter.ID_STIMULUS;


  som = nanmean(pos,3);
  ssom = nanstd(sqrt(sum(bsxfun(@minus,pos,som).^2,2)),[],3);
  
  return
  
  pos = res.pos;
  velocity = xyT.deleteInvisible(res,'velocity');
  vel = sqrt(velocity(:,:,1).^2+velocity(:,:,2).^2);
  
  mpos1 = squeeze(nanmean(pos(stmmsk,:,:),1));
  mpos2 = squeeze(nanmean(pos(nonmsk,:,:),1));

  mvel1 = nanmean(vel(stmmsk,:),1);
  mvel2 = nanmean(vel(nonmsk,:),1);
  svel1 = xy.helper.stderr(vel(stmmsk,:),1);
  svel2 = xy.helper.stderr(vel(nonmsk,:),1);

  
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
  vel = xyT.deleteInvisible(res,'velocity');
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
