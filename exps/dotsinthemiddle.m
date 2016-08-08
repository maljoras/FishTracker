
LOAD = 0
COMPUTE =0
PLOT =1;


if LOAD && ~exist('xyT','var')
  videoFile = '/data/videos/onlinelearning/dotsInTheMiddle2_4fish.avi';
  opts = [];
  opts.detector.inverted = 1;
  
  opts.nbody = 4;

  opts.stmif = 1;
  opts.display.videoHandler = true;
  opts.stimulus.screen = 1;
  opts.stimulus.presenter = 'xy.stimulus.PresenterOnlineLearningCue';
  opts.detector.adjustThresScale = 1.05;
  
  opts.bodywidth = 25;
  opts.bodylength =  90;
  
  
  xyT = xy.Tracker({0,videoFile},opts);
end


if COMPUTE

  if ~exist('sbbox','var')
    %sbbox = xyT.calibrateStimulusScreen();
    sbbox = [88,63,1584,1247];
  end
  opts = [];
  opts.avgVelocity = 5;
  opts.maxVelocity = 20;
    

  opts.display.displayEveryNFrame = 100;
  
  ostm = [];
  ostm.screenBoundingBox = sbbox;  
  ostm.usePredIdentityId = false;

  ostm.stmCol= parula(xyT.nbody);
  ostm.stmLambda= 0.0;
  
  ostm.stmSize = ones(xyT.nbody,1)*xyT.bodylength.*(1+rand(xyT.nbody,1));  

  ostm.stmType = 'bodytextures';
  ostm.stmSizeFactor = 1.2;
  ostm.stmShift = xyT.bodylength; % in px of xy.tracker frame
  ostm.stmShiftOri = 0;
  ostm.stmVelThres = 2;
  ostm.stmBorderThres = 0.075;
  
  ostm.funLRStm = @(identityIds) identityIds;
  ostm.lrif = 1;
  ostm.lrSwitchProb = 1;
  ostm.midline = 0;
  ostm.colBackground = [0,0,0];
  
  ostm.colBorder = [1,1,1];
  ostm.borderWidth = 0.075;

  ostm.nRound = 250;

  ostm.stmTime = 15; % midline is 0: thus switches back and forth
  ostm.gapTime = 0;
  ostm.adaptationTime = 30;
  
  ostm.signalTime = 0;
  ostm.testingProb = 0;
  ostm.stmBkgType = 'simpleborder';
  ostm.borderWidth = 0.05;

  opts.stimulus = ostm;
  
  xyT.setOpts(opts);

  %xyT.setDisplay(0);  
  xyT.setDisplay('tracks',false);  
  xyT.track();
  xyT.save();

  

end




if PLOT
  dagresults = 1;
  xyT.setDefaultResultType(dagresults);
  res = xyT.getTrackingResults();
  
  info = res.tracks.stmInfo;
  stmmsk = info(:,1,1)==xyT.stimulusPresenter.ID_STIMULUS_LR;
  nonmsk = info(:,1,1)==xyT.stimulusPresenter.ID_STIMULUS_RL;


  pos = xyT.deleteInvisible(res,'pos');;
  velocity = xyT.deleteInvisible(res,'velocity');
  vel = sqrt(res.tracks.velocity(:,:,1).^2+res.tracks.velocity(:,:,2).^2);
  
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
  vabs = sqrt(res.tracks.velocity(:,:,1).^2 + res.tracks.velocity(:,:,2).^2);
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
