function stimulusSimulator(stmObj,varargin)
% STIMULUSSIMULATOR(STMOBJ,...) simulatas the stimulations objects

  def.stmObj = 'fish.stimulus.PresenterFlash';
  
  def.opts.dt =  1/30;
  def.opts.tmax = 1000;
  def.opts.screen = 0;
  def.opts.windowSize = [600,400]; % PsychToolb : x,y
  def.opts.windowOrigin = [0,0];
  def.opts.nfish = 3;
  def.opts.fishSize = 20;
  def.opts.velocity = 1;
  def.opts.frameSize = [800,1000]; % Matlab: y,x
  def.opts.timeFactor = 1;
  def.opts.stimulus = []; % props ;
  doc.stimulus = 'structure for additional options of the stim obj';

  fish.helper.parseInputs;
  if HELP, return;end;
  
  % generate traces

  if ischar(stmObj)
    stmObj = eval(stmObj);
  end
  
  stmObj.screen = opts.screen;
  
  w =[opts.windowOrigin, opts.windowOrigin + opts.windowSize];

  if ~isempty(opts.stimulus)
    for f = {'screen','screenBoundingBox'}
      if isfield(opts.stimulus,f{1})
        opts.stimulus = rmfield(opts.stimulus,f{1});
      end
    end
    stmObj.setOpts(opts.stimulus);
  end

  stmObj.init(w);
  stmObj.muteAllFlipping = true;
  %stmObj.stmTime;
  
  % make fake tracks
  nt = floor(opts.tmax/opts.dt/opts.timeFactor);
  trace = subMakeTrace(nt,opts.nfish,opts.velocity*opts.dt*opts.timeFactor);


  cmap = jet(opts.nfish);
  localTimeReference = tic;
  for iframe = 1:nt

    % make fake tracks
    tdelay = tic;
    tracks = [];
    for i = 1:opts.nfish
      tracks(i).centroid = trace(iframe,:,i).*opts.frameSize([2,1]);
      tracks(i).fishId = i;
      tracks(i).predFishId = i;
    end
    pause(opts.dt-toc(tdelay));
    localTime = toc(localTimeReference);        

    % call stimulus object
    tracks = stmObj.step(tracks,opts.frameSize,localTime*opts.timeFactor);
  
    % plot fishys
    for i = 1:opts.nfish
      stmObj.plotDot(trace(iframe,1,i),trace(iframe,2,i),opts.fishSize,cmap(i,:))
    end
    stmObj.flip(true);    
    
    fprintf('StmInfo Track 1: %s\r',num2str(tracks(1).stmInfo));
  end

  stmObj.muteAllFlipping = false;
  


end


function trace = subMakeTrace(nt,nfish,cutoff);

  [Hb,Ha] = butter(4,cutoff,'low');
  rpath = filtfilt(Hb,Ha,exp(1i*4*pi*(rand(nt,2,nfish)-0.5)));
  trace = angle(rpath)/pi/2 + 0.5;
  
  
end

  



