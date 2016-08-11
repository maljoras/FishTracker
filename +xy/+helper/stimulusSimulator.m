function stimulusSimulator(stmObj,varargin)
% STIMULUSSIMULATOR(STMOBJ,...) simulatas the stimulations objects

  def.stmObj = 'xy.stimulus.PresenterFlash';
  
  def.opts.dt =  1/30;
  def.opts.tmax = 1000;
  def.opts.screen = 0;
  def.opts.windowSize = [600,400]; % PsychToolb : x,y
  def.opts.windowOrigin = [0,0];
  def.opts.nindiv = 3;
  def.opts.bodySize = 20;
  def.opts.velocity = 1;
  def.opts.frameSize = [800,1000]; % Matlab: y,x
  def.opts.timeFactor = 1;
  def.opts.stimulus = []; % props ;
  doc.stimulus = 'structure for additional options of the stim obj';

  xy.helper.parseInputs;
  if HELP, return;end;
  
  % generate traces

  if ischar(stmObj)
    stmObj = eval(stmObj);
  end

  w =[opts.windowOrigin, opts.windowOrigin + opts.windowSize];
  if ~isempty(opts.stimulus)
    for f = {'screen','screenBoundingBox'}
      if isfield(opts.stimulus,f{1})
        opts.stimulus = rmfield(opts.stimulus,f{1});
      end
    end
    stmObj.setOpts(opts.stimulus);
  end

  stmObj.screenBoundingBox = [];
  stmObj.screen = opts.screen;

  stmObj.init(w);
  stmObj.muteAllFlipping = true;
  if isinf(opts.tmax)
    stmObj.progressBar = 0;  
  end
  stmObj.tmax = opts.tmax/opts.timeFactor;
  
  %stmObj.stmTime;
  stmObj.reset();

  % make fake tracks
  nt = floor(opts.tmax/opts.dt/opts.timeFactor);
  trace = subMakeTrace(nt,opts.nindiv,opts.velocity*opts.dt*opts.timeFactor);


  cmap = jet(opts.nindiv);
  localTimeReference = tic;
  fl = opts.bodySize;
  for iframe = 1:nt

    % make fake tracks
    tdelay = tic;
    tracks = [];
    for i = 1:opts.nindiv
      tracks(i).centroid = trace(iframe,:,i).*opts.frameSize([2,1]);
      tracks(i).location = tracks(i).centroid;
      tracks(i).bbox = [tracks(i).centroid-[fl,fl]/2,fl,fl];
      tracks(i).identityId = i;
      tracks(i).predIdentityId = i;
      tracks(i).velocity = (tracks(i).centroid-trace(max(iframe-1,1),:,i))./opts.dt/opts.timeFactor;
      tracks(i).segment.Image = ones(fl,fl);
      tracks(i).segment.Orientation = -atan2(tracks(i).velocity(1), ...
                                            tracks(i).velocity(2))/pi*180 +90;
    
    end
    pause(opts.dt-toc(tdelay));
    localTime = toc(localTimeReference);        

    % call stimulus object
    lt = localTime*opts.timeFactor;
    tracks = stmObj.step(tracks,opts.frameSize,lt);
  
    % plot identityys
    for i = 1:opts.nindiv
      stmObj.plotDot(trace(iframe,1,i),trace(iframe,2,i),opts.bodySize,cmap(i,:))
    end
    stmObj.flip(true);    
    
    fprintf('StmInfo Track 1: %s\r',num2str(tracks(1).stmInfo));
  
  
    if stmObj.isFinished(lt)
      break
    end
  end

  stmObj.muteAllFlipping = false;
  


end


function trace = subMakeTrace(nt,nindiv,cutoff);

  [Hb,Ha] = butter(4,cutoff,'low');
  rpath = filtfilt(Hb,Ha,exp(1i*4*pi*(rand(nt,2,nindiv)-0.5)));
  trace = angle(rpath)/pi/2 + 0.5;
  
  
end

  



