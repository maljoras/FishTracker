function [combinedFT varargout] = trackParallel(self,inTimeRange, ...
                                                tOverlap,minDurationPerWorker)
  %  FTNEW = FT.TRACKPARALLEL() will track the results in parallel and creates a
  %  new FT object with tracks combined
  % 
  %  FT.TRACKPARALLEL(..,TOVERLAP,MINDURATIONPERWORKER) sepecifies the overlap
  %  between workers in seconds and the minmal video time per worker in
  %  seconds
  % 
  % CAUTION: if not calculated in parallel (w.g. for too short data) the
  % returned handle object will be the IDENTICAL handle (only a reference to
  % the same data).  Only in case of parallel processing the returned object
  % will have a NEW handle (and thus reference new data).

  if ~exist('inTimeRange','var')
    inTimeRange = [];
  end
  self.videoHandler.timeRange = inTimeRange;
  self.timerange = self.videoHandler.timeRange;
  self.videoHandler.reset();

  dt = 1/self.videoHandler.frameRate;
  if ~exist('minDurationPerWorker','var')
    minDurationPerWorker = 300; % 5 minutes;
  end

  if ~exist('tOverlap','var') || isempty(tOverlap)
    % allow enough time for the classifer to init in the overlapping period
    tOverlap = 5*dt*max([self.nFramesForInit,self.nFramesForUniqueUpdate]); 
    tOverlap = max(tOverlap,self.opts.tracks.costtau*self.avgTimeScale*dt);
  end

  assert(minDurationPerWorker>tOverlap);
  totalDuration = diff(self.timerange);

  if totalDuration<2*tOverlap || totalDuration < 2*(minDurationPerWorker-tOverlap)
    % just on a single computer;
    self.track();
    combinedFT = self; % SHOULD IMPLEMENT COPY OBJECT HERE !!!
    varargout{1} = [];
    return
  end
  pool = parpool(2);

  maxDuration = min(diff(self.timerange)/pool.NumWorkers,3600);
  nparts = floor((diff(self.timerange)/maxDuration));

  sec = linspace(self.timerange(1),self.timerange(2),nparts+1);
  secEnd = sec(2:end);
  secStart = max(sec(1:end-1)-tOverlap,sec(1));
  timeRanges = [secStart(:),secEnd(:)];

  nRanges = size(timeRanges,1);

  % copy objects
  xy.helper.verbose('Parallel tracking from %1.1f to %1.1fh on %d workers',...
          self.timerange(1)/3600,self.timerange(2)/3600,pool.NumWorkers);
  xy.helper.verbose('%d processes with %1.1fm duration and %1.1fsec overlap.',...
          nRanges,max(diff(timeRanges,[],2))/60,tOverlap);

  % CAUTION THIS YIELDS A REFERENCE COPY ONLY !!! HOWEVER, RUNNING THIS WITH PARFOR LOOP WILL CHANGE
  % EACH HANDLE TO A VALUE BASED OBJECT ON EACH WORKER SO WE DO NOT NEED A FORMAL COPY METHOD. SEEMS TO
  % BE A HACK IN MATLAB.
  FTs = repmat(self,[1,nRanges]);
  res = {};
  parfor i = 1:nRanges
    ft = FTs(i);
    ft.setupSystemObjects(ft.videoFile); % redo the videoHandler;
    ft.displayif = 0;
    ft.track(timeRanges(i,:));
    res{i} = ft;
  end
  combinedFT = combine(res{:});
  if nargout>1
    varargout{1} = res;
  end
end


