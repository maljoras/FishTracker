function out = trackFishInParallel(timerange,duration,varargin)


  args = varargin;

  % maybe save ?
  if duration>diff(timerange)
    % just on a single computer;
    res = trackPart(timerange,args{:});
    return
  end
  

  overlap = 2; %sec
  
  if duration<2*overlap
    error('duration of segment too short');
  end
  
  
  pool = gcp;
  sec = timerange(1):duration:timerange(2);
  res = {};
  parfor i = 1:length(sec)-1
    res{i} =  trackPart([max(sec(i)-overlap,timerange(1)), ...
                        min(sec(i+1)+overlap,timerange(2))],args);
  end
  
  for i = 1:length(res)
    out{i} = res{i};
  end

  
  
  
function res = trackPart(tr,args)


  ft = FishTracker(args{:},'displayif',0,'timerange',tr);

  ft.track();
  res = ft;
  %res = ft.getTrackingResults;

