function out = getTurningPoints(self,fishIds, plotTimeRange,res)

  if nargin<3 || isempty(plotTimeRange)
    plotTimeRange = [-inf,inf];
  end
  
  if nargin<2 || isempty(fishIds)
    fishIds = 1:self.nfish;
  end
  
  if nargin<4 || isempty(res)
    res = self.getTrackingResults();
  end
  
  t = res.t;
  plotidx = t>=plotTimeRange(1) & t<plotTimeRange(2);

  if  ~isfield(res.tracks,'centerLine') || isempty(res.tracks.centerLine)
    error('No centerLine data');
  end

  lap = getBendingLaplace(self,fishIds,[],res);
  
  mclx = self.interpolateInvisible(shiftdim(mean(lap.clx,1),1),[],res);
  mcly = self.interpolateInvisible(shiftdim(mean(lap.cly,1),1),[],res);




  n = 10;
  fidx = 1;
  idx = 600:1000;
  w = zeros(length(idx),2);
  u = zeros(length(idx),2);
  for i = 1:length(idx)-n+1
    x = lap.clx(:,idx(i):idx(i+n-1),fidx);
    y = lap.cly(:,idx(i):idx(i+n-1),fidx);
    x = x(:);
    y = y(:);
    ind = find(isnan(x)|isnan(y));
    
    if ~isempty(ind)
      x(ind) = [];
      y(ind) = [];
    end
    
    
    w(i,:) = [y,ones(size(y))]\x;
    u(i,:) = [x,ones(size(x))]\y;
  
  end  
  x = lap.clx(:,idx,fidx);
  y = lap.cly(:,idx,fidx);

  [w1,resx,u1,resy] = fish.helper.contLinearRegression(y,x,n*7);

  
  keyboard

   
