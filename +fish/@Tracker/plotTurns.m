function plotTurns(self,fishIds,plotTimeRange)
%   PLOTCENTERLINE(SELF,FISHIDS,PLOTTIMERANGE) plots the traces
%   with center line information. 
  
  if nargin<3 || isempty(plotTimeRange)
    plotTimeRange = self.timerange;
  end
  
  if nargin<2 || isempty(fishIds)
    fishIds = 1:self.nfish;
  end

  if  ~isfield(res.tracks,'centerLine') || isempty(res.tracks.centerLine)
    error('No centerLine data');
  end


  %res = self.getTrackingResults();
  %plotidx = res.t>=plotTimeRange(1) & res.t<plotTimeRange(2);


  b = getBendingLaplace(self,fishIds,plotTimeRange);
  
  mx = nanmean(b.clx(:,:,:),1);
  my = nanmean(b.cly(:,:,:),1);


  
  
  keyboard