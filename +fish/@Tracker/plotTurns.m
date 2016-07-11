function plotTurns(self,fishIds,plotTimeRange)
%   PLOTCENTERLINE(SELF,FISHIDS,PLOTTIMERANGE) plots the traces
%   with center line information. 
  
  if nargin<3 || isempty(plotTimeRange)
    plotTimeRange = self.timerange;
  end
  
  if nargin<2 || isempty(fishIds)
    fishIds = 1:self.nfish;
  end

  
  res = self.getTrackingResults(plotTimeRange);
  self.getTurningPoints(res,1);
  