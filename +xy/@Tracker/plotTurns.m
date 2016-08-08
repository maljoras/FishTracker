function plotTurns(self,identityIds,plotTimeRange)
%   PLOTCENTERLINE(SELF,IDENTITYIDS,PLOTTIMERANGE) plots the traces
%   with center line information. 
  
  if nargin<3 || isempty(plotTimeRange)
    plotTimeRange = self.timerange;
  end
  
  if nargin<2 || isempty(identityIds)
    identityIds = 1:self.nbody;
  end

  
  res = self.getTrackingResults(plotTimeRange);
  self.getTurningPoints(res,1);
  