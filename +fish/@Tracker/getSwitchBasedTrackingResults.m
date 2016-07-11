function res = getSwitchBasedTrackingResults(self,timeRange, delinvif)
% RES = GETSWITCHBASEDTRACKINGRESULTS(TIMERANGE) returns the current
% switch based tracks.  GETSWITCHBASEDTRACKINGRESULTS(..,DELINVIF)
% sets times in RES.POS where a track was lost to NaN.

  if nargin<3
    delinvif = [];
  end
  if nargin<2
    timeRange = [];
  end
    
  res = self.getTrackingResults(timeRange,0);
  if delinvif
    res.pos = self.deleteInvisible(res,'pos');
  end
end
