function res = getDagTrackingResults(self,timeRange,delinvif)
% RES = GETDAGTRACKINGRESULTS(TIMERANGE) returns the current DAG
% results in timerange. GETDAGTRACKINGRESULTS(...,DELINVIF) sets times in
% RES.POS where a track was lost to NaN.

  if nargin<3
    delinvif = [];
  end
  if nargin<2
    timeRange = [];
  end

  res = self.getTrackingResults(timeRange,1);
  if delinvif
    res.pos = self.deleteInvisible(res,'pos');
  end
  
end
