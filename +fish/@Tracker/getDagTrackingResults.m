function res = getDagTrackingResults(self,delinvif,forceif)
% RES = GETDAGTRACKINGRESULTS( ) returns the current DAG results.  
% GETDAGTRACKINGRESULTS(DELINVIF) sets times in RES.POS where a track was
% lost to NaN.
% GETTRACKINGRESULTS(..,FORCIF) forces a regeneration of the
% pos/res structure. 

  if nargin<3
    forceif = [];
  end
  if nargin<2
    delinvif = [];
  end
    
  res = self.getTrackingResults(delinvif,forceif,1);

end
