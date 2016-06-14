function res = getSwitchBasedTrackingResults(self,delinvif,forceif)
% RES = GETSWITCHBASEDTRACKINGRESULTS( ) returns the current switch
% based tracks.  GETSWITCHBASEDTRACKINGRESULTS(DELINVIF) sets times in
% RES.POS where a track was lost to NaN.
% GETSWITCHBASEDTRACKINGRESULTS(..,FORCIF) forces a regeneration of the pos/res
% structure.

  if nargin<3
    forceif = [];
  end
  if nargin<2
    delinvif = [];
  end
    
  res = self.getTrackingResults(delinvif,forceif,0);

end
