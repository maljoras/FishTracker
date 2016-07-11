function res = getTrackingResults_(self,dagif)
% RES = GETTRACKINGRESULTS_(DAGIF) returns the current
% results. For internal functions only (no verbose)
  
  if nargin<2
    dagif = self.getDefaultResultsType();
  end
  
  if isempty(self.res)
    error('No results available.');
  else
    if dagif 
      res = self.res.dag;
    else
      res = self.res.swb;
    end
  end
end

