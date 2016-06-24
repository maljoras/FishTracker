function msk = getInvisibleMsk(self,res)
% MSK = GETINVISIBLEMSK(SELF,RES) gets the msk to delete consecutive
% invisible frames from the results
  
  if nargin<2
    res = getTrackingResults(0,0);
  end
  
  c = res.tracks.consecutiveInvisibleCount;
  msk = c>0;
  %also delete the last visible (often noise detection)
  dmsk = diff(msk)==1;
  msk(1:end-1,:) = msk(1:end-1,:) | dmsk;

end
