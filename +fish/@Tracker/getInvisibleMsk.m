function msk = getInvisibleMsk(self,res)
% MSK = GETINVISIBLEMSK(SELF) gets the msk to delete consecutive
% invisible frames from the results
  
  if nargin<2
    c = self.res.tracks.consecutiveInvisibleCount;
  else
    c = res.tracks.consecutiveInvisibleCount;
  end
  
  msk = c>0;
  %also delete the last visible (often noise detection)
  dmsk = diff(msk)==1;
  msk(1:end-1,:) = msk(1:end-1,:) | dmsk;
end
