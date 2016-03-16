function msk = getInvisibleMsk(self)
% MSK = GETINVISIBLEMSK(SELF) gets the msk to delete consecutive
% invisible frames from the results
  
  c = self.res.tracks.consecutiveInvisibleCount;
  msk = c>0;
  %also delete the last visible (often noise detection)
  dmsk = diff(msk)==1;
  msk(1:end-1,:) = msk(1:end-1,:) | dmsk;
end
