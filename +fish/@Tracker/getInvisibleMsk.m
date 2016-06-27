function msk = getInvisibleMsk(self,res)
% MSK = GETINVISIBLEMSK(SELF,RES/CIC) gets the msk to delete consecutive
% invisible frames from the results
  
  if isstruct(res)
    c = res.tracks.consecutiveInvisibleCount;
  else
    c = res;
  end
  
  msk = c>0;
  %also delete the last visible (often noise detection)
  dmsk = diff(msk)==1;
  msk(1:end-1,:) = msk(1:end-1,:) | dmsk;

end
