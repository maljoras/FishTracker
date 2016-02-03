function out = deleteInvisible(self,field)
% X = DELETEINVISIBLE(SELF,filed) returns the requested field
% in res.tracks with deleted (NaN) the consecutive invisible
% frames 
  
  if isempty(self.res)
    return;
  end
  msk = self.getInvisibleMsk();
  
  if strcmp(field,'pos')
    out = self.res.pos;
    out(permute(cat(3,msk,msk),[1,3,2])) = NaN;
  elseif isfield(self.res.tracks,field)
    % always nFrames x nFish x nOther
    out = double(self.res.tracks.(field));
    sz = size(out);
    out = reshape(out,numel(msk),[]);
    out(msk,:) = NaN;
    out = reshape(out,sz);
  else
    error(sprintf('Field %s does not exist in res.tracks.',field));
  end
  
end
