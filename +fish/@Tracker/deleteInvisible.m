function out = deleteInvisible(self,res,field)
% X = DELETEINVISIBLE(SELF,RES/MAT,FIELD/CIC) returns the requested field
% in res.tracks with deleted (NaN) the consecutive invisible
% frames 
  
  if ~ischar(field)
    out = res;
    msk = self.getInvisibleMsk(field);
  else
    out = self.getResField(res,field,0);
    msk = self.getInvisibleMsk(res);
  end
  
  if ischar(field) && strcmp(field,'pos')
    out(permute(cat(3,msk,msk),[1,3,2])) = NaN;
  elseif (size(out,1)==size(msk,1)) && (size(out,2)==size(msk,2))
    % always nFrames x nFish x nOther
    sz = size(out);
    out = reshape(out,numel(msk),[]);
    out(msk,:) = NaN;
    out = reshape(out,sz);
  else
    error('dimension mismatch');
  end
  
end
