function out = interpolateInvisible(self,field,kwidth)
% X = INTERPOLATEINVISIBLE(SELF,FIELD) returns the requested field
% in res.tracks with interpolated the consecutive invisible
% frames. INTERPOLATEINVISIBLE(..,kwidth) smothes the data along
% the frame dimension with running average  kernel.
  
  if isempty(self.res)
    return;
  end
  msk = self.getInvisibleMsk();
  if strcmp(field,'pos')
    out = self.res.pos;
    out(find(permute(cat(3,msk,msk),[1,3,2]))) = NaN;
    out = reshape(out,size(self.res.pos));
  
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

  sz = size(out);
  for i = 1:prod(sz(2:end))
    x = out(:,i);
    msk1 = ~isnan(x);
    msk2 = ~msk1;
    if ~all(msk1)
      ix = interp1(find(msk1),x(msk1),find(msk2),'linear');
      x(msk2) = ix;
    end
    out(:,i) = x;
  end
  

  % maybe better change to gaussian kernel !
  if nargin>2 && kwidth>0
    sz = size(out);
    out = convn(out,ones(kwidth,1)/kwidth,'same');
  end
  
end
