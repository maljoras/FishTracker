function out = interpolateInvisible(self,field,kwidth,res)
% X = INTERPOLATEINVISIBLE(SELF,FIELD/MAT) returns the requested field
% in res.tracks with interpolated the consecutive invisible
% frames. INTERPOLATEINVISIBLE(..,kwidth) smothes the data along
% the frame dimension with running average  kernel.
  
  if nargin<4
    res = getTrackingResults(0,0);
  end
  
  msk = self.getInvisibleMsk(res);


  
  handled = 0;
  if ischar(field)
    if strcmp(field,'pos')
      out = res.pos;
      out(find(permute(cat(3,msk,msk),[1,3,2]))) = NaN;
      out = reshape(out,size(res.pos));
      handled = 1;
    elseif isfield(res.tracks,field)
      % always nFrames x nFish x nOther
      out = double(res.tracks.(field));

    else
      error(sprintf('Field %s does not exist in res.tracks.',field));
    end
  else
    if size(msk,1)~=size(field,1)
      error('Dimension mismatch');
    end
    out = field;
  end
  
  if ~handled
    sz = size(out);
    out = reshape(out,numel(msk),[]);
    out(msk,:) = NaN;
    out = reshape(out,sz);
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
  if nargin>2 && ~isempty(kwidth) && kwidth>0
    sz = size(out);
    out = convn(out,ones(kwidth,1)/kwidth,'same');
  end
  
end
