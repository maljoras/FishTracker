function out = interpolateInvisible(self,res,field,kwidth)
% X = INTERPOLATEINVISIBLE(SELF,RES,FIELD/MAT) returns the requested field
% in res.tracks with interpolated the consecutive invisible
% frames. INTERPOLATEINVISIBLE(..,KWIDTH) smoothes the data along
% the frame dimension with running average  kernel.
  

  out = self.deleteInvisible(res,field);  
  
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
  if nargin>3 && ~isempty(kwidth) && kwidth>0
    sz = size(out);
    out = xy.helper.movavg(out,kwidth);
  end
  
end
