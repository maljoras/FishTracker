function se = stderr(z,dim)
% computes std err (nan neglected)
  
  if ~exist('dim','var')
    dim = 1;
    z = shiftdim(z);
  end
  
  if size(z,dim)>1
    n = sum(~isnan(z),dim);
    izero = ~n;
    n(izero) = 1;
    % do not use nanstd...
    cz = bsxfun(@minus,z,xy.helper.nanmean(z,dim));
    s = sqrt(xy.helper.nansum(abs(cz).^2,dim)./max(n-1,1));
    se = s./sqrt(n);
    se(izero) = NaN;
  else
    se= zeros(size(z));
    se(:) = NaN;
  end
  