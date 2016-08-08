function se = stderr(z,dim);
% computes std err (nan neglected)
  
  if ~exist('dim','var')
    dim = 1;
    z = shiftdim(z);
  end
  
  if size(z,dim)>1
    tmp = sqrt(sum(~isnan(z),dim));
    izero = find(~tmp);
    tmp(izero) = 1;
    se = nanstd(z,[],dim)./tmp;
    se(izero) = NaN;
  else
    se= zeros(size(z));
    se(:) = NaN;
  end
  