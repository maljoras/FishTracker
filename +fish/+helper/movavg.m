function xm = movavg(x,nconv,dim)
% XM = MOVAVG(X,NCONV,DIM) calculates a moving average box-car filter
% of length NCONV on the DIM of X. Ignores NAN.

  if nargin<3
    dim = 1;
  end
  
  o = ones(1,max(length(size(x)),dim));
  o(dim) = nconv;
  k = ones(o);

  msk = isnan(x);
  n = convn(double(~msk),k,'same');
  x(msk) = 0;
  
  xm = convn(x,k,'same');
  
  idx = find(~n);
  xm(idx) = 0;
  n(idx) = 1;
  xm = xm./n;
  xm(idx) = NaN;