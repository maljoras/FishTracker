function index = s2i(sz,v);
% works like sub2ind but with s2i(sz,[v1,v2,v3,...]) instead of
% sub2ind(sz,v1,v2,v3,...)
% works along dimension 2 (i.e. rows specify the indices)
% CAUTION no error checking (if index out of range);
  
  index = v(:,1);
  sc=1;
  for i=2:length(sz)
    sc = sc*(sz(i-1));
    index = index + sc*(v(:,i)-1);
  end
  
   