function sub = i2s(sz,ind)
% works like ind2sub but with v = s2i(sz,ind) instead of
% [v1,v2,v3,...] = ind2sub(sz,ind)
% works along dimension 2 (i.e. rows specify the index)
% CAUTION no error checking (if index out of range);
  
  if size(ind,2)~=1
    ind= ind';
  end
  
  sub=[];
  for i=1:length(sz)-1
    sub(:,end+1) = mod(ind,sz(i));
    sub(sub(:,end)==0,end) = sz(i);
    ind = (ind - sub(:,end))/sz(i) + 1;
  end
  sub(:,end+1) = ind;
