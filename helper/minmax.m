function mm = minmax(X)

  mm = cat(2,min(X,[],2),max(X,[],2));