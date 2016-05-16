function p = getFTRoot()
% P = GETFTROOT() gets the root path of fis.Tracker (the path with
% the +fish directory)
  p = which('fish.helper.getFTRoot');
  idx = strfind(p,'+fish');
  p = p(1:idx(1)-1);
  