function p = getXYTRoot()
% P = GETXYTROOT() gets the root path of fis.Tracker (the path with
% the +xy directory)
  p = which('xy.helper.getXYTRoot');
  idx = strfind(p,'+xy');
  p = p(1:idx(1)-1);
  