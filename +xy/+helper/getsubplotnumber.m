function [r1,r2] = getsubplotnumber(n)
% [R1,R2] = GETSUBPLOTNUMBER(N) calculate appropiate subplot devisions based
% on N subplots


  r = ceil(sqrt(n));
  r2 =r;
  r1 = r;
  while r1*r2 - n >= r2
    r1=r1-1;
  end
