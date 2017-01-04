function y = nanmean(x,dim)
% like mean but omitting nans

nans = isnan(x);
x(nans) = 0;
non = sum(~nans,dim);
non(non==0)=nan; 
y = sum(x,dim)./non;