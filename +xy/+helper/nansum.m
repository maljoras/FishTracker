function y = nansum(x,dim)
% like sum but omitting nans

x(isnan(x)) = 0;
y = sum(x,dim);