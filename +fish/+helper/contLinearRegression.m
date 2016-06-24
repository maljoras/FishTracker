function [w,resx,u,resy] =contLinearRegression(y,x,nconv)
% [W,U] =CONTLINEARREGRESSION(X,Y) continous linear regression
% y=w_1*x + w_2 and x=u_1*y + u_2 with local kernel size N

 

  x = x(:);
  y = y(:);

  msk = ~isnan(x) & ~isnan(y);
  x(~msk) = 0;
  y(~msk) = 0;
  
  k = ones(nconv,1);
  n = conv(double(msk),k,'same');
  sx = conv(x,k,'same');
  sy = conv(y,k,'same');
  sigxy = conv(x.*y,k,'same');
  sigxx = conv(x.^2,k,'same');
  sigyy = conv(y.^2,k,'same');
  n(~n) = NaN;

  Dx = sigxx.*n - sx.^2;
  w = [(n.*sigxy - sx.*sy)./Dx,(-sx.*sigxy + sigxx.*sy)./Dx];

  % residuel std
  resx = sqrt(-2*(w(:,1).*sigxy+w(:,2).*sy) + sigyy + w(:,1).^2.*sigxx + ...
         2*w(:,1).*w(:,2).*sx + w(:,2).^2.*n);
  
  if nargout>2

    Dy = sigyy.*n - sy.^2;
    u = [(n.*sigxy - sx.*sy)./Dy,(-sy.*sigxy + sigyy.*sx)./Dy];
  
    resy = sqrt(-2*(u(:,1).*sigxy+u(:,2).*sx) + sigxx + u(:,1).^2.*sigyy + ...
         2*u(:,1).*u(:,2).*sy + u(:,2).^2.*n);
  
  end
  
  
