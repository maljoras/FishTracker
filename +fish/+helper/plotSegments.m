function plotSegments(seg)

  
  if ~length(seg)
    return;
  end

  clf; 
  [r1,r2] = fish.helper.getsubplotnumber(length(seg));
  t = linspace(0,2*pi);
  for i = 1:length(seg)
    
    subsubplot(r1,r2,i,2,1,1);
    x = 1:size(seg(i).FilledImage,2);
    y = 1:size(seg(i).FilledImage,1);
    imagesc(x,y,seg(i).FilledImage);
    [X,Y] = meshgrid(x,y);
    
    % plot ellipse
    phi = seg(i).Orientation/180*pi;
    a = seg(i).Size(1)/2;
    b = seg(i).Size(2)/2;
    m = seg(i).Centroid - seg(i).BoundingBox(1:2);
    x = m(1) + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
    y = m(2) + a*cos(t)*sin(phi) + b*sin(t)*cos(phi); 

    
    hold on;
    plot(x,y,'w');
% $$$ 
% $$$     e = sqrt(max(a,b)^2-min(a,b)^2);
% $$$     X = X - m(1);
% $$$     Y = Y - m(2);
% $$$     Zr= [X(:),Y(:)]*R2d(phi);
% $$$     Xr = reshape(Zr(:,1),size(X));
% $$$     Yr = reshape(Zr(:,2),size(Y));
% $$$ 
% $$$ 
% $$$     
% $$$ % $$$     
% $$$ % $$$     rA = sqrt((Xr+e).^2 + Yr.^2);
% $$$ % $$$     rB = sqrt((Xr-e).^2 + Yr.^2);
% $$$ % $$$ 
% $$$ % $$$     mu = acosh((rA + rB)/2/e);
% $$$ % $$$     eta = sign(Yr).*acos((rA-rB)/2/e);
% $$$ % $$$     
% $$$ % $$$     mu1 = mu-0.2;
% $$$ % $$$     
% $$$     testimg = X.^2./b^2 + Y.^2/a^2;
% $$$     
% $$$     mu = linspace(0,1,size(X,1));
% $$$     eta = linspace(-pi,pi,size(X,2));
% $$$     [Mu,Eta] = ndgrid(mu,eta);
% $$$     
% $$$     X1 = e*cosh(Mu).*cos(Eta)  ;
% $$$     Y1 = e*sinh(Mu).*sin(Eta) ;
% $$$ 
% $$$ % $$$     Zr= [X1(:),Y1(:)]*R2d(phi);
% $$$ % $$$     X1 = reshape(Zr(:,1),size(X1));
% $$$ % $$$     Y1 = reshape(Zr(:,2),size(Y1));
% $$$ 
% $$$     X1 = X1 + m(1);
% $$$     Y1 = Y1 + m(2);
% $$$     
% $$$     img = cv.remap(testimg,X1,Y1);
% $$$     

    
    daspect([1,1,1])
    axis off;
    hold on; 
    c = seg(i).CenterLine ;
    if ~isempty(c)
      scatter(c(:,1)-seg(i).BoundingBox(1),c(:,2)-seg(i).BoundingBox(2),100,'w','filled');
    end
    
    subsubplot(r1,r2,i,2,1,2);
    imagesc(seg(i).RotImage);
    daspect([1,1,1])
    axis off;
  end
  

    