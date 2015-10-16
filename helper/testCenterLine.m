function testCenterLine(img)

  
  
  % opening and closing
  img = cv.dilate(img);
  img = cv.erode(img);
  img = cv.erode(img);
  img = cv.dilate(img);
  
  F = cv.distanceTransform(img);
  L = cv.Laplacian(F,'KSize',5);
  L1 = L;
  L1(L1>min(L1(:))*2/3) = 0;
  

  
  if size(img,1)<size(img,2)
    [x,y] = find(L1);  
  else
    [y,x] = find(L1');
  end
  
  sz = size(img);

  %x = conv(x,ones(10,1)/10,'valid');  
  %y = conv(y,ones(10,1)/10,'valid');  
  t = linspace(0,1,length(x));
  
  probe = floor([0:0.25:1] * (length(x)-1))+1;
  xi = x(probe);
  yi = y(probe);
  xii = x(probe);
  yii = y(probe);

  %xi = interp1(t,x,[0:0.25:1],'linear');
  %yi = interp1(t,y,[0:0.25:1],'linear');
  %xii = min(max(round(xi),1),sz(1));
  %yii = min(max(round(yi),1),sz(2));

  if F(xii(1),yii(1)) + F(xii(2),yii(2)) < F(xii(end),yii(end)) + F(xii(end-1),yii(end-1))
    % thinning. reverse
    xi = xi(end:-1:1);
    yi = yi(end:-1:1);
  end
  

  v = [xi(1)-xi(2),yi(1)-yi(2)];
  d = norm(v);
  a = v(2)/d;
  b = v(1)/d;
  
  szin = sz([2,1]);
  %szout = [szin(1)*abs(a)+szin(2)*abs(b),szin(1)*abs(b)+szin(2)*abs(a)];
  szout = [120,80];

  cx = yi(1);
  cy = xi(1);
  offs = [szout(1)-10,szout(2)/2]-[cx,cy];
  T = [ a, b, (1-a)*cx-b*cy + offs(1); -b, a, b*cx+(1-a)*cy + offs(2) ];
  warpimg = cv.warpAffine(img,T,'DSize',szout,'Interpolation','Nearest','BorderValue',0);

  % rotate

  xgrad = cv.Sobel(F,'XOrder',1,'YOrder',0, 'KSize',3);
  ygrad = cv.Sobel(F,'XOrder',0,'YOrder',1, 'KSize',3);


  
  clf;
  s = 0;
  r1 = 2;
  r2 = 3;
  for z = {F,L,L1,warpimg,xgrad,img};
    s = s+1;
    subplot(r1,r2,s,'align');
   
    if s==r1*r2
      imagesc([1,sz(2)],[1,sz(1)],z{1});
      xlim([1,sz(2)])
      ylim([1,sz(1)]);

    else
      imagesc(z{1});
    end
    
    axis xy;;
    daspect([1,1,1])
    axis off;
  end
  
  hold on;
  scatter(yi,xi,60,gray(5),'o','filled');
  
