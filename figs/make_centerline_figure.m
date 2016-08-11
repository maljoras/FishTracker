

LOAD = 0;
PLOTCENTERLINE = 0;
PLOTIDENTITY = 1;
SAVEIF = 1

if LOAD || ~exist('T1','var')
  v =  '/home/malte/Videos/5Zebrafish_nocover_22min.avi';
  T1 = xy.Tracker(v,'nindiv',5,'displayif',0,'detector.fixedSize',150,'tracks.keepFullTrackStruc',true,...
                    'detector.adjustThresScale',1); 

  T1.track([0,20]);
  tracks = T1.savedTracksFull;
end


if PLOTCENTERLINE
  tracks = T1.savedTracksFull;
  
  idx = 180;
  id = 1;
  T = tracks(idx,id);
  l = T1.opts.detector.fixedSize;
  bbox = double(T.bbox);
  rbbox = bbox;
  rbbox(1:2) = bbox(1:2) - T.segment.Centroid + l/2 +1;
  rbbox = round(rbbox);

  img = false(l,l);
  sz = size(T.segment.Image);
  img(rbbox(2):rbbox(2)+sz(1)-1,rbbox(1):rbbox(1)+sz(2)-1) = T.segment.Image;

  cl = T.centerLine;
  cl(:,1) = cl(:,1) - bbox(1) + rbbox(1) ;
  cl(:,2) = cl(:,2) - bbox(2) + rbbox(2) -1;
    
  r = testCenterLine(img);
  
  clf;
  a = [];
  s = 0;
  r1 = 2;
  r2 = 3;
  xl = [20,120];
  yl = [20,100];
  
  % plot the initial image
  s = s+1;
  a(s) = subplot(r1,r2,s,'align');
  xy.helper.imagescbw(double(T.segment.FilledImageFixedSize));
  hold on;
  rectangle('position',rbbox,'edgecolor','r');
  xlabel('x [px]')
  ylabel('y [px]')
  title(sprintf('Backgroung\nsubtracted image'));
  
  
  % BW mask
  s = s+1;
  a(s) = subplot(r1,r2,s,'align');
  xy.helper.imagescbw(double(r.img));
  title('Detected blob');
  hold on;
  R = R2d(-T.segment.Orientation/180*pi + pi/2);
  S = diag(T.segment.Size([1,2]))/2;
  k = 1;
  n = 100;
  mu = (T.segment.Centroid - bbox([1,2]) + rbbox([1,2]) )';
  % Compute the points on the surface of the ellipse.
  t = linspace(0, 2*pi, n);
  u = [cos(t); sin(t)];
  w = (k * R * S) * u;
  z = repmat(mu, [1 n]) + w;
  plot(z(1, :), z(2, :),'linewidth',1);
  plot(mu(1),mu(2),'xb','linewidth',1);



  
  
  % distance tracnform
  s = s+1;
  a(s) = subplot(r1,r2,s,'align');
  imagesc(r.F);
  title('Distance transform');
  f = smallcolorbar(a(s));
  shiftaxes(f,[0.12]);
  ylabel(f,'Dist. to blob border [px]');
  
  % Laplacian
  s = s+1;
  a(s) = subplot(r1,r2,s,'align');
  imagesc(r.L);
  title('Laplace transform')
  f1 = smallcolorbar(a(s));
  shiftaxes(f1,[0.12]);
  
  m = min(r.L(:))*0.66;
  g = newaxes(f1);
  set(g,'tag','NoLabel');
  hold on;
  plot(g,[xlim(f1)],[m,m],'r--','linewidth',1)
  
  % Laplacian thresholded
  s = s+1;
  a(s) = subplot(r1,r2,s,'align');
  xy.helper.imagescbw(1-double(r.L1~=0));
  title(sprintf('Thresholded\nLaplacian'))

  % Center line points
  s = s+1;
  a(s) = subplot(r1,r2,s,'align');
  xy.helper.imagescbw(double(T.segment.FilledImageFixedSize));
  hold on;
  %y = r.xi-r.border;
  %x = r.yi-r.border;
  x = cl(:,1);
  y = cl(:,2);
  plot(x,y,'.r','linewidth',2,'MarkerSize',10);
  %plot(cl(:,1),cl(:,2),'.g','linewidth',2,'MarkerSize',10);
  xoffs = -2;
  yoffs = -7;
  for i = 1:length(r.yi)
    if ~mod(i,2)
      text(x(i) +xoffs,y(i) + yoffs,sprintf('$\\mathbf{c}_%d$',i), 'Interpreter','Latex',...
           'Horiz','left','vert','middle','color','k','fontsize',10)
    else
      text(x(i) -xoffs,y(i) - yoffs,sprintf('$\\mathbf{c}_%d$',i), 'Interpreter','Latex',...
           'Horiz','right','vert','middle','color','k','fontsize',10)
      
    end
  end
  
  title('Center line')
  xlabel('x [px]')
  ylabel('y [px]')
  set(a(s),'yaxislocation','right');
  
  set(a,'xlim',xl);
  set(a,'ylim',yl);
  set(a(2:end-1),'xticklabel',[],'yticklabel',[])
  set(a,'fontsize',8);
  set(f,'fontsize',8)

  b = labelsubplot(gcf,[],14);
  shiftaxes(b,[0.025,0.015]);
  
  
  if SAVEIF
    print -depsc2 -painters /home/malte/work/writings/papers/fishtracking/figs/centerline.eps
    saveas(gcf,'/home/malte/work/writings/papers/fishtracking/figs/centerline.fig');
  end
  
end



if PLOTIDENTITY
  clf;
  tracks = T1.savedTracksFull;
  nf = 5;
  ni = 7;
  istart = 504;
  a = [];
  r1 = 2;
  r2 = 2;
  l = T1.opts.detector.fixedSize;
  
  a(1) = subplot(r1,r2,1,'align');
  s =0;
  for i = istart:istart+ni-1
    s = s+1;
    for j = 1:nf
      seg = tracks(i,j).segment;
      
      img = double(cv.copyMakeBorder(seg.FilledImageFixedSize,3,3,3,3,'BorderType','Constant','Value',seg.mback(1)));
      imagesc([s-0.5,s+0.5],[j-0.5,j+.5],img);

      hold on
      bbox = double(seg.BoundingBox);
      rbbox = bbox;
      rbbox(1:2) = bbox(1:2) - seg.Centroid + l/2 + 1;
      
      rbbox = rbbox/l;
      rbbox(1) = rbbox(1) + s-0.5;
      rbbox(2) = rbbox(2) + j-0.5;
      rectangle('position',rbbox,'edgecolor',[0.9,0.6,0.4]);
      
    end
  end
  
  % blob
  rsz = [110,110];
  a(2) = subplot(r1,r2,2,'align');
  s = 0;
  for i = istart:istart+ni-1
    s = s+1;
    for j = 1:nf
      seg = tracks(i,j).segment;
      
      img = seg.Image;
      sz = size(img);
      w = rsz-sz;
      rimg = cv.copyMakeBorder(img,floor(w(1)/2),ceil(w(1)/2),floor(w(2)/2),ceil(w(2)/2),...
                               'BorderType','Constant','Value',0);
      xy.helper.imagescbw([s-0.5,s+0.5],[j-0.5,j+.5],double(~rimg));
      hold on
      
      % plot center line
      center = double(seg.BoundingBox(1:2)) + double(seg.BoundingBox(3:4))/2;
      cl = bsxfun(@rdivide,bsxfun(@plus,seg.CenterLine, - center + rsz/2+1),rsz);
      cl(:,1) = cl(:,1) + s-0.5;
      cl(:,2) = cl(:,2) + j-0.5;
      
      plot(cl(1,1),cl(1,2),'.r','Markersize',7,'linewidth',1);
      plot(cl(:,1),cl(:,2),'-','color',[0.8,0.2,0.2],'linewidth',0.1);
    end
  end
    
  % rotated
  a(3) = subplot(r1,r2,3,'align');
  s = 0;
  for i = istart:istart+ni-1
    s = s+1;
    for j = 1:nf
      seg = tracks(i,j).segment;
      
      img = seg.RotFilledImage;
      sz = size(img);
      w = rsz-sz;
      rimg = cv.copyMakeBorder(img,floor(w(1)/2),ceil(w(1)/2),floor(w(2)/2),ceil(w(2)/2),...
                               'BorderType','Constant','Value',seg.mback(1));
      rimg = cv.copyMakeBorder(rimg,2,2,2,2,'BorderType','Constant','Value',seg.mback(1)*0.95);
      imagesc([s-0.5,s+0.5],[j-0.5,j+.5],double(rimg));
      hold on
      
      % feature box
      fsz = size(seg.IdentityFeature);
      fbox = [rsz(1) - floor(w(2)/2) - fsz(1) ,rsz(2)/2 - fsz(2)/2+1,fsz];
      fbox([1,3]) = fbox([1,3])./rsz(1);
      fbox([2,4]) = fbox([2,4])./rsz(2);
      fbox(1) = fbox(1) + s-0.5;
      fbox(2) = fbox(2) + j-0.5;

      rectangle('position',fbox,'edgecolor',[0.9,0.6,0.4]);
    end
  end
  
  % fish feature
  a(4) = subplot(r1,r2,4,'align');
  s = 0;
  for i = istart:istart+ni-1
    s = s+1;
    for j = 1:nf
      seg = tracks(i,j).segment;
      
      img = seg.IdentityFeature;
      sz = size(img);
      b = 2;
      rimg = cv.copyMakeBorder(img,0,0,floor((sz(1)-sz(2))/2),ceil((sz(1)-sz(2))/2),...
                               'BorderType','Constant','Value',seg.mback(1));

      rimg = cv.copyMakeBorder(rimg,b,b,b,b,'BorderType', 'Constant', 'Value',seg.mback(1)*0.95);
      rimg = flipud(rimg);
      imagesc([s-0.5,s+0.5],[j-0.5,j+.5],double(rimg));
      hold on
      
    end
  end


  % formating
  border = 0.3;
  for i = 1:length(a)
    axes(a(i));
    axis xy;
    daspect([1,1,1])
    if ~mod(i-1,r2)
      ylabel('Identity [#]');
    end
    if i>r1*(r2-1)
      xlabel('Frame [#]');
    end
    

    xlim([border,ni+1-border])
    ylim([border,nf + 1-border])
  
    set(a(i),'xtick',1:ni);
    set(a(i),'ytick',1:nf);
  end
  title(a(1),'Raw image patches');
  title(a(2),'Detected blobs');
  title(a(3),'Rotated images')
  title(a(4),'Identity features');
  
  set(a,'fontsize',10);
  b = labelsubplot(gcf,[],14);
  shiftaxes(b,[0.02,0.01])

  if SAVEIF
    print -depsc2 -painters /home/malte/work/writings/papers/fishtracking/figs/identity.eps
    saveas(gcf,'/home/malte/work/writings/papers/fishtracking/figs/identity.fig');
  end

end
