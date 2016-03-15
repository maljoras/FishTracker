LOAD = 0;
PLOT = 1;
SAVEIF = 1;
if LOAD
  load ~/data/zebra/videos/longterm/Blongterm11.mat
  reader = fish.core.FishVideoReader(ft.videoFile);
  
end


if PLOT
  
  figure;
  s = 0;
  a = [];
  r1 = 2;
  r2 = 2;

  % orginal
  s = s+1;
  a(end+1) = subplot(r1,r2,s,'align');

  iframe = [15000];
  t = ft.res.tracks.t(:,1);
  frame = [];
  reader.setCurrentTime(t(iframe));
  frame = reader.readFrame();
  
  image(frame);
  axis xy;
  hold on;
  bbox = double(squeeze(ft.res.tracks.bbox(iframe,:,:)));
  cols = [rgb('red');rgb('blue');rgb('green')];

  for i = 1:size(bbox,1);
    rectangle('position',bbox(i,:),'edgecolor',cols(i,:),'linewidth',0.5);
  end
  set(a(s),'colororder',cols)

  iback = 300;
  x = squeeze(ft.res.pos(iframe-iback:iframe,1,:));
  y = squeeze(ft.res.pos(iframe-iback:iframe,2,:));
  plot(x,y,'linewidth',0.5);
  xlabel('x-Position [px]')
  ylabel('y-Position [px]');
  
  % velocity  
  s = s+1;
  a(end+1) = subplot(r1,r2,s,'align');
  
  posx = ft.res.tracks.centroid(:,:,1);
  posy = ft.res.tracks.centroid(:,:,2);
  inv = ft.res.tracks.consecutiveInvisibleCount>0;
  posx(inv) = NaN;
  posy(inv) = NaN;
  d12 = sqrt((posx(:,1) - posx(:,2)).^2 + (posy(:,1) - posy(:,2)).^2);
  d13 = sqrt((posx(:,1) - posx(:,3)).^2 + (posy(:,1) - posy(:,3)).^2);
  d23 = sqrt((posx(:,2) - posx(:,3)).^2 + (posy(:,2) - posy(:, 3)).^2);

  vel1 = abs(diff(posx)) + abs(diff(posy));
  msk =([vel1;zeros(1,3)]>50) | ([zeros(1,3);vel1]>50);
  posx(msk) = NaN;
  posy(msk) = NaN;
  d12(any(msk,2)) = NaN;
  d23(any(msk,2)) = NaN;
  d13(any(msk,2)) = NaN;
  
  %posx = nanmean(cat(3,posx,circshift(posx,[1,0]),circshift(posx,[-1,0])),3);
  %posy = nanmean(cat(3,posy,circshift(posy,[1,0]),circshift(posy,[-1,0])),3);
  vel = abs(diff(posx)) + abs(diff(posy));
  mvel = nanmean(cat(3,vel,circshift(vel,[1,0]),circshift(vel,[-1,0]),circshift(vel,[2,0]),circshift(vel,[-2,0])),3);
  d = nanmin(cat(3,d13,d23,d12),[],3);
  md = nanmean(cat(3,d,circshift(d,[1,0]),circshift(d,[-1,0]),circshift(d,[2,0]),circshift(d,[-2,0])),3);

  tspeed = 15; % 20px per frame
  tmvel = mvel>tspeed;
  inter = [];
  inter(:,1) = tmvel(:,1) & tmvel(:,2);
  inter(:,2) = tmvel(:,1) & tmvel(:,3);
  inter(:,3) = tmvel(:,2) & tmvel(:,3);

  inter = double(circshift(inter,[1,0]) & circshift(inter,[-1,0]) & inter);
  sdfinter = [];
  for i = 1:size(inter,2)
    sdfinter(:,i) = spks2sdf(find(inter(:,i)),1,[0,size(inter,1)],4*3600*2);
  end
  sdfinter =  sdfinter(1:1000:end,:)*3600;
  tt = t(1:1000:end)/3600;
  hold on;
  h = [];
  h(1) = plot(tt,sdfinter(:,1),'-','linewidth',1,'color',cols(1,:));
  plot(tt,sdfinter(:,1),'--','linewidth',1,'color',cols(2,:));
  h(3) = plot(tt,sdfinter(:,2),'-','linewidth',1,'color',cols(3,:));
  plot(tt,sdfinter(:,2),'--','linewidth',1,'color',cols(1,:));
  h(2) = plot(tt,sdfinter(:,3),'-','linewidth',1,'color',cols(2,:));
  plot(tt,sdfinter(:,3),'--','linewidth',1,'color',cols(3,:));
  xlabel('Time [h]');
  ylabel('Freq. of interaction [#/h]');
  xlim([0,t(end)/3600])
  
  [lhx,lhy,lhz,lhq] = legend(h,{'Fish ID_1','Fish ID_2','Fish ID_3'},'fontsize',8);
  
  % probmaps
  s = s+1;
  bb = subsubplot(r1,r2,s,1,1,1);
  set(bb,'visible','off');
  range = 0:3600:4*3600;%ft.videoHandler.duration;
  n = length(range)-1;
  [rr1,rr2] = fish.helper.getsubplotnumber(n);
  for i = 1:n
    aa(i) = subsubplot(r1,r2,s,rr1,rr2,i);
    ft.plotProbMap([range(i),range(i+1)]);
    title(sprintf('%dh',i),'fontsize',8);
    set(gca,'xticklabel',[],'yticklabel',[]);
    xlabel('');
    ylabel('');
    axislabel off;
  end
  
  
  % domains
  s = s+1;
  a(end+1) = subplot(r1,r2,s,'align');
  ft.plotDomains();
  colormap([0.1,0.1,0.1;cols])
  f = smallcolorbar();
  set(f,'ytick',[linspace(1+3/8,4-3/8,4)])
  set(f,'yticklabel',{'None','ID_1','ID_2','ID_3'})
  shiftaxes(f,[0.125])

  
  set(a,'fontsize',8);
  p1 = get(gcf,'position');
  p1(4) = p1(4)*1.1;
  set(gcf,'position',p1);
  
  labelsubplot(gcf,'BACD');
  
  if SAVEIF
    exportfig(gcf,'/home/malte/work/writings/papers/fishtracking/figs/longterm.eps',...
              'color','rgb','FontSizeMin',5,'FontMode','scaled','LineWidthMin',0.75)
    saveas(gcf,'/home/malte/work/writings/papers/fishtracking/figs/longterm.fig');
  end

  
end
