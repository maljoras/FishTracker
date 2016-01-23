


LOAD = 0;
PLOT = 1;
SAVEIF = 1

if LOAD || ~exist('ft1','var')
  v =  '/home/malte/Videos/5Zebrafish_nocover_22min.avi';
  ft1 = FishTracker(v,'nfish',5,'displayif',0,'detector.fixedSize',150,'tracks.keepFullTrackStruc',true,...
                    'detector.adjustThresScale',1); 

  ft1.track([0,10]);
  
  clear reader;
  
  reader = FishVideoReader(v);
  nframes = 200;
  bkgframe = 0;
  for i = 1:nframes
    bkgframe = bkgframe + double(rgb2gray(reader.readFrame()));
  end
  bkgframe = bkgframe/nframes;
  
end


if PLOT

  lbox = 'on';
  tracks = ft1.savedTracksFull;
  
  startframe = 205; % before  crossing
  endframe = 290; % within  crossing
  
  T1 = tracks(startframe:endframe,:);
  
  msk1 = ~arrayfun(@(x)isempty(x.crossedTrackIds),T1);
  crossedIds = any(msk1,1); 

  if ~any(crossedIds)
    error('no crossing');
  end
  
  T = T1(:,crossedIds);
  msk = msk1(:,crossedIds);
  
  t = ft1.res.tracks.t(startframe:endframe,1);
  iframe =  t*ft1.videoHandler.frameRate + 1;
  crossframe = find(msk(:,1),1,'first');

  igap = 11;
  n = 3;
  iframes = [crossframe-n*igap:igap:crossframe+n*igap];
  frame = [];
  for i = 1:length(iframes)
    
    reader.setCurrentTime(t(iframes(i)));
    uframe = reader.readFrame();
    frame(:,:,i) = double(rgb2gray(uframe))-bkgframe;
  
  end
  
  
  %plot the frames
  frame = max(min(frame/125,-0.06),-0.6);
  
  figure;
  s = 0;
  r1 = 2;
  r2 = 2;
  
  
  s = s+1;
  a(s) = subplot(r1,r2,s,'align');
  
  f = mean(frame(:,:,iframes<=crossframe),3) + mean(frame(:,:,iframes>crossframe),3)*0.2;
  imagescbw(f);
  hold on;
  colororder = [rgb('green');rgb('blue')];
  set(a(s),'colororder',colororder);
  pos = reshape(cat(1,T.centroid),[size(T),2]);
  pos = convn(cat(1,pos(1,:,:),pos,pos(end,:,:)),ones(3,1)/3,'same');
  pos = pos(2:end-1,:,:);
  plot(pos(1:1:crossframe-1,:,1),pos(1:1:crossframe-1,:,2),'.','linewidth',0.5,'MarkerSize',4);
  plot(pos(crossframe+igap:end,:,1),pos(crossframe+igap:end,:,2),'--','linewidth',1,'color',0.8*[1,1,1]);
  axis xy;
  
  offs = 5; %px
  
  for i = 1:length(iframes)
    if iframes(i)>=crossframe
      break;
    end
    
    for j = 1:size(T,2)
      bbox = double(T(iframes(i),j).bbox);
      rectangle('position',bbox,'edgecolor',colororder(j,:));

      text(bbox(1) + bbox(3)/2,bbox(2) + bbox(4)+offs,sprintf('t=%1.2fs',t(iframes(i))),'fontsize',5,...
           'horiz','center','vert','bottom','color',colororder(j,:));
    end
  end
  bbox = cat(1,T(crossframe,:).bbox);
  bbox = double([min(bbox(:,1:2)),max(bbox(:,3:4))]);
  rectangle('position',bbox,'edgecolor','r','linewidth',1,'curvature',1);
  text(bbox(1) + bbox(3)/2,bbox(2) + bbox(4)+offs,sprintf('t=%1.2fs',t(crossframe)),'fontsize',5,...
       'horiz','center','vert','bottom','color','r');
  
  xlabel('x-Position [px]');
  ylabel('y-Position [px]');
  
  x = pos(1:crossframe+ceil(2.9*igap),:,1);
  y = pos(1:crossframe+ceil(2.9*igap),:,2);
  xl = minmax(x(:)').*[0.9,1.1];
  xlim( xl);
  yl = minmax(y(:)').*[0.9,1.1];
  ylim(yl);

  
  % time -position plot
  s = s+1;
  a(s) = subsubplot(r1,r2,s,2,1,1);
  hold on;
  endcrossframe = max(cat(1,T(end,:).lastFrameOfCrossing))-startframe+1;
  ihorizon = endcrossframe+ft1.opts.classifier.nFramesAfterCrossing;
  thorizon = t(ihorizon);
  
  h = [];
  set(a(s),'colororder',colororder)
  h(3) = patch(t([crossframe-1,endcrossframe+1,endcrossframe+1,crossframe-1]),[xl([1,1]),xl([2,2])],[1,0.8,0.8],...
        'edgecolor','none');
  h(5) = patch(t([endcrossframe+1,ihorizon,ihorizon,endcrossframe+1]),[xl([1,1]),xl([2,2])],[0.95,0.95,0.95],...
        'edgecolor','none');

  h(1:2) = plot(t(1:crossframe-1),pos(1:crossframe-1,:,1),'-','linewidth',1);

  hh = plot(t(crossframe:endcrossframe),pos(crossframe:endcrossframe,:,1),'-','linewidth',1,'color','r');
  h(4) = hh(1);
  h(6) = plot(t(endcrossframe+1:end),pos(endcrossframe+1:end,1,1),'--','linewidth',1,'color',0.8*[1,1,1]);
  h(7) = plot(t(endcrossframe+1:end),pos(endcrossframe+1:end,2,1),'-.','linewidth',1,'color',0.8*[1,1,1]);
  ylim([xl]);

  p1 = get(a(s),'position');
  p2 = [p1(1:2)+[0,p1(4)]-[0,p1(4)]/4,p1(3:4)/4];
  p2(1) = p2(1)+0.03;
  p2(2) = p2(2)-0.04;
  p2(3) = p2(3)*0.9;
  p2(4) = p2(4)*0.9;
  lh = legend(h([3,5:7]),{'Cross','Test','Track 1', 'Track 2'},'position',p2,'fontsize',8,'box',lbox);
  p1 = get(a(1),'position');
  p2 = [p1(1:2)+p1(3:4)-p1(3:4)/4,p1(3:4)/4];
  p2(3) = p2(3)*0.5;
  p2(4) = p2(4)*0.7;
  p2(2) = p2(2)+0.015;
  legend(a(1),h([1:2]),{'Fish A','Fish B'},'position',p2,'fontsize',8,'box',lbox)

  ylabel('x-Position [px]')
  xlim(t([1,end]));
  set(a(s),'xticklabel',[]);
  aoff = 0.08;
  shiftaxes(a(s),[0,-aoff,0,aoff]);
  
  %% assignemtn cost
  s= s+1;
  a(s) = subsubplot(r1,r2,s-1,2,1,2);
  hold on;
  set(a(s),'colororder',colororder);
  c = reshape(cat(1,T.assignmentCost),size(T));
  inv = reshape(cat(1,T.consequtiveInvisibleCount),size(T));
  c(inv>0) = NaN;
  plot(t,c,'linewidth',1)
  h = plot(t(any(inv>0,2)),1,'xk','linewidth',1);
  xlim(t([1,end]));
  shiftaxes(a(s),[0,0,0,-aoff]);
  
  p1 = get(a(s),'position');
  p2 = [p1(1:2)+[0,p1(4)]-[0,p1(4)]/4,p1(3:4)/4];
  p2(1) = p2(1)+0.034;
  p2(2) = p2(2)-0.015;
  p2(3) = p2(3)*0.9;
  p2(4) = p2(4)*0.9;

  legend(h(1),{'invisible'},'position',p2)
  xlabel('Time [sec]')
  ylabel('Cost');
  
  %%classification
  s = s+1;
  a(s) = subsubplot(r1,r2,s-1,1,1,1);

  feat = {};
  pre = {};
  prob = {};
  h = [];
  for i = 1:size(T,2)
    feat{i} = T(end,i).batchFeatures(1:ft1.opts.classifier.nFramesAfterCrossing,:);
    prob{i} = ft1.fishClassifier.predict(feat{i});
    prob{i} = prob{i}(:,crossedIds);
    pre{i} = ft1.fishClassifier.preprocess(feat{i});
  end
  
  
  hold on;  
  comps = [1,2];
  mu = ft1.fishClassifier.mu(crossedIds,comps);
  sigma = ft1.fishClassifier.Sigma(comps,comps,crossedIds);
  set(a(s),'colororder',colororder)
  sym = 'xo><';

  for i = 1:length(pre)
    [~,id] = max(prob{i},[],2);
    for j = 1:size(pre{i},1);

      hh = plot(pre{i}(j,comps(1)),pre{i}(j,comps(2)),sym(id(j)),'color',colororder(i,:));
    end
    h(i) = plotGauss(mu(i,:)',sigma(:,:,i),'-','color',colororder(i,:));
  end

  xl1 = [-60,60];
  yl1 = [-50,50];
  xlim(xl1);
  ylim(yl1)
  xlabel('1st GMM component');
  ylabel('2nd GMM component');
  g = [];
  for i = 1:length(pre);
    g(i) = plot(xl1(1)-10,yl1(1)-10,sym(i),'color',0.4*[1,1,1]);
  end

  p1 = get(a(s),'position');
  p2 = [p1(1:2)+[0,p1(4)]-[0,p1(4)]/4,p1(3:4)/4];
  p2(1) = p2(1)+0.047;
  p2(2) = p2(2)-0.015;
  p2(3) = p2(3)*0.9;
  p2(4) = p2(4)*0.9;

  legend([g,h],{'Test: track 1','Test: track 2','GMM fish A','GMM fish B'},'position',p2,'fontsize',8,'box',lbox)
  
  
  xlim(xl1);
  ylim(yl1)


  %%result
  
    
  s = s+1;
  a(s) = subsubplot(r1,r2,s-1,1,1,1);
  
  f = mean(frame(:,:,iframes>crossframe),3) +  0.2*mean(frame(:,:,iframes<=crossframe),3);
  imagescbw(f);
  hold on;
  set(a(s),'colororder',colororder);
  plot(pos(:,:,1),pos(:,:,2),'.','linewidth',0.5,'MarkerSize',4);
  axis xy;
  
  for i = 1:length(iframes)
    if iframes(i)<=crossframe
      continue;
    end
    
    for j = 1:size(T,2)
      bbox = double(T(iframes(i),j).bbox);
      rectangle('position',bbox,'edgecolor',colororder(j,:));

      text(bbox(1) + bbox(3)/2,bbox(2) + bbox(4) +offs,sprintf('t=%1.2fs',t(iframes(i))),'fontsize',5,...
           'horiz','center','vert','bottom','color',colororder(j,:));
    end
  end
  bbox = cat(1,T(crossframe,:).bbox);
% $$$   bbox = double([min(bbox(:,1:2)),max(bbox(:,3:4))]);
% $$$   rectangle('position',bbox,'edgecolor','r','linewidth',1,'curvature',1);
% $$$   text(bbox(1) + bbox(3)/2,bbox(2) + bbox(4),sprintf('t=%1.2fs',t(crossframe)),'fontsize',6,...
% $$$        'horiz','center','vert','bottom','color','r');
% $$$   
  xlabel('x-Position [px]');
  ylabel('y-Position [px]');
  
  xlim( xl);
  ylim(yl);
% $$$   x = pos(crossframe-2*igap:end,:,1);
% $$$   y = pos(crossframe-2*igap:end,:,2);
% $$$   xl = minmax(x(:)').*[0.9,1.1];
% $$$   xlim( xl);
% $$$   yl = minmax(y(:)').*[0.9,1.1];
% $$$   ylim(yl);

  set(a,'fontsize',10)

  p1 = get(gcf,'position');
  p1(3:4) = p1(3:4)*1.35;
  set(gcf,'position',p1);

  set(a,'clipping','on');
  labelsubplot(gcf,'BCADE');
  
  if SAVEIF
    exportfig(gcf,'/home/malte/work/writings/papers/fishtracking/figs/crossing.eps','color','rgb','FontSizeMin',5,'FontMode','scaled')
    
    saveas(gcf,'/home/malte/work/writings/papers/fishtracking/figs/crossing.fig');
  end

end
