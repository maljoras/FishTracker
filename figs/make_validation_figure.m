LOAD = 0;
PLOT = 1;
SAVEIF = 0;
if LOAD

  %id = load('~/work/projects/zebra/trackering/trajectories.mat');
  id = load('~/Videos/segm/trajectories.mat');
  idres.pos = permute(id.trajectories,[1,3,2]);

  vid = '/home/malte/Videos/5Zebrafish_nocover_22min.avi';  
  ft = fish.Tracker(vid,'detector.adjustThresScale',1,'nfish',5,'detector.fixedSize',150);

  ft.addSaveFields('firstFrameOfCrossing', 'lastFrameOfCrossing');
  ft.setDisplay('switchFish',1,'tracks',1);
  ft.setDisplay(0);
  
  tic;
  ft.track();
  toc
  %ft.getPosFromDag();
end


if PLOT
  
  figure;
  
  r1 = 4;
  r2 = 1;
  s = 0;
  a = [];
  
  nfish = ft.nfish;
  ftres = ft.getTrackingResults();
  ftresnan = ft.getTrackingResults(1);

  dist = zeros(nfish);

  for i = 1:nfish
    for j = 1:nfish
      dist(i,j) = nanmean(sqrt(sum((ftres.pos(:,:,i) - idres.pos(1:end-1,:,j)).^2,2)));
    end
  end
  assignments = assignDetectionsToTracks(dist,1e3);

  ftpos = ftres.pos(:,:, assignments(assignments(:,1),2));
  ftposnan = ftresnan.pos(:,:, assignments(assignments(:,1),2));
  idpos = idres.pos(1:end-1,:,:);
 
  %distances
  s = s+1;
  a(end+1) = subplot(r1,r2,s,'align');
  
  t = ftres.tracks.t(:,1);
  d = sqrt(sum((ftpos(:,:,:) - idpos(:,:,:)).^2,2))/ft.fishlength;
  
  MAXDISTANCE = 0;
  if MAXDISTANCE
    plot(t,nanmax(d,[],3));
  else
    errorbarpatch(t,nanmean(d,3),stderr(d,3));
  end
  
  hold on;
  plot(t([1,end]),[1,1],'r:');
  
  set(a(s),'fontsize',8);
  if MAXDISTANCE
    ylabel(sprintf('Max distance\n[Fish length]'),'fontsize',10);
  else
    ylabel(sprintf('Avg. distance\n[Fish length]'),'fontsize',10);
  end
  
  xlim(t([1,end]));
  box off;
    
  %detections
  s = s+1;
  a(end+1) = subplot(r1,r2,s,'align');
  nconv = 75;
  ftnan = conv(sum(isnan(ftposnan(:,1,:)),3),ones(nconv,1)/nconv,'same');
  idnan = conv(sum(isnan(idpos(:,1,:)),3),ones(nconv,1)/nconv,'same');
  plot(t,[ftnan,idnan]*100/ft.nfish);
  xlim(t([1,end]))
  set(a(s),'fontsize',8);

  ylabel(sprintf('Avg. lost\n tracks [%%]'),'fontsize',10);
  
  legend({'fish.Tracker','idTracker'},'location','NorthWest','fontsize',8)
  box off;
  

  %crossings
  s = s+1;
  a(end+1) = subplot(r1,r2,s,'align');
  

  cross = diff(ftres.tracks.lastFrameOfCrossing-ftres.tracks.firstFrameOfCrossing)>0;
  ccross = conv2(sum(double(cross),2),ones(nconv,1)/nconv,'same');
  plot(t(1:end-1),ccross,'k')
  set(a(s),'fontsize',8);
  ylabel(sprintf('Avg. # of \ncrossing events'),'fontsize',10);
  xlim(t([1,end]));
  box off;


  % probability
  s = s+1;
  a(end+1) = subplot(r1,r2,s,'align');
  clp = ftres.tracks.classProb;
  dclp = zeros(size(clp,1),size(clp,2));
  mclp = zeros(size(clp,1),size(clp,2));
  for i =1:size(clp,2)
    mclp(:,i) = max(clp(:,i,setdiff(1:size(clp,2),i)),[],3);
    dclp(:,i) = clp(:,i,i);
  end
  
  percent = sum(sum(mclp>dclp))/numel(mclp);
  %fprintf('%1.1f%% correct.\n',(1-percent)*100);
  
  mclp = conv2(mclp,ones(nconv,1)/nconv,'same');
  dclp = conv2(dclp,ones(nconv,1)/nconv,'same');

  h1 = errorbarpatch(t,mean(dclp,2),stderr(mclp,2),'color',rgb('green'));
  hold on;
  h2 = errorbarpatch(t,mean(mclp,2),stderr(dclp,2),'color',rgb('orange'));
  xlim(t([1,end]));
  box off;
  ylabel(sprintf('Class\nprobability'),'fontsize',10);
  legend([h1(1),h2(1)],{'Avg.','Other'},'location','NorthWest','fontsize',8)
  set(a(end),'fontsize',8)
  
  
  
  xlabel('Time [sec]','fontsize',10);
  b = labelsubplot(gcf);
  shiftaxes(b,[0.02])
  
  
  
  if SAVEIF
    exportfig(gcf,'/home/malte/work/projects/zebra/validation.eps',...
              'color','rgb','FontSizeMin',5,'FontMode','scaled','LineWidthMin',0.75)
    saveas(gcf,'/home/malte/work/projects/zebra/validation.fig');
  end

  
end
