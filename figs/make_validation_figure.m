LOAD = 1;
PLOT = 1;
SAVEIF = 1;
if LOAD

  %id = load('~/work/projects/zebra/trackering/trajectories.mat');
  id = load('~/Videos/segm/trajectories.mat');
  idres.pos = permute(id.trajectories,[1,3,2]);

  vid = '/home/malte/Videos/5Zebrafish_nocover_22min.avi';  
  ft = FishTracker(vid,'detector.adjustThresScale',1,'nfish',5,'classifier.reassignProbThres',0.3,...
                   'classifier.handledProbThres',0.4);
  ft.addSaveFields('firstFrameOfCrossing', 'lastFrameOfCrossing');
  ft.setDisplay(1);
  ft.setDisplay('switchFish',1,'tracks',1);
  
  tic;
  ft.track();
  toc
end


if PLOT
  
  figure;
  
  r1 = 3;
  r2 = 1;
  s = 0;
  a = [];
  
  nfish = ft.nfish;
  ftres = ft.getTrackingResults();
  ftresnan = ft.getTrackingResults(1);

  dist = zeros(nfish);

  for i = 1:nfish
    for j = 1:nfish
      dist(i,j) = nanmean(sqrt(sum((ftres.pos(:,:,i) - idres.pos(1:end,:,j)).^2,2)));
    end
  end
  assignments = assignDetectionsToTracks(dist,1e3);

  ftpos = ftres.pos(:,:, assignments(assignments(:,1),2));
  ftposnan = ftresnan.pos(:,:, assignments(assignments(:,1),2));
  idpos = idres.pos(1:end,:,:);
 
  %distances
  s = s+1;
  a(end+1) = subplot(r1,r2,s,'align');
  
  t = ftres.tracks.t(:,1);
  d = sqrt(sum((ftpos(:,:,:) - idpos(:,:,:)).^2,2));
  
  errorbarpatch(t,nanmean(d,3),stderr(d,3));
  
  set(a(s),'fontsize',8);
  ylabel('Avg distance [px]','fontsize',10);
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
  
  legend({'FishTracker','idTracker'},'location','NorthWest','fontsize',8)
  box off;
  

  %crossings
  s = s+1;
  a(end+1) = subplot(r1,r2,s,'align');
  

  cross = diff(ftres.tracks.lastFrameOfCrossing-ftres.tracks.firstFrameOfCrossing)>0;
  ccross = conv2(sum(double(cross),2),ones(nconv,1)/nconv,'same');
  plot(t(1:end-1),ccross,'k')
  set(a(s),'fontsize',8);
  ylabel(sprintf('Avg. # of \ncrossing events'),'fontsize',10);
  xlabel('Time [sec]','fontsize',10);
  xlim(t([1,end]));
  box off;
  
  b = labelsubplot();
  shiftaxes(b,[0.02])
  
  if SAVEIF
    exportfig(gcf,'/home/malte/work/writings/papers/fishtracking/figs/validation.eps',...
              'color','rgb','FontSizeMin',5,'FontMode','scaled','LineWidthMin',0.75)
    saveas(gcf,'/home/malte/work/writings/papers/fishtracking/figs/validation.fig');
  end

  
end
