LOAD = 0;
PLOT = 1;
SAVEIF = 1;
COMPUTE = 0;

NEWTRACK = 0
LOADALL = 0;

if LOAD || ~exist('ftstm','var')

  dname = '~/data/zebra/videos/halfplane';
  fname = 'five2onefish2';
  

  vid = [dname filesep fname '.avi'];
  ftstm = load([dname filesep fname '.mat']);

  stmpos = ftstm.ft.res.pos;
  
end

if LOADALL
  dname = '~/data/zebra/videos/halfplane';
  fnames = {{'onefish1','onefish2','onefish3'},{'twofish1','twofish2','twofish3'},{'testThreefish','testThreefish2','testThreefish3','testThreefish4','testThreefish5'},{},{'fivefish1','fivefish2','fivefish3'}};
  
  for i = 1:length(fnames)

    for j = 1:length(fnames{i})
      ftstmall{i}{j} = load([dname filesep fnames{i}{j} '.mat']);
    end
  end
    
end

  


if COMPUTE && NEWTRACK % for tracking again
  
  ft = FishTracker(vid,'detector.adjustThresScale',1,'nfish',2,'detector.inverted',1);

  ft.addSaveFields('firstFrameOfCrossing', 'lastFrameOfCrossing');
  ft.setDisplay(0);
  ft.setDisplay('switchFish',0,'tracks',0,'level', 3);
  
  tic;
  ft.track();
  toc;
end

if NEWTRACK && ~exist('stmres1','var')
  %% get the adjusted stminfo 
  if ~ft.res.tracks.t(1,1)
    tfile = dlmread([vid '.txt'])  ;
    t = tfile(1,3) + ft.res.tracks.t(:,1);
  else
    t =  ft.res.t(:,1);
  end
  tc = ceil(t*1e4)/1e4;
  tf = floor(t*1e4)/1e4;
  
  stmres = ftstm.ft.res.tracks;
  ftres = ft.res.tracks;
  ftpos = ft.res.pos;


  tstm = stmres.t(:,1);
  tstm = round(tstm*1e4)/1e4;

  %tpos is the position of tstm in t of the video
  [~,tpos1] = ismember(tstm,tf);
  [~,tpos2] = ismember(tstm,tc);
  tpos = max(tpos1,tpos2);
  % some frames might be lost during saving (but used in the tracking), so tpos has still some zeros

  for f = fieldnames(stmres)'
    idx = ~~tpos;
    % will be somewhat weired to interp for vectors but do not care..
    stmres1.(f{1}) = interp1(tstm(idx),double(stmres.(f{1})(idx,:,:,:,:)),t,'nearest',NaN);
  end

  stmpos1 = interp1(tstm(idx),stmpos(idx,:,:),t,'nearest',NaN);

  
  % make ordering correct
  dist = [];
  for i = 1:ft.nfish
    for j = 1:ft.nfish
      dist(i,j) = nanmean(sqrt(sum((ftpos(:,:,i) - stmpos1(:,:,j)).^2,2)));
    end
  end

  assignments = assignDetectionsToTracks(dist,1e3);
  change = assignments(assignments(:,1),2);
  ftpos  = ftpos(:,:,change);
  for f = fieldnames(ftres)'
    ftres.(f{1}) = ftres.(f{1})(:,change,:,:,:);
  end
  ftres.stmInfo = stmres1.stmInfo;

  % check tracking differences
  %figure;
  %plot(t,squeeze(stmpos1(:,1,:)-ftpos(:,1,:)));
  
end


if PLOT
  figure;

  if NEWTRACK
    res = ftres; % take new
    pos = ftpos;
    tt = t;
  else
    %res = ftstm.ft.res.tracks;
    %pos = stmpos;
    %tt = ftstm.ft.res.tracks.t(:,1);

    ftstm2 = ftstmall{5}{1};
    ftstm2.ft.videoFile
    pos = ftstm2.ft.res.pos;
    tt = ftstm2.ft.res.tracks.t(:,1);    
    res = ftstm2.ft.res.tracks;
  end

  c = res.consequtiveInvisibleCount;
  cc = permute(cat(3,c,c),[1,3,2]);
  pos(cc>0) = NaN;
  
  r1  = 2;
  r2 = 1;
  s = 0;
  a = [];
  
  s = s+1;
  a(end+1) = subplot(r1,r2,s,'align');


  if size(pos,3)<=2 % SOMEHOW FOR VIDEO WITH 1 or 2 fish stimulus
                    % was reversed!!
    left = 2;
    right = 1;
  else
    left = 1;
    right = 2;
  end
  
  tt = tt/60;

  info = res.stmInfo(:,1,2);
  blackdinfo = [diff(info==right);0];
  whitedinfo = [diff(info==left);0];

  xhalf = ftstm.ft.videoHandler.frameSize(2)/2;
  xend = ftstm.ft.videoHandler.frameSize(2);
  gray = [1,1,1]*0.6;
  for idx = find(blackdinfo>0)'
    tstart = tt(idx);
    idx1 = find(blackdinfo<0);
    tend = tt(idx1(find(idx1>idx,1,'first')));
    patch([tstart,tend,tend,tstart,tstart],[1,1,xhalf,xhalf,1],gray,'edgecolor','none');
  end
  for idx = find(whitedinfo>0)'
    tstart = tt(idx);
    idx1 = find(whitedinfo<0);
    tend = tt(idx1(find(idx1>idx,1,'first')));
    patch([tstart,tend,tend,tstart,tstart],[xend,xend,xhalf,xhalf,xend],gray,'edgecolor','none');
  end

  hold on;
  set(a(end),'colororder',jet(size(res.t,2)))
  p = squeeze(pos(:,1,:));
  for i = 1:size(p,2);
    pp = p(:,i);
    plot(tt(~isnan(pp)),pp(~isnan(pp)),'linewidth',1);
  end
  
  set(a(end),'fontsize',8);
  
  tstmstart = tt(find([diff(info==0);0]<0,1,'first'));
  tstmbreakend = tt(find([diff(info==3);0]<0));
  tstmbreakstart = tt(find([diff(info==3);0]>0));
  tspan = tstmbreakend(1) - tstmbreakstart(1);
  xlim([tstmstart-tspan,tstmbreakend(2)-tspan/2])
  ylim([1,xend]);
  ylabel('X-Position [px]','fontsize',10);    
  xlabel('Time [min]','fontsize',10);    
  title('Example: Group size 5 fish','fontsize',10,'fontweight','normal')
  
  %population
  s = s+2;
  a(end+1) = subplot(r1,r2+1,s,'align');
  
  Nl = [];
  Nr = [];
  nbins = 50;

  pos = ftstmall{5}{1}.ft.res.pos; % estimate from 5 fish
  x = pos(:,1,:);
  edges = linspace(quantile(x(:),0.001),quantile(x(:),0.999),nbins);

  for nfish = 1:length(ftstmall)
    for i = 1:length(ftstmall{nfish})
      pos = ftstmall{nfish}{i}.ft.res.pos;
      res = ftstmall{nfish}{i}.ft.res.tracks;

      if size(pos,3)<=2 % SOMEHOW FOR VIDEO WITH 1 or 2 fish stimulus
                        % was reversed!!
        left = 2;
        right = 1;
      else
        left = 1;
        right = 2;
      end
  
      c = res.consequtiveInvisibleCount;
      cc = permute(cat(3,c,c),[1,3,2]);
      pos(cc>0) = NaN;

      tt = res.t(:,1);
      idx = res.stmInfo(:,1,2) == left ;
      
      Nl{nfish}(:,:,i) = histc(squeeze(pos(idx,1,:)),edges)/sum(idx);
      Nl{nfish}(:,:,i) = bsxfun(@rdivide,Nl{nfish}(:,:,i),sum(Nl{nfish}(:,:,i),1));

      idx = res.stmInfo(:,1,2) == right ;
      Nr{nfish}(:,:,i) = histc(squeeze(pos(idx,1,:)),edges)/sum(idx);
      Nr{nfish}(:,:,i) = bsxfun(@rdivide,Nr{nfish}(:,:,i),sum(Nr{nfish}(:,:,i),1));
    end
  end

  mNr = [];
  mNl = [];
  sNr = [];
  sNl = [];
  for i = 1:length(Nr)
    if isempty(Nr{i})
      mNr(:,i) = NaN;
      sNr(:,i) = NaN;
      mNl(:,i) = NaN;
      sNl(:,i) = NaN;
    else
      mNr(:,i) = nanmean(Nr{i}(:,:),2);
      sNr(:,i) = stderr(Nr{i}(:,:),2);
      mNl(:,i) = nanmean(Nl{i}(:,:),2);
      sNl(:,i) = stderr(Nl{i}(:,:),2);
    end
  end
  tstart = 0;
  tend = 0.5; 
  xend = 0;
  xhalf = 1;
  patch([tstart,tend,tend,tstart,tstart],[xend,xend,xhalf,xhalf,xend],gray,'edgecolor','none');
  hold on;
  h = errorbarpatch(linspace(0,1,length(edges)),(mNr + mNl(end:-1:1,:))/2,(sNr + sNl(end:-1:1,:))/2,'linewidth',1);
  legend(h([1,2,3,5]),{'1 fish','2 fish','3 fish','5 fish'});
  set(a(end),'fontsize',8);
  xlabel(sprintf('Relative location\n during stimulation'),'fontsize',10)
  ylabel('Probability','fontsize',10)
  box off;
  ylim([0,0.25])

  s = s+1;
  a(end+1) = subplot(r1,r2+1,s,'align');
  
  y = nan(1,length(Nl));
  sy = nan(1,length(Nl));
  x = 1:length(Nl);
  for i = 1:length(Nl)
   if isempty(Nl{i}) 
     continue;
   end
   m =[sum(Nr{i}(floor(nbins/2):end,:),1),sum(Nl{i}(1:floor(nbins/2),:),1)];
   y(i) = mean(m);
   sy(i) = stderr(m);
  end
  
  errorbarline(x,y*100,sy*100,'x');
  set(a(end),'fontsize',8);
  ylabel(sprintf('Fraction of time on\n dark side [%%]'), 'fontsize',10);
  xlabel('Group size [# fish]','fontsize',10);
  
  box off;



  b = labelsubplot(gcf);
  shiftaxes(b(1),[0.06])
  shiftaxes(b(2:3),[-0.01,-0.01])

  if SAVEIF
    exportfig(gcf,'/home/malte/work/projects/zebra/halfplane.eps',...
              'color','rgb','FontSizeMin',5,'FontMode','scaled','LineWidthMin',0.75)
    saveas(gcf,'/home/malte/work/projects/zebra/halfplane.fig');
  end



end

