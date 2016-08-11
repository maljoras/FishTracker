LOAD = 0;
PLOT = 1;
SAVEIF = 1;
COMPUTE = 0;

NEWTRACK = 0
LOADALL = 0;

if LOAD || ~exist('xystm','var')

  dname = '~/data/zebra/videos/halfplane';
  fname = 'five2onefish2';
  

  vid = [dname filesep fname '.avi'];
  xystm = load([dname filesep fname '.mat']);

  stmpos = xystm.T.res.pos;
  
end

if LOADALL
  dname = '~/data/zebra/videos/halfplane';
  fnames = {{'onefish1','onefish2','onefish3'},{'twofish1','twofish2','twofish3'},{'testThreefish','testThreefish2','testThreefish3','testThreefish4','testThreefish5'},{},{'fivefish1','fivefish2','fivefish3'}};
  
  for i = 1:length(fnames)

    for j = 1:length(fnames{i})
      xystmall{i}{j} = load([dname filesep fnames{i}{j} '.mat']);
    end
  end
    
end

  


if COMPUTE && NEWTRACK % for tracking again
  
  T = xy.Tracker(vid,'detector.adjustThresScale',1,'nindiv',2,'detector.inverted',1);

  T.addSaveFields('firstFrameOfCrossing', 'lastFrameOfCrossing');
  T.setDisplay(0);
  T.setDisplay('switchIdentity',0,'tracks',0,'level', 3);
  
  tic;
  T.track();
  toc;
end

if NEWTRACK && ~exist('stmres1','var')
  %% get the adjusted stminfo 
  if ~T.res.t(1,1)
    tfile = dlmread([vid '.txt'])  ;
    t = tfile(1,3) + T.res.t(:,1);
  else
    t =  T.res.t(:,1);
  end
  tc = ceil(t*1e4)/1e4;
  tf = floor(t*1e4)/1e4;
  
  stmres = xystm.T.res.tracks;
  xyres = T.res.tracks;
  xypos = T.res.pos;


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
  for i = 1:T.nindiv
    for j = 1:T.nindiv
      dist(i,j) = nanmean(sqrt(sum((xypos(:,:,i) - stmpos1(:,:,j)).^2,2)));
    end
  end

  assignments = assignDetectionsToTracks(dist,1e3);
  change = assignments(assignments(:,1),2);
  xypos  = xypos(:,:,change);
  for f = fieldnames(xyres)'
    xyres.(f{1}) = xyres.(f{1})(:,change,:,:,:);
  end
  xyres.stmInfo = stmres1.stmInfo;

  % check tracking differences
  %figure;
  %plot(t,squeeze(stmpos1(:,1,:)-xypos(:,1,:)));
  
end


if PLOT
  figure;
  r1  = 2;
  r2 = 2;
  s = 1;
  a = [];
  n_plot = 2;

  s = s+1;
  for i_plot =1:n_plot
    
    aa = subsubplot(r1,r2,s,n_plot,1,i_plot);

    if i_plot==2
      if NEWTRACK
        res = xyres; % take new
        pos = xypos;
        tt = t;
      else
        %res = xystm.T.res.tracks;
        %pos = stmpos;
        %tt = xystm.T.res.t(:,1);
        xystm2 = xystmall{5}{1};
        xystm2.T.videoFile
        pos = xystm2.T.res.pos;
        tt = xystm2.T.res.t(:,1);    
        res = xystm2.T.res.tracks;
        tstr = 'Group size 5 fish';

      end
    elseif i_plot==1
        xystm2 = xystmall{2}{1};
        tstr = 'Group size 2 fish';
        xystm2.T.videoFile
        pos = xystm2.T.res.pos;
        tt = xystm2.T.res.t(:,1);    
        res = xystm2.T.res.tracks;
    else
      error('Do not know what to plot');
    end

        
    c = res.consecutiveInvisibleCount;
    cc = permute(cat(3,c,c),[1,3,2]);
    pos(cc>0) = NaN;
    

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

    xhalf = xystm.T.videoHandler.frameSize(2)/2;
    xend = xystm.T.videoHandler.frameSize(2);
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
    set(aa,'colororder',jet(size(res.t,2)))
    p = squeeze(pos(:,1,:));
    nconv = 24;
    for i = 1:size(p,2);
      pp = p(:,i);
      ipp = interp1(tt(~isnan(pp)),pp(~isnan(pp)),tt(isnan(pp)),'linear','extrap');
      pp(isnan(pp)) = ipp;
      cpp = conv(pp,ones(nconv,1)/nconv,'same');
      plot(tt,cpp,'linewidth',0.5);
    end
    
    set(aa,'fontsize',8);
    
    tstmstart = tt(find([diff(info==0);0]<0,1,'first'));
    tstmbreakend = tt(find([diff(info==3);0]<0));
    tstmbreakstart = tt(find([diff(info==3);0]>0));
    tspan = tstmbreakend(1) - tstmbreakstart(1);
    xlim([tstmstart-tspan,tstmbreakend(2)-tspan/2])
    ylim([1,xend]);

    if i_plot==n_plot
      xlabel('Time [min]','fontsize',10);    
    else
      set(aa,'xticklabel',[]);
    end
    ttl = title(tstr,'fontsize',8,'fontweight','normal','vert','middle');
    p1 = get(ttl,'position');
    p1(2) = p1(2)+80;
    set(ttl,'position',p1);
    axislabel off;
  end
  a(end+1) = subsubplot(r1,r2,s,1,1,1);
  ylim([1,xend]);
  xlim([tstmstart-tspan,tstmbreakend(2)-tspan/2])
  yl = ylabel('X-Position [px]','fontsize',10);    
  p1 = get(yl,'position');
  set(a(end),'visible','off');
  text(p1(1),p1(2),'X-Position [px]','fontsize',10,'rotation',90,'vert','bottom','horiz','center');
  
  %population
  s = s+1;
  a(end+1) = subplot(r1,r2,s,'align');
  
  Nl = [];
  Nr = [];
  nbins = 50;

  pos = xystmall{5}{1}.T.res.pos; % estimate from 5 fish
  x = pos(:,1,:);
  edges = linspace(quantile(x(:),0.001),quantile(x(:),0.999),nbins);

  for nindiv = 1:length(xystmall)
    for i = 1:length(xystmall{nindiv})
      pos = xystmall{nindiv}{i}.T.res.pos;
      res = xystmall{nindiv}{i}.T.res.tracks;

      if size(pos,3)<=2 % SOMEHOW FOR VIDEO WITH 1 or 2 fish stimulus
                        % was reversed!!
        left = 2;
        right = 1;
      else
        left = 1;
        right = 2;
      end
  
      c = res.consecutiveInvisibleCount;
      cc = permute(cat(3,c,c),[1,3,2]);
      pos(cc>0) = NaN;

      tt = res.t(:,1);
      idx = res.stmInfo(:,1,2) == left ;
      
      Nl{nindiv}(:,:,i) = histc(squeeze(pos(idx,1,:)),edges)/sum(idx);
      Nl{nindiv}(:,:,i) = bsxfun(@rdivide,Nl{nindiv}(:,:,i),sum(Nl{nindiv}(:,:,i),1));

      idx = res.stmInfo(:,1,2) == right ;
      Nr{nindiv}(:,:,i) = histc(squeeze(pos(idx,1,:)),edges)/sum(idx);
      Nr{nindiv}(:,:,i) = bsxfun(@rdivide,Nr{nindiv}(:,:,i),sum(Nr{nindiv}(:,:,i),1));
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
  a(end+1) = subplot(r1,r2,s,'align');
  
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
  xlim([0.5,length(y)+0.5])
  box off;



  b = labelsubplot(gcf);
  shiftaxes(b(1),[0.06])
  shiftaxes(b(2:3),[-0.01,-0.01])

  if SAVEIF
    exportfig(gcf,'/home/malte/work/writings/papers/fishtracking/figs/halfplane.eps','color','rgb','FontSizeMin',5,'FontMode','scaled')
    
    saveas(gcf,'/home/malte/work/writings/papers/fishtracking/figs/halfplane.fig');
  end



end

