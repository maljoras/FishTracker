

LOAD = 0;
COMPUTE = 0;
PLOT = 1;
SAVEIF = 1;

VELTRIGGERPLOT = 0;

if LOAD 
 load ~/data/zebra/videos/textures/texturesinregionsF4-1-1.mat
 T.getTrackingResults([],[],1); % force to generate new DAG
end


if COMPUTE

  timeRange = [T.stimulusPresenter.adaptationTime, 10*3600];  
  dagif = 1;
  
  b = T.stimulusPresenter.borderWidth;

  res = T.getTrackingResults(timeRange,dagif);
  turn = T.getTurningStats(res);

  % we first have to match the current fishID to the stimIdentityID at the time of
  % stimulus presentations
  stmInfo = res.tracks.stmInfo; % this is sorted accroding to gloval estimate
                                % of fishIDs (might include some ID switches)
  stmIdentityId = stmInfo(:,:,T.stimulusPresenter.IDX_IDENTITYID);

  pos = T.interpolateInvisible(res,'pos',3);
  velocity = [zeros(1,2,T.nbody);bsxfun(@rdivide,diff(pos,1),diff(res.tabs))];
  avelocity = squeeze(sqrt(sum(velocity.^2,2)));
  maxvel = quantile(avelocity(:),0.99);
  dlmsk = avelocity>maxvel;
  avelocity(dlmsk) = NaN;
  pos(dlmsk) = NaN;
  

  % we need to raw data to reconstruct the velocity at stimulus
  % onset
  rawres = T.getRawTrackData();
  
  velfishid = sqrt(sum(T.getResField(res,'velocity').^2,3));
  trackid  = T.getResField(res,'id');

  rawfishid  = rawres.stmInfo(:,:,T.stimulusPresenter.IDX_IDENTITYID);
  velfishraw = sqrt(sum(rawres.velocity.^2,3));
  
  % msk includes boundary off (same as ~isnan(bbox(1)))'
  stmmsk = stmInfo(:,:,T.stimulusPresenter.IDX_MSK); 
  stmbbox = zeros([size(stmmsk),length(T.stimulusPresenter.IDX_BBOX)]);
  for i = 1:size(stmInfo,2)
    stmbbox(:,i,:) = T.stimulusPresenter.fromScreenBbox( squeeze(stmInfo(:,i,T.stimulusPresenter.IDX_BBOX))); 
  end
  
  stmcenter = stmbbox(:,:,1:2) + stmbbox(:,:,3:4)/2;
  stmxy = stmInfo(:,:,T.stimulusPresenter.IDX_XY); 
  stmdx = stmcenter - permute(pos,[1,3,2]);

  stmsf = stmInfo(:,:,T.stimulusPresenter.IDX_SIZEFACTOR);
  stmshift = stmInfo(:,:,T.stimulusPresenter.IDX_SHIFT);
  stmshiftori = stmInfo(:,:,T.stimulusPresenter.IDX_SHIFTORI);
  theory_stmdx = cat(3,sin(stmshiftori).*stmshift,cos(stmshiftori).*stmshift);

  avelregion = []; 
  xreg = T.stimulusPresenter.xRegions;
  for j = 1:length(xreg)
    msk = stmxy(:,:,1)>b & stmxy(:,:,1)<=xreg(j);
    tmp = avelocity;
    tmp(~msk) = NaN;
    avelregion(:,j) = nanmedian(tmp);
  end
  tmp = avelocity;
  msk = stmxy(:,:,1)<1-b & stmxy(:,:,1)>xreg(end);
  tmp(msk) = NaN;
  avelregion(:,j+1) = nanmedian(tmp);
  avelregion =  avelregion/T.bodylength
  
  stmstate = stmInfo(:,:,T.stimulusPresenter.IDX_STATE); 
  tt  = 1:size(stmstate,1);
  msk = isnan(res.tabs);
  stmstate(msk,:) = interp1(tt(~msk),stmstate(~msk,:),tt(msk),'nearest');
  stmstate(isnan(stmstate)) = 0;
  
  %sf = stmInfo(:,:,T.stimulusPresenter.IDX_SIZEFACTOR);
  %stmstate = stmstate.*(sf>0); % stmmsk sets sf to zero

  % NOTE: some of the stmstates correspond to sf==0. These can be used for
  % control comparison (same length statistics)
  
  r1 = 3;
  r2 = 3;
  
  % do not use stmmsk because occasonially has boundary effect 
  
  for i = 1:length(turn)
    stmstatei = stmstate(:,i);
    
    stmtidx_start = find(diff([0;stmstatei])==1);
    stmtidx_stop = find(diff([stmstatei;0])==-1);

    delmsk =  isnan(res.iframe(stmtidx_start)) | isnan(res.iframe(stmtidx_stop));
    stmtidx_start(delmsk) = [];
    stmtidx_stop(delmsk) = [];
    
    turn(i).stm.tidx = stmtidx_start;

    nstm = length(stmtidx_start);
    turn(i).stm.stoptidx = stmtidx_stop;
    turn(i).stm.tlen =res.t(stmtidx_stop)- res.t(stmtidx_start);
    
    
    
    accummsk = zeros(size(stmstate,1),1);
    accummsk(stmtidx_start) = 1;
    accummsk = cumsum(accummsk);
    accummsk(~accummsk) = nstm+1;

    accummsk_nogap = accummsk;
    turn(i).stm.accummsk_nogap = accummsk_nogap;

    accummsk_gap  =  accummsk;
    accummsk_gap(logical(stmstatei)) = nstm+1;
    turn(i).stm.accummsk_gap = accummsk_gap;
    
    accummsk(~stmstatei) = nstm+1;
    turn(i).stm.accummsk = accummsk;

    for f = {'IDX_SIZEFACTOR','IDX_SHIFTORI','IDX_SHIFT','IDX_XY','IDX_IDENTITYID'}
      turn(i).stm.(lower(f{1}(5:end))) = stmInfo(stmtidx_start,i,T.stimulusPresenter.(f{1}));
    end
    turn(i).stm.acc_identityId = accumarray(accummsk,stmIdentityId(:,i)==i,[nstm+1,1],@nanmean);
    turn(i).stm.acc_len = accumarray(accummsk,1,[nstm+1,1],@sum);

    
    ind = res.iframe(turn(i).stm.tidx);
    v = velfishraw(ind,:)';
    fidx = find(bsxfun(@eq,rawfishid(ind,:),turn(i).stm.fishid)');
    turn(i).stm.velid = v(fidx);

    
    turnmsk = zeros(size(accummsk));
    turnmsk(turn(i).tidx) = 1;
    turn(i).stm.acc_nturns = accumarray(accummsk,turnmsk,[nstm+1,1],@sum);

    turnmsk(turn(i).tidx) = 1:length(turn(i).tidx);
    turnmsk(~turnmsk) = NaN;
    turn(i).stm.acc_firstturn_idx = accumarray(accummsk,turnmsk,[nstm+1,1],@nanmin,NaN);
    turn(i).stm.acc_lastturn_idx = accumarray(accummsk,turnmsk,[nstm+1,1],@nanmax,NaN);
    turn(i).stm.acc_gapfirstturn_idx = accumarray(accummsk_gap,turnmsk,[nstm+1,1],@nanmin,NaN);
    turn(i).stm.acc_gaplastturn_idx = accumarray(accummsk_gap,turnmsk,[nstm+1,1],@nanmax,NaN);

    turn(i).stm.acc_nogapfirstturn_idx = accumarray(accummsk_nogap,turnmsk,[nstm+1,1],@nanmin,NaN);
    turn(i).stm.acc_nogaplastturn_idx = accumarray(accummsk_nogap,turnmsk,[nstm+1,1],@nanmax,NaN);
  
  end
  

end




if PLOT

  figure;
  
  %% turning
  clf;

  a = subplot(1,1,1);
  
  ifish = 2; % 2 % 2
  iidx = 20; % 15 % 20
  toffs = 22; %index 
  
  idx = find(turn(ifish).stm.sizefactor>1 & turn(ifish).stm.sizefactor<1.3);

  istart = -toffs + turn(ifish).stm.tidx(idx(iidx));
  tstart = res.t(istart);
  istop = turn(ifish).stm.stoptidx(idx(iidx))+toffs;
  tend = res.t(istop);

  trange = [tstart tend];
  rest = T.getTrackingResults(trange,dagif);
  centerLine = T.getResField(rest,'centerLine',0);
  clx = permute(centerLine(:,:,1,:),[4,1,2,3]);
  cly = permute(centerLine(:,:,2,:),[4,1,2,3]);
  cic = rest.tracks.consecutiveInvisibleCount;
  mclx = T.interpolateInvisible(shiftdim(mean(clx,1),1),cic);
  mcly =  T.interpolateInvisible(shiftdim(mean(cly,1),1),cic);
  turnt = T.getTurningPoints(rest);


  
  col = jet(T.nbody);
  for i = 1:T.nbody
    g  =[0.5,0.5,0.5];
    plot(clx(:,:,i),cly(:,:,i),'color',col(i,:)*0.8);
    hold on;
    mx = mclx(:,i);
    my = mcly(:,i);
    x = turnt(i).x;
    y = turnt(i).y;
    plot(mx,my,'color',col(i,:),'linewidth',2);
    plot(x,y,'o','linewidth',2,'markersize',4,'color',col(i,:))
    plot(x,y,'.','linewidth',2,'markersize',10,'color','k')
    %plot(mx(idx_thr),my(idx_thr),'xb','linewidth',2)
      
    L = 50;
    xx = [x, x + L*cos(turnt(i).ori(:))]';
    yy = [y, y + L*sin(turnt(i).ori(:))]';
    plot(xx,yy,'-','linewidth',2,'color',col(i,:))
    plot(xx,yy,'-','linewidth',1,'color','k')
    
    bbox = squeeze(stmbbox(istart:istop,i,:));
    ss = 0;
    m = sum(~isnan(bbox(:,1)));
    for ii = 1:size(bbox,1)
      
      if ~isnan(bbox(ii,1))
        ss = ss+1;
        rectangle('position',bbox(ii,:),'Edgecolor',(1-ss/m)*ones(1,3));
      end
    end
  end

  fs = T.videoHandler.frameSize;
  xlim([1,fs(2)])
  ylim([1,fs(1)]);
  
  rectangle('position',T.stimulusPresenter.screenBoundingBox,'edgecolor','r')
  sbbox = T.stimulusPresenter.screenBoundingBox;

  bw = T.stimulusPresenter.borderWidth;
  border = T.stimulusPresenter.fromScreenBbox(T.stimulusPresenter.toScreenRect([bw,bw,1-bw*2,1-bw*2]));
  rectangle('position',border,'edgecolor','r','linestyle',':')
  
  xreg = T.stimulusPresenter.xRegions;
  l1 = T.stimulusPresenter.fromScreenBbox(T.stimulusPresenter.toScreenRect([xreg(1),bw,1-bw*2,1-bw*2]));
  plot([l1(1),l1(1)],[l1(2),l1(2)+l1(4)],'r--')
  l2 = T.stimulusPresenter.fromScreenBbox(T.stimulusPresenter.toScreenRect([xreg(2),bw,1-bw*2,1-bw*2]));
  plot([l2(1),l2(1)],[l2(2),l2(2)+l2(4)],'r--')
  xlabel('X-position [px]')
  ylabel('Y-position [px]')
  

  regscale = T.stimulusPresenter.regSizeFactorScale;
  g = [];amount = 10000;
  for ii = 1:length(regscale)
    gmu = T.stimulusPresenter.stmSizeFactor*ones(amount,1)*regscale(ii);
    gCV = T.stimulusPresenter.stmSizeFactorCV;
    g(:,ii) = T.stimulusPresenter.grand(gmu,gCV);
  end

  box on;

  col = get(gca,'colororder');
  for i = 1:3
    text((l1(1)-sbbox(3)/6) + (i-1)*sbbox(3)/3, 1.4*l1(2),sprintf('Region %d',i),'color',col(i,:),'horiz','center');
  end
  xoffset = -.16;
  shiftaxes(a,[-0.03,0,xoffset]);
  
  
  a = xy.helper.inset(gca,'switchsides',1,'margin',0.15,'left',0.75,'lower',0.75);
  set(a,'visible','on')
  [N,edges] = hist(g,50);
  h = plot(edges,N/amount,'.','linewidth',1);
  ylim([0,0.2])
  xlim(a,[0,5])
  xlabel('Size factor [BL]');
  ylabel('Prob.');
  box on;
  set(a,'fontsize',10)
  %shiftaxes(a,[0,-0.35])
  
  

  %% plot velocity probmap
  szFrame = T.videoHandler.frameSize;
  for i = 1:T.nbody
    aa = xy.helper.subsubplot(1,4,4,4,1,i);
    P1 = T.plotVelocityMap(timeRange,i,[0,0.5]);
    P2 = T.plotVelocityMap(timeRange,i,[0.5,0.99]);
    
    Z = imfilter(P2-P1,fspecial('Gaussian',[5,5],1),'same');
    imagesc(1:szFrame(2),1:szFrame(1),Z,[-0.6,0.6]);
    hold on;
    fs = T.videoHandler.frameSize;
    xlim([1,fs(2)])
    ylim([1,fs(1)]);
    
    rectangle('position',T.stimulusPresenter.screenBoundingBox,'edgecolor','r')
    sbbox = T.stimulusPresenter.screenBoundingBox;
    
    bw = T.stimulusPresenter.borderWidth;
    border = T.stimulusPresenter.fromScreenBbox(T.stimulusPresenter.toScreenRect([bw,bw,1-bw*2,1-bw*2]));
    rectangle('position',border,'edgecolor','r','linestyle',':')
    
    xreg = T.stimulusPresenter.xRegions;
    l1 = T.stimulusPresenter.fromScreenBbox(T.stimulusPresenter.toScreenRect([xreg(1),bw,1-bw*2,1-bw*2]));
    plot([l1(1),l1(1)],[l1(2),l1(2)+l1(4)],'r--')
    l2 = T.stimulusPresenter.fromScreenBbox(T.stimulusPresenter.toScreenRect([xreg(2),bw,1-bw*2,1-bw*2]));
    plot([l2(1),l2(1)],[l2(2),l2(2)+l2(4)],'r--')
    
    ylabel('');
    xlabel('');
    set(aa,'xticklabel',[]);
    shiftaxes(aa,[0.04,0,-0.02]);
    axis xy;
    title(sprintf('#%d',i),'color',col(i,:));
    set(aa,'fontsize',8);
    if i>1
      set(aa,'tag','nolabel');
    end

  end
  
  xlabel('X-position [px]');
  set(aa,'xticklabelmode','auto');
  h = xy.helper.smallcolorbar();
  set(h,'fontsize',8)
  ylabel(h,'Diff. probability','fontsize',10)


  if SAVEIF

    a = xy.helper.labelsubplot(gcf,'BA');
    shiftaxes(a(2),0.05)
    %set(gcf,'paperunits','centimeters')
    exportfig(gcf,'~/work/writings/papers/fishtracking/figs/turning1.eps','Color','rgb','width',10,'fontmode','fixed');
  end
  
  
  
  
  if VELTRIGGERPLOT
    figure;
    r1 = 2;
    r2 = T.nbody;
    for i = 1:length(turn)

      
      %% stim average turning velocity with position
      subplot(r1,r2,i);cla;
      
      tuidx = turn(i).stm.acc_firstturn_idx;
      delmsk = isnan(tuidx);
      tuidx(delmsk) = [];
      delmsk = delmsk(1:end-1);
      dori = turn(i).dori(tuidx(1:end-1));
      vel0 = turn(i).locvel0(tuidx(1:end-1));
      len = turn(i).len(tuidx(1:end-1));
      ind = find(~delmsk);
      sf = turn(i).stm.sizefactor(ind);
      xy = turn(i).stm.xy(ind,:);
      shift = turn(i).stm.shift(ind,:);
      shiftori = turn(i).stm.shiftori(ind,:);

      stmstatei = stmstate(:,i);
      
      nbins = 20;
      clu = 1.5;
      cl = [-clu-0.1,clu];
      edges = linspace(-60,60,nbins);
      stmdxi = [shift.*sin(shiftori),shift.*cos(shiftori)];
      
      [N,ibin] = histc(stmdxi,edges);
      ibin(~ibin) = 1;
      
      ibin_sf0 = ibin;
      ibin_sf0(sf~=0,:) = 1;
      ibin_sf1 = ibin;
      ibin_sf1( sf==0,:) = 1;
      

      len(len>500) = NaN;
      z = vel0/T.bodylength;

      z( xy(:,1)<b | xy(:,1)>1-b) = NaN;
      
      imori0 = accumarray(ibin_sf0,z,[nbins,nbins],@nanmedian,NaN);
      imori1 = accumarray(ibin_sf1,z,[nbins,nbins],@nanmedian,NaN);
      imori0(1,1) = NaN;
      
      Z = imori1-imori0;
      dmsk = isnan(Z); 
      Z = max(imori1-imori0,-clu);
      Z(dmsk) = -clu-0.1;

      
      x = edges + diff(edges(1:2))/2;
      imagesc(x,x,Z,cl);
      C = parula;
      C(1,:) = [1,1,1];
      colormap(C);
      
      hold on;
      fl = T.bodylength*0.5;
      fw = T.bodywidth*0.5;
      rectangle('position',[-fw/2,-fl/2,fw,fl],...
                'curvature',1,'facecolor',[1,1,1]*.5,'edgecolor','none')
      
      xlabel('X rel. shift [px]','fontsize',10);
      if i==1
        ylabel('Y  rel. shift [px]','fontsize',10);
      end
      
      col = jet(T.nbody);
      title(sprintf(' #%d',i),'color',col(i,:));

      if i==T.nbody
        h = xy.helper.smallcolorbar;
        xy.helper.shiftaxes(h,[0.1]);
        ylabel(h,'Difference speed [BL/sec]')
        %title(sprintf('Stimulus induced speed\n at first turning'))
        set(h,'ylim',[-clu,clu]);
      end
      
      
      
      %% velocity ETA
      subplot(r1,r2,r2+i);cla;

      sfthres = 1;
      vel = avelocity(:,i)/T.bodylength;

      % interpolate
      dt = 1/T.videoHandler.frameRate;
      msk = ~isnan(vel);
      %ivel = interp1(res.tabs(msk),vel(msk),(res.tabs(1):dt:res.tabs(end))','next');
      %istmx = interp1(res.tabs(msk),stmxy(msk,i,1),(res.tabs(1):dt:res.tabs(end))','next');
      
      sf = turn(i).stm.sizefactor;
      mxt = 40;
      tidx =   turn(i).stm.tidx;
      tidx_stop =   turn(i).stm.stoptidx;
      xy = turn(i).stm.xy(:,:);
      velid = turn(i).stm.velid;
      toffs = -10;
      mvel = [];svel = [];
      tlen = turn(i).stm.tlen;
      s = 0;
      b = T.stimulusPresenter.borderWidth;
      for indmsk = {sf==0,sf>0 & sf<sfthres, sf>=sfthres}
        s = s+1;
        
        msk = indmsk{1} & ~(xy(:,1)<b | xy(:,1)>1-b | velid< ...
                            T.stimulusPresenter.stmVelThres| tlen<0.5);
        
        tidx1 = tidx(msk)+toffs;
        tidx1_stop = tidx_stop(msk)+toffs;

        % convert to interp1 idx
        %itidx1 = floor((res.t(tidx1)-res.t(1))/dt)+1;
        %itidx1_stop = floor((res.t(tidx1_stop)-res.t(1))/dt)+1;

        accmsk = ones(size(vel,1),1);      
        accmsk(1:tidx1(1)-1) = 0;
        accmsk(tidx1_stop(1:end-1)+1) = -(tidx1(2:end)-tidx1(1:end-1))+1;
        accmsk = cumsum(accmsk);
        
        accmsk(accmsk<=0 | accmsk>mxt ) = mxt+1;
        
        % avg turn-triggered velocity
        mvel(:,s) = accumarray(accmsk,vel,[mxt+1,1],@nanmedian);
        svel(:,s) = accumarray(accmsk,vel,[mxt+1,1],@xy.helper.stderr);
        
        %mposx(:,s) = accumarray(accmsk,istmx,[mxt+1,1],@nanmedian);
      end
      
      %mvel = bsxfun(@minus,mvel,avelregion(i,:));
% $$$     xreg = [0,T.stimulusPresenter.xRegions,1];
% $$$     for j = 1:length(xreg)-1
% $$$       msk1 = (accmsk==mxt+1 & istmx>xreg(j) & istmx<=xreg(j+1) & istmx>b & istmx<1-b);
% $$$       mvel(:,j) = mvel(:,j) - nanmedian(ivel(msk1));
% $$$     end
      
      h = xy.helper.errorbarpatch(((0:mxt-1)+toffs)*dt,mvel(1:end-1,:),svel(1:end-1,:));
      xlabel('Rel time [s]','fontsize',10)
      if i==1
        ylabel('Avg. velocity [BL/s]','fontsize',10);
      end
      %xlim([0,1.5])
      %ylim([0,5])
      box off;
      if i==2
        legend('size factor = 0',sprintf('size factor < %1.1f', sfthres),sprintf(['size factor ' '> %1.1f'],sfthres))
      end
      
    end

    if SAVEIF
      exportfig(gcf,'~/work/projects/zebra/turning1.pdf','Color','rgb')
    end
  end
  
end
