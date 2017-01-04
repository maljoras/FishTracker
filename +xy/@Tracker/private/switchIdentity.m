function switchIdentity(self,trackIndices,assignedIdentityIds,crossingflag)
% changes the IdentityIds of the tracks. expects a premutations. updates all relevant identityId
% dependet varibeles. Resets the data to avoid the learning of wrong features. Also
% resets the uniqueIdentityFrame counter
%
% if the crossing flag is set, the switching point is calculated based on the
% distance within the last crossing. If not set, the switchingpoint is calculated
% based on the classProbHistory. 

  
  minDagFrames =  25;
  

  assert(length(trackIndices)==length(assignedIdentityIds));
  oldIdentityIds = [self.tracks(trackIndices).identityId];

  if length(trackIndices)==1 || all(oldIdentityIds==assignedIdentityIds)
    return; % no need to do anything
  end

  [memb,oldloc] = ismember(assignedIdentityIds,oldIdentityIds);
  assert(all(memb));  
  

  f2t = self.identityId2TrackId;
  orgpos = self.pos;
  self.resetBatchIdx(trackIndices); % reset for all
  self.uniqueIdentityFrames = 0;
  velocity = cat(1,self.tracks(trackIndices).velocity);

  if self.opts.display.switchIdentity && self.displayif && self.opts.display.tracks;
    self.displayCurrentTracks();
  end
  
  
  % start end time definition
  tfirst =[self.tracks(trackIndices).firstFrameOfCrossing]; 
  tlast = [self.tracks(trackIndices).lastFrameOfCrossing]; 
  tcurrent = self.currentFrame;
      
  if crossingflag
    tstart = max(min(tfirst) - 2*self.nFramesAfterCrossing,1);
    tend = min(self.currentFrame,max(tlast) + self.nFramesAfterCrossing);
  else
    tstart = max(max(tlast) - 2*self.nFramesAfterCrossing,1) ; % put about last crossing in...
    tend = tcurrent;
  end

  
  idxOldIdentityIds = [self.tracks(trackIndices).identityId]; 
  idxTrackIndices = trackIndices;
  idxAssignedIdentityIds = assignedIdentityIds;
  [tminfinallost,delmskall] = subCalcLostTrackBasedSwitchPoint();

  %% get the dag pos for the whole crossing
  [dagpos,dagf2t,errorflag]= subGetDagBasedPos();
  errorflag = errorflag | (length(assignedIdentityIds)~=size(dagf2t,2));  % could this happen?
  
  if ~errorflag && (crossingflag  || length(trackIndices)>2|| (self.currentFrame-tstart > minDagFrames))
    % use dag pos
    if self.verbosity>1
      xy.helper.verbose('Use DAG based switching');
    end
    
    t = self.currentFrame-size(dagf2t,1)+1:self.currentFrame;   
    self.identityId2TrackId(t,assignedIdentityIds) = dagf2t;
    self.pos(:,assignedIdentityIds,t) = dagpos;
    

    steps = -length(t)+1;
    for j = 1:length(trackIndices)
      self.tracks(trackIndices(j)).identityId = assignedIdentityIds(j);
      self.tracks(trackIndices(j)).switchedFrames = steps;

      % switch other stuff... seems not working with velocity.
      self.tracks(trackIndices(j)).velocity = velocity(oldloc(j),:);
    end
    
  else
  
    %% switch the tracks for each connected component  
    if self.verbosity>1
      xy.helper.verbose('Use conventional switching');
    end
    
    [~, ncidx]= self.connectedComponents(trackIndices,assignedIdentityIds);
    for i_ncidx = 1:length(ncidx)  % actually do not need to check again.
      idx = ncidx{i_ncidx};
      
      if length(idx)==1
        continue
      end
      
      % define common pars for subfunctions
      idxOldIdentityIds = [self.tracks(trackIndices(idx)).identityId]; 
      idxTrackIndices = trackIndices(idx);
      idxAssignedIdentityIds = assignedIdentityIds(idx);

      % calculate switch points & delmsk
      [tminfinallost,delmsk] = subCalcLostTrackBasedSwitchPoint();
      tminfinaldist = subCalcDistanceBasedSwitchPoint();
      
      
      
      %% choose one of the points..
      if crossingflag
        if ~isempty(tminfinallost)
          tminfinal = tminfinallost;
        else
          tminfinal = tminfinaldist;
        end
      else 
        % identityClassUpdate
        tminprob = subCalcClassProbBasedSwitchPoint();
        
        if tcurrent-tminprob> 2*self.nFramesAfterCrossing
          tminfinal = tminprob;
        else
          if length(idx)>2
            tminfinal = subCalcDistanceBasedSwitchPoint();
          else
            tminfinal = tminfinallost;
          end
        end
      end
      
      
      %% rather swap using the tminfinal
      t = tminfinal:self.currentFrame;    
      steps = tminfinal-self.currentFrame;
      for j = idx(:)'
        self.identityId2TrackId(t,assignedIdentityIds(j)) = f2t(t,oldIdentityIds(j));
        self.pos(:,assignedIdentityIds(j),t) = orgpos(:,oldIdentityIds(j),t);
        self.tracks(trackIndices(j)).identityId = assignedIdentityIds(j);
        self.tracks(trackIndices(j)).switchedFrames = steps;

        % switch other stuff
        self.tracks(trackIndices(j)).velocity = velocity(oldloc(j),:);
      end


    end % for ncidx

    
    if crossingflag 
      % delete some of the critical points
      tcross = tstart:tend;
      p = self.pos(:,assignedIdentityIds,tcross);
      p(:,delmskall') = NaN;
      self.pos(:,assignedIdentityIds,tcross) = p;
    end

    if self.opts.display.switchIdentity && self.displayif
      subPlotCrossings();
    end

  end
  
    
    

  function subPlotCrossings()

    %define just for the plotting  (might be confused by the loop
    %above. Just ignore)
    [tminfinalclprob,pdifference,tfinalsearchstart,tfinalsearchend] = subCalcClassProbBasedSwitchPoint();
    tminfinaldist = subCalcDistanceBasedSwitchPoint();

    figure(124);
    clf;
    subplot(1,2,1);
    imagesc(self.leakyAvgFrame);
    hold on;
    t1 = tstart-20;
    tt = t1:self.currentFrame;
    plotIdentityIds = assignedIdentityIds; 
    pold =orgpos(:,plotIdentityIds,tt); % old pos
    plot(squeeze(pold(1,:,:))',squeeze(pold(2,:,:))','x-','linewidth',1.5);
    title('before + DAG (square)');
    pdag = dagpos;
    plot(squeeze(pdag(1,:,:))',squeeze(pdag(2,:,:))','s:','linewidth',1);
    

    subplot(1,2,2);
    imagesc(self.leakyAvgFrame);
    hold on;
    pnew = self.pos(:,plotIdentityIds,tt); % new pos
    plot(squeeze(pnew(1,:,:))',squeeze(pnew(2,:,:))','x-') ;
    title('after');

    figure(125);
    clf;
    xlim(tt([1,end])-t1+1)
    subplot(2,1,1);
    plot(pdifference,'-x')
    hold on;
    plot(tminfinal([1,1]) - t1 + 1,ylim(),'-k');
    plot(tminfinalclprob([1,1]) - t1 + 1,ylim(),'--r');
    plot(tminfinaldist([1,1]) - t1 + 1,ylim(),'--m');

    plot([tfinalsearchstart,tfinalsearchend]-t1 + 1,[0,0],'r-o','linewidth',2);
    if crossingflag
      plot([min(tfirst),max(tlast)]-t1 + 1,[0,0],'m-x','linewidth',2);        
    end

    subsubplot(2,1,2,3,1,1);
    plot(tt-t1+1,squeeze(pold(1,:,:))','-','linewidth',2)
    title('old');
    ylabel('x');

    subsubplot(2,1,2,3,1,2);
    plot(tt-t1+1,squeeze(pnew(1,:,:))','-','linewidth',2)
    title('new');
    ylabel('x');
    xl = xlim;

    subsubplot(2,1,2,3,1,3);
    plot((self.currentFrame-size(pdag,3)+1:self.currentFrame)-t1+1,squeeze(pdag(1,:,:))','-','linewidth',2);
    title('DAG')
    ylabel('x');
    xlabel('rel. Frame');
    xlim(xl);


    if crossingflag
      % first frames of crossings are always switch points
      tswitch =  unique(tfirst);
      % as well as the last (which should be the same for all because
      % crossings get extended for all participating tracks);
      tswitch =  [tswitch,unique(tlast)];
    else
      tswitch = [];
    end

    % add switch points based on distance
    oldpos = orgpos;
    localpos = permute(oldpos(:,oldIdentityIds,tstart:tend),[3,2,1]);

    xd = bsxfun(@minus,localpos(1:end-1,:,1),permute(localpos(2:end,:,1),[1,3,2]));
    yd = bsxfun(@minus,localpos(1:end-1,:,2),permute(localpos(2:end,:,2),[1,3,2]));
    d = sqrt(xd.^2 + yd.^2);
    [~,crosses] = min(d,[],3);
    crossesmsk = crosses - ones(size(crosses,1),1)*(1:size(crosses,2));
    crosses(~crossesmsk) = 0;

    tswitch = [tswitch, find(any(crosses,2))'+tstart-1];

    tswitch = unique(tswitch);

    ttt = tswitch - t1 + 1;

    subplot(2,1,1);
    yl = ylim;
    hold on;
    plot([ttt(:),ttt(:)]',yl'*ones(1,length(ttt)),':b');
    

  end


  function [dagpos,dagf2t,flag] = subGetDagBasedPos()
  % get the backtrace from  tstart (if possible). 
  
    flag = 0;
    stepsback = self.currentFrame-tstart+1;
    [dagpos,dagf2t] = self.daGraph.backtrace(assignedIdentityIds,trackIndices,[],stepsback);
    dagpos = permute(dagpos,[1,3,2]);

    % object index does not correspond to the objindex in dag if
    % deletions happen. Need to save trackids in dag to handle this
    assert(~self.opts.tracks.withTrackDeletion)
    
    % CAUTION: might involve other identity, so that f2t has not unique
    % tracks for aeach time frame anymore. 
    u = unique([unique(dagf2t(:))',assignedIdentityIds]);
    Nindiv = accumarray(dagf2t(:),1,[self.nindiv,1]);
    if length(u)>length(assignedIdentityIds)
      extra = setdiff(u,assignedIdentityIds);
      if sum(Nindiv(extra))/sum(Nindiv)/length(assignedIdentityIds) > 0.1
        if self.verbosity>2
          xy.helper.verbose('CAUTION: More tracks involved in crossing !')
        end
        
        flag=0; % ignore
      end
    end
    
    if min(Nindiv(assignedIdentityIds))<stepsback/2
      if self.verbosity>2
        xy.helper.verbose('CAUTION: Tracks seem to have heavy overlap in DAG !');
      end
      flag = 2;
    end
    if ~all(self.identityId2TrackId(tcurrent-stepsback+1,assignedIdentityIds)==dagf2t(1,:))
      if self.verbosity>2
        xy.helper.verbose('CAUTION: DAG initial ordering is different!');
      end
      flag = 3;
    end

  end
  
  
  
  function [tminOut,pdiff,tsearchstart,tsearchend] = subCalcClassProbBasedSwitchPoint()
  % get switch time time from classification accuracy 
    
    tt = tstart:tcurrent;
    t1 = tstart;

    pdiff = [];
    tspan = tcurrent-tstart + 1;

    for ii = 1:length(idxTrackIndices)
      trIdx = idxTrackIndices(ii);
      if tspan>self.tracks(trIdx).classProbHistory.nHistory
        tminOut = tcurrent; % cannot compute
        tsearchstart = t1; % dummy
        tsearchend = tcurrent;
        return;
      end
      [p,w] = self.tracks(trIdx).classProbHistory.getData(tspan);
      if isempty(pdiff)
        pdiff = zeros(length(w),length(idxTrackIndices));
      end
      pdiff(:,ii) = p(:,idxOldIdentityIds(ii)) - p(:,idxAssignedIdentityIds(ii));
      %rt =self.tracks(trIdx).classProbHistory.reasonableThres;
      pdiff(w==0,ii) = NaN;
    end


    if crossingflag
      % currently crossing
      % if currently crossing restrict on the region between
      % the crossings + nearby lost points
      
      tfirstlocal = min(tfirst);
      tlastlocal = max(tlast);
      srel = ones(1,3);
      crossingMsk = (tt>=(tfirstlocal) & tt<=(tlastlocal))';
      nanmsk = any(isnan(pdiff),2);
      crossingMsk = crossingMsk | imerode(imdilate(nanmsk,srel),srel);

      tlast1 = tlastlocal-t1+1;
      tfirst1 = tfirstlocal -t1 + 1;

      % find continuous region around the last crossings
      findidx = find(~crossingMsk(tlast1+1:end),1,'first');
      if ~isempty(findidx)
        tsearchend1 = findidx+ tlast1-1;
      else
        tsearchend1 = tcurrent - t1 + 1;
      end
      
      findidx = find(~crossingMsk(tfirst1-1:-1:1),1,'first');
      if ~isempty(findidx)
        tsearchstart1 = tfirst1 - findidx + 1;
      else
        tsearchstart1 = tstart - t1 + 1;
      end
      
    else
      tsearchstart1 = 1;
      tsearchend1 = tcurrent - t1 + 1;
    end
    tsearchstart = tsearchstart1 + t1 -1;
    tsearchend = tsearchend1 + t1 -1;
    
    % determine where last time above zero within region
    probmsk = nanmean(pdiff,2)>0 | isnan(nanmean(pdiff,2)); % or all nan
    srel = [1,1,1]; % delete single frames above zero;
    probmsk = imdilate(imerode(probmsk,srel),srel);
    
    lastidx = find(probmsk(tsearchstart1:tsearchend1),1,'last');
    if isempty(lastidx)
      tminOut = tsearchstart1 + t1 - 1;
    else
      tminOut = lastidx + tsearchstart1 + t1 - 1; %
    end
    if self.verbosity>2
      xy.helper.verbose('N[P] = %d',tminOut-tcurrent);
    end
    
  end % nested function
    
    
  function [tminOut,nanmskOut] = subCalcLostTrackBasedSwitchPoint()
    % candidate switch points where tracks were lost

    localpos = permute(self.pos(:,idxOldIdentityIds,tstart:tend),[3,2,1]);
    localpos(isnan(localpos)) = 0;
    sz = size(localpos);
    sz(1) = 1;
    localpos1 = cat(1,zeros(sz),localpos,-ones(sz));
    msk1 = all(~diff(localpos1,[],1),3);
    msk = diff(msk1==1);
    tcandidates = find(any(msk,2))' + 1 + tstart -1;
    
    tc1 = min(tcandidates-tstart+1,size(localpos,1));
    [~,change] = ismember(idxOldIdentityIds,idxAssignedIdentityIds);
    d = localpos(tc1,:,:) - localpos(tc1,change,:);
    [m, tidx1] = nanmin(sum(sqrt(d(:,:,1).^2 + d(:,:,2).^2),2));
    nanmskOut = msk | msk1(1:end-1,:);

    if ~isnan(m)
      tswitch1 = tc1(tidx1);
      %localpos(cat(3,msk,msk)~=0 | cat(3,msk1(1:end-1,:),msk1(1:end-1,:) )) = NaN; 
      %finalpos = cat(1,localpos(1:tswitch31-1,:,:),localpos(tswitch31:end,change,:));
      tminOut = tswitch1 +tstart -1;
      if self.verbosity>2
        xy.helper.verbose('N[L] = %d',tminOut-tcurrent);
      end
      
    else 
      tminOut = [];
    end

  end
    
  
  function tminOut = subCalcDistanceBasedSwitchPoint()
  % to be called from switchIdentity: calculated the end point based on the distance of
  % the tracks. It is assumed that the given indices are only ONE permutation


  % build Gaussian kernel
  %sigma = self.clpMovAvgTau; % same for dist/clprob based
  %dt = 1; % in frames
  %l = 0:dt:sigma*8;
  %l = [l(end:-1:2) l];
  %G = (exp(-(l.^2/2/sigma^2))/sqrt(2*pi)/sigma*dt)';
  %g2 = ceil(length(G)/2);
    

  % check min distance within crossing
    tt = tstart:tend;
    localpos = permute(self.pos(:,idxOldIdentityIds,tt),[3,2,1]);
    %some positions might be nan. 
    localpos(cat(3,delmsk,delmsk)) = NaN;

    for i_inter = 1:size(localpos,2)
      inter_msk = isnan(localpos(:,i_inter,1));
      if sum(~inter_msk)>2
        if any(inter_msk) 
          localpos(inter_msk,i_inter,:) = ...
              interp1(tt(~inter_msk),localpos(~inter_msk,i_inter,:),tt(inter_msk),'linear','extrap');
        end     
      else
        % all nan...!?!
        if all(inter_msk)
          localpos(inter_msk,i_inter,:) = 0;
          if length(idxOldIdentityIds)==2
            % cannot be distance based for two identity. just take start
            tminOut =  tstart;
            return;
          end
        else
          % just one non-nan
          for i_xy = 1:2
            localpos(:,i_inter,i_xy) = nanmean(localpos(:,i_inter,i_xy));
          end
        end
      end
    end
    
    [~,change] = ismember(idxOldIdentityIds,idxAssignedIdentityIds);
    dist = sum((localpos - localpos(:,change,:)).^2,3);
    %dist = squeeze(sum((self.pos(:,idxOldIdentityIds,tstart:tend) - self.pos(:,idxAssignedIdentityIds,tstart:tend)).^2,1))';
    %dist = [repmat(dist(1,:),g2,1);dist;repmat(dist(end,:),g2,1)];
    %cdist = conv2(dist,G);
    %cdist = cdist(2*g2:end-2*g2+1,:);
    mdist = mean(dist,2); %just one component and symmetric distance. (DOES NOT
                          %ACCOUNT FOR MULTIPLE CROSSINGS)
    [~,tminsub] = min(mdist);
    tminOut = tstart + tminsub-1;
    
    if self.verbosity>2
      xy.helper.verbose('N[D] = %d',tminOut-tcurrent);
    end
    
  end % nested function
    
end
  