function switchFish(self,trackIndices,assignedFishIds,crossingflag)
% changes the FishIds of the tracks. expects a premutations. updates all relevant fishId
% dependet varibeles. Resets the data to avoid the learning of wrong features. Also
% resets the uniqueFishFrame counter
%
% if the crossing flag is set, the switching point is calculated based on the
% distance within the last crossing. If not set, the switchingpoint is calculated
% based on the classProbHistory. 


  assert(length(trackIndices)==length(assignedFishIds));
  oldFishIds = [self.tracks(trackIndices).fishId];
  assert(all(sort(oldFishIds)==sort(assignedFishIds)));

  if length(trackIndices)==1 || all(oldFishIds==assignedFishIds)
    return; % no need to do anything
  end

  self.uniqueFishFrames = 0;

  %  it should be guranteed that it is a true permutation !! Do it anyway...
  [nc ncidx]= self.connectedComponents(trackIndices,assignedFishIds);


  %% switch the tracks for each connected component
  f2t = self.fishId2TrackId;
  orgpos = self.pos;
  self.resetBatchIdx(trackIndices); % reset for all

  if self.opts.display.switchFish && self.displayif && self.opts.display.tracks;
    self.displayCurrentTracks();
  end

  for i_ncidx = 1:length(ncidx)  % actually do not need to check again.
    idx = ncidx{i_ncidx};

    % define common pars for subfunctions
    idxOldFishIds = [self.tracks(trackIndices(idx)).fishId]; 
    idxTrackIndices = trackIndices(idx);
    idxAssignedFishIds = assignedFishIds(idx);
    tfirst =[self.tracks(idxTrackIndices).firstFrameOfCrossing]; 
    tlast = [self.tracks(idxTrackIndices).lastFrameOfCrossing]; 
    tcurrent = self.currentFrame;
    
    if crossingflag
      tstart = max(min(tfirst) - self.opts.classifier.nFramesAfterCrossing,1);
      tend = min(self.currentFrame,max(tlast) + self.opts.classifier.nFramesAfterCrossing);
    else
      tstart = max(max(tlast) - 2*self.opts.classifier.nFramesAfterCrossing,1) ; % put about last crossing in...
      tend = tcurrent;
    end
    
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
      % fishClassUpdate
      tminprob = subCalcClassProbBasedSwitchPoint();

      if tcurrent-tminprob> 2*self.opts.classifier.nFramesAfterCrossing
        tminfinal = tminprob;
      else
        if length(idx)>2
          tminfinal = subCalcDistanceBasedSwitchPoint();
        else
          tminfinal = tminfinallost;
        end
      end
    end
    
    
    if crossingflag 
      % delete some of the critical points
      tcross = tstart:tend;
      p = orgpos(:,oldFishIds(idx),tcross);
      p(:,delmsk') = NaN;
      orgpos(:,oldFishIds(idx),tcross) = p;
    end
    

    % get the backtrace up to tminfinal (if possible). Same dims as pos
    % order of oldfishID for plotting
    rtrace = permute(self.daGraph.backtracePositions(assignedFishIds(idx),orgpos(:,oldFishIds(idx),self.currentFrame,self.currentFrame-tfirst+1)),[1,3,2]);
    
    
    %% swap
    t = tminfinal:self.currentFrame;        
    for j = idx(:)'
      self.fishId2TrackId(t,assignedFishIds(j)) = f2t(t,oldFishIds(j));
      self.pos(:,assignedFishIds(j),t) = orgpos(:,oldFishIds(j),t);
      self.tracks(trackIndices(j)).fishId = assignedFishIds(j);
      self.tracks(trackIndices(j)).switchedFrames = tminfinal-self.currentFrame;
    end
    
    %if length(idx)>2
    %warning(['Higher order permutation! Tracks might got mixed up during the time ' ...
    %         'of the crossing.'])
    %tcross = tstart:tend;
    %self.pos(:,assignedFishIds(idx),tcross) = NaN; % too conservative
    %end
    
    if self.opts.display.switchFish && self.displayif
      subPlotCrossings();
    end
    
  end % for ncidx


  function subPlotCrossings();

    [tminfinalclprob,pdifference,tfinalsearchstart,tfinalsearchend] = subCalcClassProbBasedSwitchPoint();
    tminfinaldist = subCalcDistanceBasedSwitchPoint();

    figure(124);
    clf;
    subplot(1,2,1);
    imagesc(self.leakyAvgFrame);
    hold on;
    t1 = tstart-20;
    tt = t1:self.currentFrame;
    plotFishIds = assignedFishIds(idx); 
    pold =orgpos(:,plotFishIds,tt); % old pos
    plot(squeeze(pold(1,:,:))',squeeze(pold(2,:,:))','x-','linewidth',1.5);
    title('before + DAG (square)');
    pdag = rtrace;
    plot(squeeze(pdag(1,:,:))',squeeze(pdag(2,:,:))','s:','linewidth',1);
    

    subplot(1,2,2);
    imagesc(self.leakyAvgFrame);
    hold on;
    pnew = self.pos(:,plotFishIds,tt); % new pos
    plot(squeeze(pnew(1,:,:))',squeeze(pnew(2,:,:))','x-') ;
    title('after');

    figure(125);
    clf;
    xlim([tt([1,end])]-t1+1)
    subplot(2,1,1);
    plot(pdifference,'-x')
    hold on;
    plot([tminfinal([1,1])] - t1 + 1,ylim(),'-k');
    plot([tminfinalclprob([1,1])] - t1 + 1,ylim(),'--r');
    plot([tminfinaldist([1,1])] - t1 + 1,ylim(),'--m');

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



    % TODO: 
    % - determine multiple switching points within a crossing based on
    % distance permutation points.
    % -  problem is: some of the not changed tracks might have multiple
    % switching !!!! how to account for this ?
    % - search for multiple switchpoints within crossing


    %function tswitch = subGetSwitchPoints()
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
    localpos = permute(oldpos(:,idxOldFishIds,tstart:tend),[3,2,1]);

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
      pdiff(:,ii) = p(:,idxOldFishIds(ii)) - p(:,idxAssignedFishIds(ii));
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

    fish.helper.verbose('N[P] = %d',tminOut-tcurrent);

  end % nested function
    
    
  function [tminOut,nanmskOut] = subCalcLostTrackBasedSwitchPoint()
    % candidate switch points where tracks were lost

    localpos = permute(self.pos(:,idxOldFishIds,tstart:tend),[3,2,1]);;
    localpos(isnan(localpos)) = 0;
    sz = size(localpos);
    sz(1) = 1;
    localpos1 = cat(1,zeros(sz),localpos,-ones(sz));
    msk1 = all(~diff(localpos1,[],1),3);
    msk = diff(msk1==1);
    tcandidates = find(any(msk,2))' + 1 + tstart -1;
    
    tc1 = min(tcandidates-tstart+1,size(localpos,1));
    [~,change] = ismember(idxOldFishIds,idxAssignedFishIds);
    d = localpos(tc1,:,:) - localpos(tc1,change,:);
    [m, tidx1] = nanmin(sum(sqrt(d(:,:,1).^2 + d(:,:,2).^2),2));
    nanmskOut = msk | msk1(1:end-1,:);

    if ~isnan(m)
      tswitch1 = tc1(tidx1);
      %localpos(cat(3,msk,msk)~=0 | cat(3,msk1(1:end-1,:),msk1(1:end-1,:) )) = NaN; 
      %finalpos = cat(1,localpos(1:tswitch31-1,:,:),localpos(tswitch31:end,change,:));
      tminOut = tswitch1 +tstart -1;
      
      fish.helper.verbose('N[L] = %d',tminOut-tcurrent);
    else 
      tminOut = [];
    end

  end
    
  
  function tminOut = subCalcDistanceBasedSwitchPoint()
  % to be called from switchFish: calculated the end point based on the distance of
  % the tracks. It is assumed that the given indices are only ONE permutation


  % build Gaussian kernel
  %sigma = self.opts.classifier.clpMovAvgTau; % same for dist/clprob based
  %dt = 1; % in frames
  %l = 0:dt:sigma*8;
  %l = [l(end:-1:2) l];
  %G = (exp(-(l.^2/2/sigma^2))/sqrt(2*pi)/sigma*dt)';
  %g2 = ceil(length(G)/2);
    

  % check min distance within crossing
    tt = tstart:tend;
    localpos = permute(self.pos(:,idxOldFishIds,tt),[3,2,1]);
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
          if length(idxOldFishIds)==2
            % cannot be distance based for two fish. just take start
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
    
    [~,change] = ismember(idxOldFishIds,idxAssignedFishIds);
    dist = sum((localpos - localpos(:,change,:)).^2,3);
    %dist = squeeze(sum((self.pos(:,idxOldFishIds,tstart:tend) - self.pos(:,idxAssignedFishIds,tstart:tend)).^2,1))';
    %dist = [repmat(dist(1,:),g2,1);dist;repmat(dist(end,:),g2,1)];
    %cdist = conv2(dist,G);
    %cdist = cdist(2*g2:end-2*g2+1,:);
    mdist = mean(dist,2); %just one component and symmetric distance. (DOES NOT
                          %ACCOUNT FOR MULTIPLE CROSSINGS)
    [~,tminsub] = min(mdist);
    tminOut = tstart + tminsub-1;
    
    
    fish.helper.verbose('N[D] = %d',tminOut-tcurrent);
    
  end % nested function
    
end
  