function generateResults(self)
% generates the correct  "res" structure 
  
  if isempty(self.savedTracks.id) 
    fish.helper.verbose('WARNING: cannot generate results!')
    fish.helper.verbose('WARNING: not all fish detected. Maybe adjust "nfish" setting.');
    return;
  end

  if ~isempty(self.savedTracksFull) 
    warning(['Will NOT process the full track structure. Switching of tracks are thus not ' ...
             'considered. The full tracks results can be accessed in the field ' ...
             '"savedTracksFull"'])
  end

  
  
  nFrames = size(self.savedTracks.id,3)/self.nfish;
  if self.currentFrame~=nFrames
    fish.helper.verbose(['WARNING: %d frames got lost (premature abort while ' ...
             'tracking?)'],self.currentFrame-nFrames);
  end
  self.currentFrame = nFrames;
  
  if self.opts.classifier.onlyDAGMethod
    subGetPosFromDag([],1); % overwrites pos
  end

  fishId2TrackId = self.fishId2TrackId(1:nFrames,:)';
  self.res.swb = subGenerateTracks(fishId2TrackId);
  self.res.swb.pos = permute(self.pos(:,:,1:nFrames),[3,1,2]);      
  self.res.t = self.res.tracks.t(:,1);
  self.res.tabs = self.tabs(1:nFrames,:);

  if isfield(self.savedTracks,'stmInfo')
    sz = size(self.savedTracks.stmInfo);
    d = length(sz);
    self.res.stmInfo = permute(reshape(self.savedTracks.stmInfo,[sz(1:end-1),self.nfish,nFrames]),[d+1,d,2:d-1,1]);
    
    if self.stimulusPresenter.usePredFishId && isfield(self.savedTracks,'predFishId')
      self.res.stmFishId = reshape(self.savedTracks.predFishId, [self.nfish,nFrames]);
    elseif ~self.stimulusPresenter.usePredFishId && isfield(self.savedTracks,'fishId')
      self.res.stmFishId = reshape(self.savedTracks.fishId,[self.nfish,nFrames]);
    end

  end


  if ~self.opts.classifier.onlyDAGMethod
    % also generate dag
    clear fishId2TrackId
    [pos,dagf2t] = subGetPosFromDag();

    % check for big overlaps and correct them with switch-based
    dagf2t = subCorrectDagOverlaps(dagf2t);

    %need to do again
    self.res.dag = subGenerateTracks(dagf2t'); 
    pos = self.res.dag.tracks.centroid;
    self.res.dag.pos = permute(pos,[1,3,2]);      
    
  end
  
  
  function df2t = subCorrectDagOverlaps(df2t)
  
    
    % ASSUMES IDX MAT IS SAME AS ID MSK (no track deletion)

    minOverlap = 21;
    eqmsk = bsxfun(@eq,df2t,permute(df2t,[1,3,2]));

    n = size(eqmsk,1);
    se = ones(minOverlap,1);
    tmp = imdilate(imerode(eqmsk(:,:),se),se);
    tmp = imerode(imdilate(eqmsk(:,:),se),se);
    eqmsk = reshape(tmp,[],self.nfish,self.nfish);

    idx = find(tril(ones(self.nfish),-1));
    ieqmsk = eqmsk(:,idx);
    noverlap = sum(ieqmsk,2);    
    doverlap = sum(eqmsk,3)>2; % double overlaps    
    

    % get the positions that were lost in DAG
    lostmsk = bsxfun(@eq,1:self.nfish,permute(df2t,[1,3,2]));
    tmp = imdilate(imerode(lostmsk(:,:),se),se);
    tmp = imerode(imdilate(tmp,se),se);
    lostmsk = reshape(tmp,[],self.nfish,self.nfish);
    lostmsk = all(~lostmsk,3);

    
    % note that their might be overlaps in DAG too because we used
    % the DAG switching methods... however, they should not be too
    % long (except for fishupdate type of corrections...)
    

    % ignore double and strange losses (more than 1) for now. 
    onelost = sum(lostmsk,2)==1;
    [~,idxlost] = max(lostmsk,[],2);
    %idxlost(~onelost) = self.nfish+1; % no need. only between start:stop
    
    msk = bsxfun(@and,ieqmsk,onelost & all(~doverlap,2));

    % get the cl prob of the lost id
    indl = (1:n)' + (idxlost-1)*n;
    idl = idxlost; % assume that track ID is same as IDX (correct with
                   % no deletion)    
    
    % get the cl prob of the overlapping id
    [~,idxeq] = max(ieqmsk,[],2);
    subeq = fish.helper.i2s([self.nfish,self.nfish],idxeq);
    indo = (1:n)' + (subeq(:,1)-1)*n;
    ido  = df2t(indo);
    
    
    % get prob in the original TrackID(idx) order 
    f = {'classProb'};
    sz = size(self.savedTracks.(f{1}));
    d = length(sz); % at least 3
    tmp = permute(self.savedTracks.(f{1}),[d,2,1,3:d-1]);
    tmp = reshape(tmp,[self.nfish,n,sz(2),sz(1),sz(3:d-1)]);
    classProb = reshape(permute(tmp,[2,1,3:d+1]),[],self.nfish);

    cll = classProb(indl,:);
    %idl = res.swb.tracks.id(indl);

    clo = classProb(indo,:);    
    %ido = res.dag.tracks.id(indo);
    
    for i = 1:length(idx)
      mski = msk(:,i);
      [fid1,fid2] = ind2sub(self.nfish([1,1]),idx(i));
      
      % find onsets and offsets
      d = diff([0;mski;0]);
      stop = find(d==-1)-1;
      start = find(d==1);

      % get prob data
      accmsk = zeros(n,1);
      accmsk(start) = 1;
      accmsk = cumsum(accmsk);
      accmsk(~mski) = length(start)+1;
      
      mclo1 = accumarray(accmsk,clo(:,fid1),[length(start)+1,1],@mean);
      mclo2 = accumarray(accmsk,clo(:,fid2),[length(start)+1,1],@mean);
      mcll1 = accumarray(accmsk,cll(:,fid1),[length(start)+1,1],@mean);
      mcll2 = accumarray(accmsk,cll(:,fid2),[length(start)+1,1],@mean);

      s12 = mclo1 + mcll2;
      s21 = mclo2 + mcll1;
      order12 = s12(1:end-1) >= s21(1:end-1);
      order21 = ~order12;
      %diffprob = abs(s12-s21);
      
      % re-order the results. Just redefine the 
      msk12 = zeros(n+1,1);
      msk12(start(order12)) = 1;
      msk12(stop(order12)+1) = -1;
      msk12 = cumsum(msk12);
      idx12 = find(msk12(1:end-1));

      msk21 = zeros(n+1,1);
      msk21(start(order21)) = 1;
      msk21(stop(order21)+1) = -1;
      msk21 = cumsum(msk21);
      idx21 = find(msk21(1:end-1));

      
      df2t(idx12,fid1) = ido(idx12);
      df2t(idx12,fid2) = idl(idx12);

      df2t(idx21,fid1) = idl(idx21);
      df2t(idx21,fid2) = ido(idx21);
    end
    
  
  
  end
  
  function [res] = subGenerateTracks(f2t);
    res = [];
    trackIdMat = reshape(self.savedTracks.id,self.nfish,nFrames);
    for j = 1:self.nfish
      u = unique(trackIdMat(j,:));
      n = find(~isnan(u));
      L(j) = length(n);
      U(1,j) = u(n(1));
    end

    if all(L==1)
      % short-cut to avoid the loop
      [~,Loc] = ismember(f2t,U);
    else
      % HOW CAN THIS DONE A BIT MORE EFFICIENTLY?
      Loc = zeros(size(f2t));
      for i = 1:nFrames
        [~,Loc(:,i)] = ismember(f2t(:,i),trackIdMat(:,i));
      end
    end

    % fill in the gaps
    order = (1:self.nfish)';
    msk = ~Loc;
    idx = find(any(msk,1));
    for i = 1:length(idx)
      loc = Loc(:,idx(i));
      rest = setdiff(order,loc(~~loc));
      Loc(~loc,idx(i)) = rest;
    end

    fridx = ones(1,self.nfish)' * (1:nFrames);
    idx = fish.helper.s2i(size(Loc),[Loc(:),fridx(:)]);
    for f = fieldnames(self.savedTracks)'
      if isempty(self.savedTracks.(f{1}))
        continue;
      end
      if strcmp(f{1},'stmInfo')
        continue;
      end
      sz = size(self.savedTracks.(f{1}));
      d = length(sz); % at least 3
      trackdat = permute(self.savedTracks.(f{1}),[d,2,1,3:d-1]);
      fishdat = reshape(trackdat(idx,:),[self.nfish,nFrames,sz(2),sz(1),sz(3:d-1)]);
      res.tracks.(f{1}) = permute(fishdat,[2,1,3:d+1]);
    end

    fishid =  (1:self.nfish)' * ones(1,nFrames);
    fishid(msk) = NaN;
    res.tracks.fishId = fishid';
    
  end
  


  function [postrace,trackIdxMat] = subGetPosFromDag(assignedFishId,force)
    
    if nargin>0 && ~isempty(assignedFishId)
      predFishIds = assignedFishId;
    else
      predFishIds = [self.tracks.predFishId]; % use predFish. Can do Better ?!?
    end
    if nargin<2
      force = 0;
    end
    
    %self.daGraph.checkOverlap([],1);
    
    %trackIdx  and trackids SHOULD be the same! (if with no deletion)
    % at laest assert for last tracks (if 1:nfish, all previous should be too)
    assert(all([self.tracks.id] == 1:self.nfish));
    
    % backtrace. 
    [postrace,trackIdxMat] = self.daGraph.backtrace(1:self.nfish,predFishIds);
    postrace = permute(postrace,[1,3,2]);
    
    mt = size(postrace,3);
    t = self.currentFrame-mt+1:self.currentFrame;
    
    if force
      self.fishId2TrackId(t,:) = trackIdxMat;
      self.pos(:,:,t) = postrace;
    else
      pos = self.pos(:,:,1:self.currentFrame);
      pos(:,:,t) = postrace;
      postrace = pos;
      
      f2t = self.fishId2TrackId(1:self.currentFrame,:);
      f2t(t,:) = trackIdxMat;
      trackIdxMat = f2t;
    end
    
    
  end
  

end
  
  

