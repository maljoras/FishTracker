function generateResults(self)
% generates the correct  "res" structure 
  
  if isempty(self.savedTracks.id) 
    xy.helper.verbose('WARNING: cannot generate results!')
    xy.helper.verbose('WARNING: not all identity detected. Maybe adjust "nbody" setting.');
    return;
  end

  if ~isempty(self.savedTracksFull) 
    warning(['Will NOT process the full track structure. Switching of tracks are thus not ' ...
             'considered. The full tracks results can be accessed in the field ' ...
             '"savedTracksFull"'])
  end

  
  for f = fieldnames(self.savedTracks)'
    if iscell(self.savedTracks.(f{1}))
      d =  length(size(self.savedTracks.(f{1}){1}));% at least 3
      self.savedTracks.(f{1}) = cat(d,self.savedTracks.(f{1}){:});
    end
  end
  
  nFrames = size(self.savedTracks.id,3)/self.nbody;
  if self.currentFrame~=nFrames
    xy.helper.verbose(['WARNING: %d frames got lost (premature abort while ' ...
             'tracking?)'],self.currentFrame-nFrames);
  end
  
  self.currentFrame = nFrames;
  self.res = [];
  
  identityId2TrackId = self.identityId2TrackId(1:nFrames,:)';
  self.res.swb = subGenerateTracks(identityId2TrackId);
  self.res.swb.pos = subGeneratePos(self.res.swb);
  tabs = self.tabs(1:nFrames,:);


  % also generate dag
  clear identityId2TrackId
  [pos,dagf2t] = subGetPosFromDag();
  
  % check for big overlaps and correct them with switch-based
  dagf2t = subCorrectDagOverlaps(dagf2t);

  %need to do again
  self.res.dag = subGenerateTracks(dagf2t'); 
  self.res.dag.pos = subGeneratePos(self.res.dag);


  % correct the time for unique distance dt
  dt = 1/self.videoHandler.frameRate;
  tidx = round((tabs-tabs(1))/dt)+1;
  frames = (1:nFrames)';
  t = (0:tidx(end)-1)'*dt + tabs(1);
  for f = {'swb','dag'}
    for f2 = fieldnames(self.res.(f{1}).tracks)'
      field = self.res.(f{1}).tracks.(f2{1});
      sz = size(field);
      sz(1) = length(t);
      tmp = nan(sz);
      assert(length(sz)<7)
      tmp(tidx,:,:,:,:,:) = field; 
      self.res.(f{1}).tracks.(f2{1}) = tmp;
    end

    self.res.(f{1}).t = t;
    self.res.(f{1}).tabs = nan(size(t));
    self.res.(f{1}).tabs(tidx) = tabs;

    self.res.(f{1}).iframe = nan(size(t));
    self.res.(f{1}).iframe(tidx) = frames;
    
    pos = self.res.(f{1}).pos;
    tmp = nan(size(pos));
    tmp(tidx,:,:) = pos;
    self.res.(f{1}).pos = tmp;
  end
  
  
  
  function pos = subGeneratePos(res)
  % gets a new pos from the re-ordered tracks. Do not use the original
  % pos (which might contain Kalman predictions) but centerLine if available
  
    if isfield(res.tracks,'centerLine')
      cl = mean(res.tracks.centerLine,4);
      ce = res.tracks.centroid;
      idx = find(isnan(cl));
      cl(idx) = ce(idx);
      pos = permute(cl,[1,3,2]);
    else
      ce = res.tracks.centroid;
      pos = permute(ce,[1,3,2]);
    end
  
    % delete beyond border pixels
    posx = squeeze(pos(:,1,:));
    posy = squeeze(pos(:,2,:));

    sz = self.videoHandler.frameSize;
    posx(posx>sz(2) | posx<1) = NaN;
    posy(posy>sz(1) | posy<1) = NaN;

    pos(:,1,:) = posx;
    pos(:,2,:) = posy;
  
  end
    
    
  
  function df2t = subCorrectDagOverlaps(df2t)
  
    % ASSUMES IDX MAT IS SAME AS ID MSK (no track deletion)

    MINOVERLAP = ceil(self.videoHandler.frameRate/2); 
    PROBTHRES = 0;%self.maxClassificationProb*self.opts.tracks.probThresForIdentity;
    
    eqmsk = bsxfun(@eq,df2t,permute(df2t,[1,3,2]));

    n = size(eqmsk,1);
    se = ones(MINOVERLAP,1);
    tmp = imdilate(imerode(eqmsk(:,:),se),se);
    %tmp = imerode(imdilate(eqmsk(:,:),se),se);
    eqmsk = reshape(tmp,[],self.nbody,self.nbody);

    idx = find(tril(ones(self.nbody),-1));
    ieqmsk = eqmsk(:,idx);
    noverlap = sum(ieqmsk,2);    
    doverlap = sum(eqmsk,3)>2; % double overlaps    
    

    % get the positions that were lost in DAG
    lostmsk = bsxfun(@eq,1:self.nbody,permute(df2t,[1,3,2]));
    tmp = imdilate(imerode(lostmsk(:,:),se),se);
    %tmp = imerode(imdilate(tmp,se),se);
    lostmsk = reshape(tmp,[],self.nbody,self.nbody);
    lostmsk = all(~lostmsk,3);

    
    % note that their might be overlaps in DAG too because we used
    % the DAG switching methods... however, they should not be too
    % long (except for bodyupdate type of corrections...)
    

    % ignore double and strange losses (more than 1) for now. 
    onelost = sum(lostmsk,2)==1;
    [~,idxlost] = max(lostmsk,[],2);
    %idxlost(~onelost) = self.nbody+1; % no need. only between start:stop
    
    msk = bsxfun(@and,ieqmsk,onelost & all(~doverlap,2));

    % get the cl prob of the lost id
    indl = (1:n)' + (idxlost-1)*n;
    idl = idxlost; % assume that track ID is same as IDX (correct with
                   % no deletion)    
    
    % get the cl prob of the overlapping id
    [~,idxeq] = max(ieqmsk,[],2); % this is now bodyID
    subeq = xy.helper.i2s([self.nbody,self.nbody],idxeq);
    findo = (1:n)' + (subeq(:,1)-1)*n;
    ido  = df2t(findo); % also idx in cl
    indo = (1:n)' + (ido(:,1)-1)*n;
    indo(isnan(indo)) = 1; % nan's will not enter anyway

    % get prob in the original TrackID(idx) order 
    f = {'classProb'};
    sz = size(self.savedTracks.(f{1}));
    d = length(sz); % at least 3
    tmp = permute(self.savedTracks.(f{1}),[d,2,1,3:d-1]);
    tmp = reshape(tmp,[self.nbody,n,sz(2),sz(1),sz(3:d-1)]);
    classProb = reshape(permute(tmp,[2,1,3:d+1]),[],self.nbody);

    cll = classProb(indl,:);
    %idl = res.swb.tracks.id(indl);

    clo = classProb(indo,:);    
    %ido = res.dag.tracks.id(indo);
    
    for i = 1:length(idx)
      mski = msk(:,i);
      [fid1,fid2] = ind2sub(self.nbody([1,1]),idx(i));
      
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
      diffprob = abs(s12-s21);
      probmsk = diffprob<PROBTHRES;
      probmsk = probmsk | max(mcll1,mcll2)<self.maxClassificationProb*self.opts.tracks.probThresForIdentity;
      order12(probmsk) = false;
      order21(probmsk) = false;

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
    trackIdMat = reshape(self.savedTracks.id,self.nbody,nFrames);
    for j = 1:self.nbody
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
    order = (1:self.nbody)';
    msk = ~Loc;
    idx = find(any(msk,1));
    for i = 1:length(idx)
      loc = Loc(:,idx(i));
      rest = setdiff(order,loc(~~loc));
      Loc(~loc,idx(i)) = rest;
    end

    fridx = ones(1,self.nbody)' * (1:nFrames);
    idx = xy.helper.s2i(size(Loc),[Loc(:),fridx(:)]);
    for f = fieldnames(self.savedTracks)'
      if isempty(self.savedTracks.(f{1}))
        continue;
      end
% $$$       if strcmp(f{1},'stmInfo')
% $$$         continue;
% $$$       end
      sz = size(self.savedTracks.(f{1}));
      d = length(sz); % at least 3
      trackdat = permute(self.savedTracks.(f{1}),[d,2,1,3:d-1]);
      bodydat = reshape(trackdat(idx,:),[self.nbody,nFrames,sz(2),sz(1),sz(3:d-1)]);
      res.tracks.(f{1}) = permute(bodydat,[2,1,3:d+1]);
    end

    bodyid =  (1:self.nbody)' * ones(1,nFrames);
    bodyid(msk) = NaN;
    res.tracks.identityId = bodyid';
    
  end
  


  function [postrace,trackIdxMat] = subGetPosFromDag(assignedIdentityId,force)
    
    if nargin>0 && ~isempty(assignedIdentityId)
      predIdentityIds = assignedIdentityId;
    else
      predIdentityIds = [self.tracks.predIdentityId]; % use predBody. Can do Better ?!?
    end
    if nargin<2
      force = 0;
    end
    
    %self.daGraph.checkOverlap([],1);
    
    %trackIdx  and trackids SHOULD be the same! (if with no deletion)
    % at laest assert for last tracks (if 1:nbody, all previous should be too)
    assert(all([self.tracks.id] == 1:self.nbody));
    
    % backtrace. 
    [postrace,trackIdxMat] = self.daGraph.backtrace(1:self.nbody,predIdentityIds);
    postrace = permute(postrace,[1,3,2]);
    
    mt = size(postrace,3);
    t = self.currentFrame-mt+1:self.currentFrame;
    
    if force
      self.identityId2TrackId(t,:) = trackIdxMat;
      self.pos(:,:,t) = postrace;
    else
      pos = self.pos(:,:,1:self.currentFrame);
      pos(:,:,t) = postrace;
      postrace = pos;
      
      f2t = self.identityId2TrackId(1:self.currentFrame,:);
      f2t(t,:) = trackIdxMat;
      trackIdxMat = f2t;
    end
    
    
  end
  

end
  
  

