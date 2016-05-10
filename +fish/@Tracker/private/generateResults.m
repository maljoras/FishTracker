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

  nFrames = self.currentFrame;
  

  if self.opts.classifier.onlyDAGMethod
    self.getPosFromDag([],1); % overwrites pos
  end

  fishId2TrackId = self.fishId2TrackId(1:nFrames,:)';
  self.res = subGenerateTracks(fishId2TrackId);
  self.res.pos = permute(self.pos(:,:,1:nFrames),[3,1,2]);      

  
  if ~self.opts.classifier.onlyDAGMethod
    % also generate dag
    clear fishId2TrackId
    [pos,fishId2TrackId] = self.getPosFromDag();
    self.res.dag = subGenerateTracks(fishId2TrackId');
    
    self.res.dag.pos = permute(pos(:,:,1:nFrames),[3,1,2]);      

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
  
end
  
  




