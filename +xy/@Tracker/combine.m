function combinedObj = combine(self,varargin)
  % Combination of two (or more) tracking objects with overlap

  % make sure that the objects do not refer to the same data (i.e. have
  % different handle references)

  if any(self==cat(1,varargin{:}))
    error('cannot combine handles with identical data references');
  end

  dt = 1/self.videoHandler.frameRate;
  tbinsize = max(self.bodylength/self.maxVelocity,2*dt);

  objs = {self,varargin{:}};

  % check for empty results
  emptymsk = cellfun(@(x)isempty(x.getTrackingResults()),objs);
  objs(emptymsk) = [];

  if isempty(objs)
    error('no results available. First track');
  end


  %first sort according to the starting time
  starttimes = cellfun(@(x)x.timerange(1),objs);
  [~,sidx] = sort(starttimes,'ascend');
  objs = objs(sidx);

  % combine results, that is results are appended. 
  combinedObj = objs{1};
  combinedRes = getTrackingResults(combinedObj);
  for i = 2:length(objs)
    obj = objs{i};
    
    % assert that two video files are the same
    assert(strcmp(combinedObj.videoHandler.videoFile,obj.videoHandler.videoFile));
    assert(combinedObj.nindiv==obj.nindiv);

    res = getTrackingResults(obj);

    if obj.timerange(1)> combinedObj.timerange(2)
      % no overlap
      xy.helper.verbose('WARNING: no overlap. Identity IDs might get mixed up!!');
      keyboard
      % just append
      combinedRes.tracks = cat(1,combinedRes.tracks,res.tracks);
      combinedRes.pos = cat(1,combinedRes.pos,res.pos);

    else
      % some overalp
      toverlap = [obj.timerange(1),combinedObj.timerange(2)];
      
      
      % CAUTION : DISCARD THE TIMES WHERE THE ID classifier is not yet initialized !!!!!!!!!!!!1
      
      track1 = res.tracks(:,1);
      tObj = cat(1,track1.t);
      
      track1 = combinedRes.tracks(:,1);
      tCombined = cat(1,track1.t);
      
      % do not consider positions where assignments costs are high (see getTrackingResults())
      overlapedIdx = find(tCombined(end)>= tObj &  ~any(isnan(res.pos(:,1,:)),3));
      overlapedCombinedIdx = find(tCombined>= tObj(1) & ~any(isnan(combinedRes.pos(:,1,:)),3));

      if isempty(overlapedIdx)|| isempty(overlapedCombinedIdx)
        warning(['no valid indices found for overlap. Identity IDs might get mixed up!!. Take all ' ...
                 'available. ']);

        % simply take all res
        res = obj.res;
        combinedRes = combinedObj.res;
        
        overlapedIdx = find(tCombined(end)>= tObj &  ~any(isnan(res.pos(:,1,:)),3));
        overlapedCombinedIdx = find(tCombined>= tObj(1) & ~any(isnan(combinedRes.pos(:,1,:)),3));
      end
      
      overlappedPos = res.pos(overlapedIdx,:,:);
      overlappedCombinedPos = combinedRes.pos(overlapedCombinedIdx,:,:);
      
      dim = size(overlappedPos,2); % 2-D for now;
      tstart = tCombined(overlapedCombinedIdx(1));
      tidxCombined = floor((tCombined(overlapedCombinedIdx) - tstart)/tbinsize)+1;
      [X,Y] = ndgrid(tidxCombined,1:self.nindiv*dim);
      overlappedCombinedPosInterp = reshape(accumarray([X(:),Y(:)],overlappedCombinedPos(:),[],@mean),[],obj.nindiv);
      
      tidx = min(max(floor((tObj(overlapedIdx) - tstart)/tbinsize)+1,1),tidxCombined(end));
      [X,Y] = ndgrid(tidx,1:self.nindiv*dim);
      overlappedPosInterp = reshape(accumarray([X(:),Y(:)],overlappedPos(:),[tidxCombined(end),self.nindiv*dim],@mean),[],obj.nindiv);

      % now the two position vectors can be compared. 
      cost = pdist2(overlappedCombinedPosInterp',overlappedPosInterp','correlation');

      % use the hungarian matching assignments from the vision toolbox
      assignments = xy.helper.assignDetectionsToTracks(cost,1e4);
      [~,sidx] = sort(assignments(:,1),'ascend');
      assignments = assignments(sidx,:); % make sure that the (combined) tracks are in squential order
                                         % append;
      combinedRes.pos = cat(1,combinedRes.pos,res.pos(overlapedIdx(end)+1:end,:,assignments(:,2)));
      combinedRes.tracks = cat(1,combinedRes.tracks,res.tracks(overlapedIdx(end)+1:end,assignments(:,2)));
    end

    combinedObj.res = combinedRes;
    combinedObj.timerange= [combinedObj.timerange(1),obj.timerange(2)];
    combinedObj.currentTime = self.timerange(1);
    combinedObj.pos = []; % not valid 
    combinedObj.identityId2TrackId = []; % not valid 
  end

end

