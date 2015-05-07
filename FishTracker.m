classdef FishTracker < handle;
  
  
  
  properties 


    
    nfish = 1;
    timerange= [];
    visreader = 0;
    scalecost = [10,10,0,0,0,0];
   
    %options for the parts
    opts = struct([]);
    displayif = 2;
    dt = 1;

    fishlength = [];
    fishwidth = [];
    %
    res = []; % output structure
    pos = [];
  end

  properties (SetAccess = private)

    medianCost = ones(5,1);
    costinfo= {'Location','Pixel overlap','Classifier','MeanIntensity','Orientation'};
    detector = [];
    blobAnalyzer = [];
    fishClassifier = [];
    reader = [];
    videoPlayer = [];
    currentTime = 0;
    tracks = [];
    fishId2TrackId = [];
  end
  
  properties (SetAccess = private, GetAccess=private);
  
    nextId = 1; % ID of the next track

  
    uframe = [];
    frame = [];
    
    segments = [];
    features = [];
    bboxes = [];
    centroids = [];
    idfeatures = [];
    classProb = [];
    crossings = [];

    
    assignments = [];
    unassignedTracks = [];
    unassignedDetections = [];
    cost = [];
    X = [];
    y = [];
    uniqueFishFrames = 0;
    currentFrame = [];
  end
  
  
  methods (Access=private)
    
  

    
    function [reader] = newVideoReader(self,vidname)
      
      if ~isempty(self.timerange)  && self.visreader
        warning('Timerange not supported with visreader');
        self.visreader = 0; % NOT SUPPORTED
      end
      if isempty(self.timerange)
        self.timerange = [0,Inf];
      end
      

      if self.visreader
        reader = vision.VideoFileReader(vidname);
        isDone(reader); % check
        self.visreader = 1;
      else
        reader = VideoReader(vidname);
        self.visreader = 0;
      end
      
      if ~self.visreader
        self.dt = 1/double(reader.FrameRate);
        assert(length(self.timerange)==2 && self.timerange(1)<reader.Duration/self.dt);
        reader.CurrentTime = self.timerange(1);
        self.currentTime = self.timerange(1);

      else
        tmp = info(reader);
        self.dt = 1/double(tmp.VideoFrameRate);
      end

      
    end

    
    function predictor = newPredictor(self,centroid);
      % Create a Kalman filter object.
      predictor = configureKalmanFilter('ConstantVelocity', centroid, [100, 10], [50, 5], 10);
    end
  
    function detector = newForegroundDetector(self,varargin);
      % The foreground detector is used to segment moving objects from
      % the background. It outputs a binary mask, where the pixel value
      % of 1 corresponds to the foreground and the value of 0 corresponds
      % to the background.

      % might become overloaded
      detector = MyForegroundDetector(varargin{:});
    end
    
    
   
    
    function stepVideoPlayer(self,uframe);
      self.videoPlayer.step(uframe);
    end
    
    
    
    function blobAnalyzer=newBlobAnalyzer(self);
      % Connected groups of foreground pixels are likely to correspond to moving
      % objects.  The blob analysis system object is used to find such groups
      % (called 'blobs' or 'connected components'), and compute their
      % characteristics, such as area, centroid, and the bounding box.

      blobAnalyzer = MyBlobAnalysis('fishlength',self.fishlength,'fishwidth',self.fishwidth,...
                                    'maxextent',2*(self.fishlength+self.fishwidth),...
                                    'minextent',0.1*(self.fishlength+self.fishwidth));
        
    end
      
        
    function self=setupSystemObjects(self,vid)
      % Initialize Video I/O
      % Create objects for reading a video from a file, drawing the tracked
      % objects in each frame, and playing the video.

        
      %% Create a video file reader.
      self.reader = self.newVideoReader(vid);
      nfirstframes = 1;
      verbose('Reading first %d,frames of the video to establish background',nfirstframes);
      for i = 1:nfirstframes
        self.readFrame();
        if i==1
          firstframe = zeros([size(self.frame),nfirstframes],'like', self.frame);
        end
        firstframe(:,:,:,i) = self.frame;
      end
      delete(self.reader); % detroy again
      self.reader = self.newVideoReader(vid);

      if self.displayif
        self.videoPlayer = vision.VideoPlayer();
      end
      
      %% get new detector
      args = {};
      for f = fieldnames(self.opts.detector)'
        args{end+1} = f{1};
        args{end+1} = self.opts.detector.(f{1});
      end
      self.detector = newForegroundDetector(self,args{:});
      
      [fishlength,fishwidth] = self.detector.init(firstframe,self.nfish);

      if isempty(self.fishlength) 
        self.fishlength = ceil(1.5*fishlength);
      end
      if isempty(self.fishwidth) 
        self.fishwidth = ceil(2*fishwidth);
      end
      if isempty(self.opts.classifier.maxCrossingDistance)
        self.opts.classifier.maxCrossingDistance = self.fishlength; 
      end
      if isempty(self.opts.classifier.minCrossingDistance)
        self.opts.classifier.minCrossingDistance = self.fishwidth/2; 
      end
      
      %% get new blobAnalyzer
      self.blobAnalyzer = newBlobAnalyzer(self);
      for f = fieldnames(self.opts.blob)'
        if isprop(self.blobAnalyzer,f{1})
          self.blobAnalyzer.(f{1}) = self.opts.blob.(f{1});
        end
      end
    

      %% get new fish classifier 
      self.fishClassifier = newFishClassifier(self);
      % set proporties
      for f = fieldnames(self.opts.classifier)'
        if isprop(self.fishClassifier,f{1})
          self.fishClassifier.(f{1}) = self.opts.classifier.(f{1});
        end
      end
      
      self.fishId2TrackId = nan(250,self.nfish);
      
      
    end

    
    
    function setCurrentTime(self,time)
      if self.visreader
        self.reader.reset();
      else
        self.reader.CurrentTime = time;
      end
      self.currentTime = time;
    end
    
    function readFrame(self);      
 
      if self.visreader
        frame = self.reader.step();
        self.uframe = im2uint8(frame);
      else
        frame = self.reader.readFrame();
        self.uframe = frame;
        frame = im2double(frame);
      end
      self.frame = frame;
      self.currentTime = self.currentTime + self.dt;
    end
    
    function bool = hasFrame(self);      

      if self.currentTime>self.timerange(2) 
        bool = false;
      else
        
        if ~self.visreader
          bool = hasFrame(self.reader);
        else
          bool = ~isDone(self.reader);
        end
      end

    end
    
    
    function  detectObjects(self)
    % calls the detectors
    
      [mask,Iframe] = self.detector.step(self.frame);

      segm = self.blobAnalyzer.step(mask,Iframe,self.frame);
      
      if ~isempty(segm)
        self.centroids = cat(1,segm.Centroid);
        for i = 1:length(segm)
          % better take MSER regions if existent
          if ~isempty(segm(i).MSERregions)
            self.centroids(i,:) = double(segm(i).MSERregionsOffset + segm(i).MSERregions(1).Location);
          end
        end
        
        self.bboxes = int32(cat(1,segm.BoundingBox));
        self.features = [cat(1,segm.MeanIntensity), cat(1,segm.Orientation)];
        self.segments = segm;
        self.idfeatures =  permute(cat(4,segm.fishFeature),[4,1,2,3]);
        self.classProb = self.computeClassProb(self.idfeatures(:,:));
      
      else
        self.centroids = [];
        self.bboxes = [];
        self.segments = [];
        self.idfeatures = [];
      end
      
        

    end
  
    %%% ClASSIFICATION
      
    function cla = newFishClassifier(self,varargin);
      % Create Classifier objects
      featdim = [self.blobAnalyzer.featureheight,self.blobAnalyzer.featurewidth];
      cla = FishBatchClassifier(self.nfish,featdim,varargin{:});
    end


    function clprob = computeClassProb(self,idfeatures)
      % just calls the FishClassifer 
      clprob = predict(self.fishClassifier,idfeatures);
    end
    
    
    function feat = getFeatureDataFromTracks(self,trackidx)
      s = 0;
      for i = trackidx
        s = s+1;
        track = self.tracks(i);
        %take data
        dataidx = 1:track.nextBatchIdx-1;
        feat{s} = track.batchFeatures(dataidx,:);
      end
    end

    
    function updateTrackFishIds(self,trackIndex,newFishIds,varargin)
    % updates all ids since the last crossing

      oldFishIds = cat(2,self.tracks(trackIndex).fishId);
      
      if any(isnan(trackIndex)) || any(isnan(newFishIds)) 
        keyboard
      end
      
      if any(newFishIds~=oldFishIds)     

        if ~all(sort(newFishIds)==sort(oldFishIds)) || length(newFishIds)==1
          % not a permutation. 
          return;
        end

        f2t = self.fishId2TrackId;
        pos = self.pos;

        if nargin>3
          startframe = varargin{1};
        else
          startframe = [];
        end

        for i = 1:length(trackIndex)
          if ~isnan(newFishIds(i))
            oldFishId = oldFishIds(i);
            newFishId = newFishIds(i);

            if oldFishId~=newFishId
              verbose('switching fish %d->%d',oldFishId,newFishId)
              self.tracks(trackIndex(i)).fishId = newFishId;
              
              % update 
              t = self.currentFrame;
              if isempty(startframe)
                startframe = self.tracks(trackIndex(i)).frameOfLastCrossing;
              end
              
              trackId = self.tracks(trackIndex(i)).id;

              % swap
              f2t(startframe:t,newFishId) = self.fishId2TrackId(startframe:t,oldFishId);
              pos(:,newFishId,startframe:t) = self.pos(:,oldFishId,startframe:t);
            
              self.uniqueFishFrames = 0; % reset (if there is an update order is compromised) 
              self.resetBatchIdx(trackIndex(i)); % to make sure that features do not get mixed
            end
          end
        end
        self.fishId2TrackId = f2t;
        self.pos = pos;
      end
    end
    
    
    
    function resetBatchIdx(self,trackidx)
      [self.tracks(trackidx).nextBatchIdx] = deal(1);
    end
    
    
    
    function updateClassifier(self)
    % updates the identity classifier and tracks the identity to
    % trackID assisgments

     % save current track->fishID
      t = self.currentFrame;
      if size(self.fishId2TrackId,1)< t
        self.fishId2TrackId = cat(1,self.fishId2TrackId,nan(250,self.nfish));
      end
      
      fishIds = cat(2,self.tracks.fishId);      
      trackIds = cat(2,self.tracks.id);      
      self.fishId2TrackId(t,fishIds) = trackIds;
        
      % save track positions
      if size(self.pos,3)< t
        if isempty(self.pos)
          self.pos = nan(2,self.nfish,250);
        else
          self.pos = cat(3,self.pos,nan(2,self.nfish,250));
        end
      end
      self.pos(1:2,fishIds,t) = cat(1,self.tracks.centroid)';

      if isempty(trackIds) % no tracks
        return;
      end
      
      % check whether exist crossing or time to update
      crossif = any(self.crossings(:));
      
      if ~crossif && length(self.tracks)==self.nfish
        self.uniqueFishFrames = self.uniqueFishFrames + 1;
      else
        self.uniqueFishFrames = 0;
      end
      nfeat = cat(2,self.tracks.nextBatchIdx)-1;
      
      if  ~isInit(self.fishClassifier) && min(nfeat)< self.opts.classifier.minBatchN && self.uniqueFishFrames<self.opts.classifier.nFramesForInit
        %wait regardless of crossing
        return
      end

      if ~isInit(self.fishClassifier) 
        % First crossing but not yet initialized. Force to init! 
        feat = self.getFeatureDataFromTracks(1:length(self.tracks));
        [~,sidx] = sort(fishIds);
        feat = feat(sidx);
        batchsample = cell(1,self.nfish);
        [batchsample{:}] = deal([]);
        batchsample(1:length(feat)) = feat;
        
        self.fishClassifier.init(batchsample);
        newassignments = nan;
        self.resetBatchIdx(find(~isnan(newassignments)));
        return
      end
      
      
      if self.uniqueFishFrames-1 > self.opts.classifier.nFramesForUniqueUpdate  ...
        
        feat = self.getFeatureDataFromTracks(1:length(self.tracks));
        % uniqueFishFrames> thres
        % INITILIZE NEW CLASSIFIER!
        [newassignments prob] = batchUpdate(self.fishClassifier,fishIds,feat,1); % force?
    % $$$           verbose('Learn new classifiers!')
% $$$           [~,sidx] = sort(fishIds);
% $$$           feat = feat(sidx);
% $$$           batchsample = cell(1,self.nfish);
% $$$           [batchsample{:}] = deal([]);
% $$$           batchsample(1:length(feat)) = feat;
% $$$ 
% $$$           self.fishClassifier.init(batchsample);
% $$$           newassignments = nan;
% $$$           self.uniqueFishFrames = 0;

        self.resetBatchIdx(1:length(self.tracks));
        newassignments = NaN; 
      end

      % handle reseting after crosses (but currently no cross)
      % in this case the LastIds should already be handled (because currently no crossing, thus the last crossing is
      % in the current fields. 
      for i = 1:length(self.tracks)
        cross =self.crossings(i,:); 
        if length(self.tracks(i).crossedCurrentTrackIds)>0 ...
            &&  self.currentFrame-self.tracks(i).frameOfCurrentCrossing > self.opts.classifier.nFramesAfterCrossing ...
            && ~any(cross)

          %handle all current 
          crossedTracks = self.tracks(i).crossedCurrentTrackIds;
            
          % look if the tracks still exist
          trackidx= find(ismember(trackIds,crossedTracks));
          crossedTracks = trackIds(trackidx);
          [~,crossedFishIds] = ismember(crossedTracks,self.fishId2TrackId(self.tracks(i).frameOfCurrentCrossing,:));
            
          if any(~crossedFishIds)
            %some tracks got deleted...
            trackidx(~crossedFishIds)=[];
            crossedFishIds(~crossedFishIds) = [];
          end
            
          % test assignments
          feat = self.getFeatureDataFromTracks(trackidx);
          [assignedFishIds prob] = self.fishClassifier.predictPermutedAssignments(feat,crossedFishIds);
            
          if min(prob)>self.opts.classifier.reassignProbThres
            idx = ~isnan(assignedFishIds);
            self.updateTrackFishIds(trackidx(idx),assignedFishIds(idx),self.tracks(i).frameOfCurrentCrossing);
          end
          %self.resetBatchIdx(i);  do not reset. if updated get's reset anyways
          self.tracks(i).crossedCurrentTrackIds = []; % show that it was handled
        end
      end

      
      % handles crossings
      for i = 1:length(self.tracks)
        cross =self.crossings(i,:); 
        if any(cross) % currently at crossing. do handle the last crossing
            
        
          if length(self.tracks(i).crossedLastTrackIds)==1 
            % self crossing. max data reacjed Update
            feat = self.getFeatureDataFromTracks(i);
            batchUpdate(self.fishClassifier,self.tracks(i).fishId,feat);
            self.tracks(i).crossedLastTrackIds = []; % show that it was handled
            self.resetBatchIdx(i);
            
          elseif length(self.tracks(i).crossedLastTrackIds)>1 
            % multiple crossings. 
            
            %handle all the previous
            crossedTracks = self.tracks(i).crossedLastTrackIds;
            
            % look if the tracks still exist
            trackidx= find(ismember(trackIds,crossedTracks));
            crossedTracks = trackIds(trackidx);
            [~,crossedFishIds] = ismember(crossedTracks,self.fishId2TrackId(self.tracks(i).frameOfLastCrossing,:));
            
            if any(~crossedFishIds)
              %some tracks got deleted...
              trackidx(~crossedFishIds)=[];
              crossedFishIds(~crossedFishIds) = [];
            end
            
            % test assignments
            feat = self.getFeatureDataFromTracks(trackidx);
            [assignedFishIds prob] = self.fishClassifier.predictPermutedAssignments(feat,crossedFishIds);
            
            if min(prob)>self.opts.classifier.reassignProbThres
              idx = ~isnan(assignedFishIds);
              self.updateTrackFishIds(trackidx(idx),assignedFishIds(idx),self.tracks(i).frameOfLastCrossing);
            end
            self.resetBatchIdx(i);
            self.tracks(i).crossedLastTrackIds = []; % show that it was handled
            %self.tracks(i).frameOfLastCrossing = NaN; % DONOT RESET THIS BECAUSE OF LATER UPDATE
          else
            % seems to be the first crossing. Wait. 
          end
        end
      end
      

      if self.displayif>1
        if any(self.crossings(:))
          figure(2)
          self.fishClassifier.plotMeans();
        end
      end
      

    end
        
   
    %% TRACKS

    function initializeTracks(self)
    % create an empty array of tracks
      self.tracks = struct(...
        'id', {}, ...
        'fishId', {}, ...
        'bbox', {}, ...
        'centroid',{},...
        'features',{},...
        'segment', {},...
        'lastSegment', {},...
        'classProb',{},...
        'batchFeatures',{},...
        'nextBatchIdx',{},...
        'predictor', {}, ...
        'frameOfCurrentCrossing', {}, ...
        'crossedCurrentTrackIds',{},...
        'frameOfLastCrossing', {}, ...
        'crossedLastTrackIds',{},...
        'age', {}, ...
        'totalVisibleCount', {}, ...
        'consecutiveInvisibleCount', {},...
        'assignmentCost', {} );
    end
    
    
    function predictNewLocationsOfTracks(self)
      
      for i = 1:length(self.tracks)
        bbox = self.tracks(i).bbox;
        
        % Predict the current location of the track.
        predictedCentroid = predict(self.tracks(i).predictor);
        
        % Shift the bounding box so that its center is at
        % the predicted location.
        predictedCentroid(1:2) = int32(predictedCentroid(1:2)) - bbox(3:4) / 2;
        self.tracks(i).bbox = [predictedCentroid(1:2), bbox(3:4)];
        self.tracks(i).centroid = predictedCentroid;
      end
    end
      
    
    function fcost = computeCostMatrix(self)
      
      nTracks = length(self.tracks);
      nDetections = size(self.centroids, 1);
      fcost = zeros(nTracks, nDetections);
      features = self.features;
      idfeatures = self.idfeatures;
      centroids = self.centroids;
      tracks = self.tracks;
      medtau = self.opts.tracks.medtau;
      segments = self.segments;
      clprob = self.classProb;
      
      if isempty(self.X)
        self.X = cell(1,self.nfish);
        self.y = cell(1,self.nfish);
      end
      
      
      if nDetections>0
        %mdist = min(pdist(double(centroids)));
        
        % Compute the cost of assigning each detection to each track.
        for i = 1:nTracks
          cost = {};
          % center position 
          center = tracks(i).centroid; % kalman predict/pointtracker
          cost{end+1} = pdist2(double(center),double(centroids),'Euclidean');

          % cost of 

          n = zeros(1,length(segments));
          for jj = 1:length(segments)
            memb = ismember(tracks(i).segment.PixelIdxList(:)',segments(jj).PixelIdxList(:)');
            n(jj) = sum(memb)/length(memb);
          end
          cost{end+1} = 1-n + (n==0)*10;
          

          if ~isempty(tracks(i).classProb) && sum(tracks(i).classProb)>1e-5
            % id distance 
            cost{end+1} = pdist2(tracks(i).classProb,clprob,'correlation');
            cost{end}(isnan(cost{end})) = 1;
          else
            % no prob. just constant;
            cost{end+1} = ones(1,nDetections);
          end
          
          % other features
          for j = 1:size(features,2)
            cost{end+1}= pdist2(double(tracks(i).features(j)),double(features(:,j)),'Euclidean');
          end

        
          fcost(i,:) = 0;
          for k = 1:length(cost)
            fcost(i,:) = fcost(i,:) + cost{k}/self.medianCost(k)*self.scalecost(k); % scale cost
            self.medianCost(k) = self.medianCost(k)*(medtau-1)/medtau  + 1/medtau*median(cost{k}(:));
          end
       
        end
      end

    end
    
    
      
    function detectionToTrackAssignment(self)

      nDetections = size(self.centroids, 1);

      self.cost = [];
      self.assignments= [];
      self.unassignedTracks = [];
      self.unassignedDetections = [];
      
      if nDetections
        
        % determine cost of each assigment to tracks
        fcost = computeCostMatrix(self);
        
        % Solve the assignment problem.
        [self.assignments, self.unassignedTracks, self.unassignedDetections] = ...
            assignDetectionsToTracks(fcost, self.opts.tracks.costOfNonAssignment);
        
        outcost = Inf*zeros(1,size(fcost,2));
        if ~isempty(self.assignments)
          for i = 1:size(self.assignments,1)
            outcost(self.assignments(i,2)) = fcost(self.assignments(i,1),self.assignments(i,2));
          end
          for i = 1:length(self.unassignedDetections)
            outcost(self.unassignedDetections(i)) = min(fcost(:,self.unassignedDetections(i)));
          end
        end
        
        self.cost = outcost;
      end
    
    end
      
    
    
    function updateTracks(self)

      %% update assigned
      assignments = self.assignments;
      tracks = self.tracks;
      mdist = min(pdist(self.centroids));
      numAssignedTracks = size(assignments, 1);
      for i = 1:numAssignedTracks
        trackIdx = assignments(i, 1);
        detectionIdx = assignments(i, 2);
        centroid = self.centroids(detectionIdx, :);
        bbox = self.bboxes(detectionIdx, :);
        segm = self.segments(detectionIdx);
        feat = self.features(detectionIdx,:);
        cost = self.cost(detectionIdx);
        idfeat = self.idfeatures(detectionIdx,:);
        if ~isempty(self.classProb)
          classprob = self.classProb(detectionIdx,:);
        else
          classprob = [];
        end
        
        % Correct the estimate of the object's location
        % using the new detection.
        correct(tracks(trackIdx).predictor, [centroid]);
        
        % update classifier 
        % if isempty(tracks(trackIdx).classProb) || isempty(classprob)
        tracks(trackIdx).classProb =  classprob;
          %    tracks(trackIdx).nClassProbUpdates = 1;
          % else
          %  n = tracks(trackIdx).nClassProbUpdates;
          % tracks(trackIdx).classProb =  n/(n+1) * tracks(trackIdx).classProb  + 1/(n+1)*classprob;
          % tracks(trackIdx).nClassProbUpdates = tracks(trackIdx).nClassProbUpdates + 1;
          % end
        
        % Replace predicted bounding box with detected
        % bounding box.
        tracks(trackIdx).bbox = bbox;
        
        % Update track's age.
        tracks(trackIdx).age = tracks(trackIdx).age + 1;

        %segment/features
        tracks(trackIdx).lastSegment = tracks(trackIdx).segment;
        tracks(trackIdx).segment = segm;
        tracks(trackIdx).features = feat;

        if tracks(trackIdx).segment.bendingStdValue<self.opts.classifier.bendingThres
          tracks(trackIdx).batchFeatures(tracks(trackIdx).nextBatchIdx,:)  = idfeat;
          tracks(trackIdx).nextBatchIdx = tracks(trackIdx).nextBatchIdx+1;
        end
        
        tracks(trackIdx).centroid = centroid;
        tracks(trackIdx).assignmentCost = cost;

        % Update visibility.
        tracks(trackIdx).totalVisibleCount = ...
            tracks(trackIdx).totalVisibleCount + 1;
        tracks(trackIdx).consecutiveInvisibleCount = 0;
      end

      %% update unassigned
      for i = 1:length(self.unassignedTracks)
        ind = self.unassignedTracks(i);
        tracks(ind).age = tracks(ind).age + 1;
        tracks(ind).consecutiveInvisibleCount = ...
            tracks(ind).consecutiveInvisibleCount + 1;
      end

      
      self.tracks = tracks;
      
    end

    
    function deleteLostTracks(self)
      if isempty(self.tracks)
        return;
      end
      
      % Compute the fraction of the track's age for which it was visible.
      ages = [self.tracks(:).age];
      totalVisibleCounts = [self.tracks(:).totalVisibleCount];
      visibility = totalVisibleCounts ./ ages;
      
      % Find the indices of 'lost' tracks.
      lostInds = (ages < self.opts.tracks.ageThreshold & visibility < 0.6) | ...
          [self.tracks(:).consecutiveInvisibleCount] >= self.opts.tracks.invisibleForTooLong;
      
      % Delete lost tracks.
      self.tracks = self.tracks(~lostInds);
    end
    
    
    function createNewTracks(self)
    % create new tracks for unassigndetections (if less the number of fish)

      availablefishids = setdiff(1:self.nfish,cat(2,self.tracks.fishId));
      
      s =0;
      for i = self.unassignedDetections(:)'
        s = s+1;
        if s<=length(availablefishids)
          newfishid = availablefishids(s);
        else
          % just take the first detections. Maybe sort which to take?
          break
        end
        
        pred = self.newPredictor(self.centroids(i,:));

        %save the features to later update the batchClassifier
        idfeat = self.idfeatures(i,:);
        bfeat = zeros(self.opts.classifier.maxFramesPerBatch,size(idfeat,2));
        bfeat(1,:) = idfeat;
        if ~isempty(self.classProb)
          clprob = self.classProb(i,:);
        else
          clprob = [];
        end
        
        % Create a new track.
        newTrack = struct(...
          'id',        self.nextId,       ...
          'fishId',   newfishid,...
          'bbox',      self.bboxes(i,:),  ...
          'centroid',  self.centroids(i,:),  ...
          'features',  self.features(i,:),...
          'segment',   self.segments(i),  ... 
          'lastSegment',[],               ... 
          'classProb',clprob,...
          'batchFeatures',bfeat,          ...
          'nextBatchIdx', 1,              ...
          'predictor', pred,              ...
          'frameOfCurrentCrossing', NaN, ...
          'crossedCurrentTrackIds',[],...
          'frameOfLastCrossing', NaN, ...
          'crossedLastTrackIds',[],...
          'age',       1,                 ...
          'totalVisibleCount', 1,         ...
          'consecutiveInvisibleCount', 0, ...
          'assignmentCost',Inf);

        % Add it to the array of tracks.
        self.tracks(end + 1) = newTrack;
        
        % Increment the next id.
        self.nextId = self.nextId + 1;
        

      end
    end
    

    function detectTrackCrossings(self)
      
      self.crossings = false(length(self.tracks),length(self.tracks));
      if length(self.tracks)<2
        return
      end

      centroids = cat(1,self.tracks.centroid);
      dist = squareform(pdist(centroids));
      msk = dist<self.opts.classifier.maxCrossingDistance;
      
      for i = 1:length(self.tracks)
        track = self.tracks(i);
        dist(i,i) = Inf;

        if sum(msk(i,:))>1   % self-dist always smaller, thus >1
          % MAYBE ADD SOMETHING TO TEST WHETHER TRACKS ARE REALLY LOST
          
          if track.consecutiveInvisibleCount>0 || track.assignmentCost>self.opts.classifier.crossCostThres ...
            || any(dist(i,:)<self.opts.classifier.minCrossingDistance)
            self.crossings(i,:) = msk(i,:);
          end
          
        end

        % check whether storage of the feawtures is full and thus needs an update
        if track.nextBatchIdx > self.opts.classifier.maxFramesPerBatch;
          self.crossings(i,i) = true; %"self" crossing to update the classifier
        end

      end
      
      %make symmetric
      self.crossings = self.crossings | self.crossings';
    
    end
    
    
    function updateTrackCrossings(self)
    % detects and updated cross counter in tracks  

      self.detectTrackCrossings();
      trackIds = cat(2,self.tracks.id);

      for i = 1:length(self.tracks)
        cross = self.crossings(i,:); 

        if any(cross)

          % make sure that it is actually not the same crossing spanning multiple frames
          dt = self.currentFrame - self.tracks(i).frameOfCurrentCrossing;
          if dt<self.opts.classifier.minCrossingFrameDist
            % seems to be the same crossing. Might be others involved, but their tracks will get handled anyway
            %self.tracks(i).frameOfCurrentCrossing = self.currentFrame;
            self.tracks(i).crossedCurrentTrackIds =  union(self.tracks(i).crossedCurrentTrackIds,trackIds(find(cross)));
          else  
            %new crossing. 
            self.tracks(i).frameOfLastCrossing = self.tracks(i).frameOfCurrentCrossing;
            self.tracks(i).frameOfCurrentCrossing = self.currentFrame;
            
            self.tracks(i).crossedLastTrackIds = self.tracks(i).crossedCurrentTrackIds;
            self.tracks(i).crossedCurrentTrackIds = trackIds(find(cross));
          end
        end
      end
    end
    
    
    function handleTracks(self)
    % handles the assignment of the tracks etc
      
      self.predictNewLocationsOfTracks();
      self.detectionToTrackAssignment();
      
      self.updateTracks();
      %self.deleteLostTracks(); % ? MAYBE DO NOT DELETE
      
      if length(self.tracks)<self.nfish
        self.createNewTracks();
      end

      self.updateTrackCrossings();
    end
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  methods


    function self = FishTracker(vid,varargin) 
    % constructor
      self = self@handle();

      if ~exist('emclustering') && exist('helper')
        addpath('helper');
      end
      
      % some additional defaults [can be overwritten by eg 'detector.thres' name]
      self.opts(1).detector(1).thres = 'auto';
      self.opts.detector.dtau = 10; % 
      self.opts.detector.mtau= 1e4;
      self.opts.detector.inverse= 0;
      self.opts.detector.rgbchannel= [];
      
      % blob anaylser
      self.opts(1).blob(1).overlapthres= 0.8; % just for init str
      self.opts.blob.minArea = 100;
      
      % classifier 
      self.opts(1).classifier.crossCostThres = 3;
      self.opts(1).classifier.reassignProbThres = 0.2;
      self.opts(1).classifier.maxFramesPerBatch = 200; 
      self.opts(1).classifier.minBatchN = 25; 
      self.opts(1).classifier.nFramesForInit = 100; 
      self.opts(1).classifier.minCrossingDistance = [];  % will set based on fishlength;
      self.opts(1).classifier.maxCrossingDistance = [];  % will set based on fishlength
      self.opts(1).classifier.minCrossingFrameDist = 5;
      self.opts(1).classifier.nFramesAfterCrossing = 30;
      self.opts(1).classifier.nFramesForUniqueUpdate = 50; % all simultanously...
      self.opts(1).classifier.bendingThres = 1.5;
      % tracks
      self.opts(1).tracks.medtau = 30;
      self.opts(1).tracks.costOfNonAssignment =  2*mean(self.medianCost) + sum(self.scalecost);   

      %lost tracks
      self.opts(1).tracks.invisibleForTooLong = 5;
      self.opts(1).tracks.ageThreshold = 10;


      args = varargin;
      nargs = length(args);
      if nargs>0 && mod(nargs,2)
        error('expected arguments of the type ("PropName",pvalue)');
      end
      for i = 1:2:nargs
        if ~ischar(args{i})
          error('expected arguments of the type ("PropName",pvalue)');
        else
          idx = findstr(args{i},'.');
          if ~isempty(idx)
            optsname = args{i}(1:idx-1); 
            assert(any(strcmp(optsname,fieldnames(self.opts))))
            self.opts.(optsname).(args{i}(idx+1:end)) = args{i+1};
          else
            assert(~strcmp(args{i},'opts'));
            assert(isprop(self,args{i}))
            self.(args{i}) = args{i+1};
          end
        end
      end
      
      self.setupSystemObjects(vid);
      self.initializeTracks(); % Create an empty array of tracks.
    end

    
    function trackinfo = getCurrentTracks(self,allif)
      % gets the track info 
      
      trackinfo = [];
      trackinfo = self.tracks;
      if ~allif
        trackinfo = rmfield(trackinfo,{'segment','predictor','nextBatchIdx','batchFeatures','age','totalVisibleCount','consecutiveInvisibleCount','lastSegment'});
      end
      
    end
  
       
    function displayTrackingResults(self)

      cjet = jet(self.nextId);

      minVisibleCount = 2;
      tracks = self.tracks;
      uframe = self.uframe;
      if ~isempty(tracks)
        
        % Noisy detections tend to result in short-lived tracks.
        % Only display tracks that have been visible for more than
        % a minimum number of frames.
        reliableTrackInds = ...
            [tracks(:).totalVisibleCount] > minVisibleCount;
        reliableTracks = tracks(reliableTrackInds);
        
        % Display the objects. If an object has not been detected
        % in this frame, display its predicted bounding box.
        if ~isempty(reliableTracks)
          % Get bounding boxes.
          bboxes = cat(1, reliableTracks.bbox);
          
          % Get ids.
          ids = int32([reliableTracks(:).id]);

          % Create labels for objects indicating the ones for
          % which we display the predicted rather than the actual
          % location.
          labels = cellstr(int2str(ids'));
          predictedTrackInds = [reliableTracks(:).consecutiveInvisibleCount] > 0;
          isPredicted = cell(size(labels));
          isPredicted(predictedTrackInds) = {' predicted'};

          for i = 1:length(reliableTracks)
            if ~predictedTrackInds(i) & reliableTracks(i).assignmentCost>150
              isPredicted{i} = sprintf('  %1.0d',round(reliableTracks(i).assignmentCost));
            end
          end
          

          labels = strcat(labels, isPredicted);
          cols = jet(self.nfish);
          
          clprob = cat(1,reliableTracks.classProb);
          [cl,cli] = max(clprob,[],2);
          if length(cl)==length(reliableTracks)
            clabels = uint8(cols(cli,:)*255);

            % Draw the objects on the frame.
            uframe = insertObjectAnnotation(uframe, 'rectangle', bboxes, labels,'Color',clabels);
            
          else
            % Draw the objects on the frame.
            uframe = insertObjectAnnotation(uframe, 'rectangle', bboxes, labels,'Color', uint8([0.5,0.5,0.5]*255));
          end
          
          center = cat(1, reliableTracks.centroid);
          center = cat(2,center,max(self.fishlength/2,self.fishwidth)*ones(size(center,1),1));
          cross = find(cat(1,reliableTracks.frameOfCurrentCrossing)==self.currentFrame);
          uframe = insertObjectAnnotation(uframe, 'circle', center(cross,:), 'Crossing!','Color','red');
          
        
          %% insert some  features

          for i = 1:length(reliableTracks)

            id = reliableTracks(i).fishId;
            segm = reliableTracks(i).segment;

            
            if isfield(segm,'MSERregions') && ~isempty(segm.MSERregions)
              for j = 1:length(segm.MSERregions)
                pos = segm.MSERregionsOffset + segm.MSERregions(j).Location;
                uframe = insertMarker(uframe, pos, '*', 'color', uint8(cols(id,:)*255), 'size', 4);
              end
            end
            pos = reliableTracks(i).centroid;
            uframe = insertMarker(uframe, pos, 'o', 'color', uint8(cols(id,:)*255), 'size', 5);
          
            
            
            
            %% insert more markers
            howmany = 20;
            trackpos = self.pos(:,:,max(self.currentFrame-howmany,1):self.currentFrame);
            trackpos(:,:,any(any(isnan(trackpos),1),2)) = [];
            if ~isempty(trackpos)
              for ii = 1:self.nfish
                uframe = insertMarker(uframe, squeeze(trackpos(:,ii,:))', 'o', 'color',  uint8(cols(ii,:)*255), 'size', 2);
              end
            end
          end
        end


        self.stepVideoPlayer(uframe);

      end
      
    end

    
    
    function [pos,t,savedTracks] = track(self)
    % Detect objects, and track them across video frames.
      s = 0;

      self.setCurrentTime(self.timerange(1));
      
      verbose('Start Tracking...')
      t = [];
      while  hasFrame(self)
        s = s+1;
        self.currentFrame = s; 

        self.readFrame();

        self.detectObjects();

        self.handleTracks();


        self.updateClassifier();

        
        savedTracks(s).tracks = self.getCurrentTracks(0);
        t(s) = self.currentTime - self.dt; % current time is actually the next frame's time
        
        
        if self.displayif
          self.displayTrackingResults();
        else
          fprintf(' computed frame %d...\r',s);
        end
      
      end
      pos = permute(self.pos(:,:,1:self.currentFrame),[3,1,2]);
    end

    
    
    
    
  end
  
end


      

    