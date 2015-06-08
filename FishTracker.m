classdef FishTracker < handle;
  
  
  
  properties 
    
    scalecost = [20,5,1,0,1];
    opts = struct([]);
    saveFields = {'centroid','classProb','bbox','assignmentCost','consecutiveInvisibleCount', ...
                  'segment.Orientation','segment.bendingStdValue'};
    fishlength = [];
    fishwidth = [];
    maxVelocity = [];
    displayif = 0;
    timerange= [];
    nfish = 1;
    useGpu = 0;

  end

  properties (SetAccess = private)


    medianCost = ones(5,1);
    costinfo= {'Location','Pixel overlap','Classifier','MeanIntensity','Orientation'};
    detector = [];
    blobAnalyzer = [];
    fishClassifier = [];
    videoReader = [];
    videoPlayer = [];
    videoWriter = [];
    currentTime = 0;
    tracks = [];
    fishId2TrackId = [];
    visreader = 0;

    dt = 1;
    writefile = '';
    res = [];
    pos = [];
    
    videoFile = '';
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
    uniqueFishFrames = 0;
    currentFrame = [];
  end
  

  methods (Access=private)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % Gereral & init functions   
    %---------------------------------------------------

    function closeVideoObjects(self)
      if (self.visreader)
        release(self.videoReader); % close the input file
        if ~isempty(self.videoPlayer)
          release(self.videoPlayer);
        end
        
        if ~isempty(self.videoWriter)
          release(self.videoWriter);
        end
      else
        %close(self.videoReader); % close the input file
        if ~isempty(self.videoWriter)
          close(self.videoWriter);
        end
        if ~isempty(self.videoPlayer)
          release(self.videoPlayer);
        end
      end
    end
    
    function [writer] = newVideoWriter(self,vidname)
      
      if self.visreader
        writer = vision.VideoFileWriter(vidname,'FrameRate',1/self.dt);
      else
        writer = VideoWriter(vidname);
        open(writer);
      end
    end
    
    function [player] = newVideoPlayer(self,vidname)
      player = vision.VideoPlayer();
    end
    
    
    function resetVideoReader(self)
      if (self.visreader)
        release(self.videoReader); % close the input file
      else
        %close(self.videoReader); % close the input file
      end
      delete(self.videoReader);
      self.videoReader = newVideoReader(self,self.videoFile);
      self.setCurrentTime(self.currentTime);
    end
    
      
      
      
    function [reader, dt, timerange] = newVideoReader(self,vidname,timerange)
    % note: does not set the current time. Initialize with self.setCurrentTime
      if ~exist('timerange','var')
        timerange = [];
      end
      
      
      if ~isempty(timerange)  && self.visreader
        error('Timerange not supported with visreader');
      end
      
      if self.visreader
        reader = vision.VideoFileReader(vidname);
        isDone(reader); % check
      else
        reader = VideoReader(vidname);
      end

      if isempty(timerange)
        timerange = [0,Inf];
      end
      

      
      if ~self.visreader
        dt = 1/double(reader.FrameRate);
        assert(length(timerange)==2 && timerange(1)<reader.Duration/self.dt);

        if isinf(timerange(2))
          timerange(2) = reader.duration;
        end
        
      else
        tmp = info(reader);
        dt = 1/double(tmp.VideoFrameRate);
      
        if isinf(timerange(2))
          timerange(2) = tmp.duration;
        end

      end

      
    end

    function setDisplayType(self)

      if self.displayif>2
        self.fishClassifier.plotif = 1;
        self.blobAnalyzer.plotif = 1;
      end
      if ~self.displayif
        self.fishClassifier.plotif = 0;
        self.blobAnalyzer.plotif = 0;
      end
      
      if self.displayif> 1
        self.opts.tracks.displayEveryNFrame = 1;
      end
      
      if self.displayif && isempty(self.videoPlayer)
        self.videoPlayer = self.newVideoPlayer();
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
    
    function stepVideoWriter(self,uframe);
      if ~isempty(self.videoWriter)
        if self.visreader
          self.videoWriter.step(uframe);
        else
          self.videoWriter.writeVideo(uframe);
        end
      end
    end
    
    
    function blobAnalyzer=newBlobAnalyzer(self);
      % Connected groups of foreground pixels are likely to correspond to moving
      % objects.  The blob analysis system object is used to find such groups
      % (called 'blobs' or 'connected components'), and compute their
      % characteristics, such as area, centroid, and the bounding box.

      blobAnalyzer = MyBlobAnalysis('fishlength',self.fishlength,'fishwidth',self.fishwidth,...
                                    'maxextent',3*(self.fishlength+self.fishwidth),...
                                    'minextent',0.2*(self.fishlength+self.fishwidth));
        
    end
      
        
    function self=setupSystemObjects(self,vid)
      % Initialize Video I/O
      % Create objects for reading a video from a file, drawing the tracked
      % objects in each frame, and playing the video.

        
      %% Create a video file reader.
      [self.videoReader,self.dt,self.timerange] = self.newVideoReader(vid,self.timerange);
      self.videoFile = vid;
      
      % get some frames to initialize the detector
      n = max(self.opts.detector.nAutoFrames,1);
      verbose('Reading %d frames of the video to establish background',n);              
      timePoints = linspace(self.timerange(1),self.timerange(2),n+1);

      for i = 1:length(timePoints)-1
        self.setCurrentTime(timePoints(i));
          
        self.readFrame();
        if i==1
          autoFrames = zeros([size(self.frame),n],'like', self.frame);
        end
        autoFrames(:,:,:,i) = self.frame;
      end

      
      if ~isempty(self.writefile) && self.displayif
        self.videoWriter = self.newVideoWriter(self.writefile);
      else
        if ~isempty(self.writefile)
          warning('set displayif>0 for writing a video file!');
        end
        self.videoWriter = [];
      end
            
      %% get new detector
      args = {};
      for f = fieldnames(self.opts.detector)'
        if ~strcmp(f{1},'nAutoFrames')
          args{end+1} = f{1};
          args{end+1} = self.opts.detector.(f{1});
        end
      end
      self.detector = newForegroundDetector(self,args{:});


      [fishlength,fishwidth, fv, mu] = self.detector.init(autoFrames,self.nfish);

       
      if isempty(self.fishlength) 
        self.fishlength = ceil(1.2*mean(fishlength));
      end
      if isempty(self.fishwidth) 
        self.fishwidth = ceil(1.2*fishwidth);
      end
      
      self.fishwidth = max(self.fishwidth,8);
      assert(self.fishlength>self.fishwidth);
      
      
      
    end

    
    function initTracking(self)

    % MAYBE RESET THE CLASSIFIER ? 
                            
      self.setCurrentTime(self.timerange(1));
      
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

      self.setDisplayType();
      
      self.fishId2TrackId = nan(250,self.nfish);
      self.pos = [];
      self.res = [];
      self.currentFrame = 0;
      self.initializeTracks(); % Create an empty array of tracks.
      self.nextId = 1;
      self.uniqueFishFrames = 0;
      self.medianCost(:) = 1;
    
      self.uframe = [];
      self.frame = [];
    
      self.segments = [];
      self.features = [];
      self.bboxes = [];
      self.centroids = [];
      self.idfeatures = [];
      self.classProb = [];
      self.crossings = [];
      self.medianCost = ones(1,length(self.scalecost));
    
      self.assignments = [];
      self.unassignedTracks = [];
      self.unassignedDetections = [];
      self.cost = [];

      if isempty(self.maxVelocity)
        self.maxVelocity = self.fishlength*3/self.dt;
      end
    end
    
    
    
    function setCurrentTime(self,time)
      if self.visreader
        self.videoReader.reset();
      else
        self.videoReader.CurrentTime = time;
      end
      self.currentTime = time;
    end
    
    
    
    function readFrame(self);      
 
      if self.visreader
        if self.useGpu
          frame = gpuArray(self.videoReader.step());
        else
          frame = self.videoReader.step();
        end
        self.uframe = im2uint8(frame);
      else
        if self.useGpu
          frame = gpuArray(self.videoReader.readFrame());
        else
          frame = self.videoReader.readFrame();
        end
        self.uframe = frame;
        frame = im2single(frame);
      end
      self.frame = frame;
      self.currentTime = self.currentTime + self.dt;
      self.currentFrame = self.currentFrame + 1; 

      if ~mod(self.currentFrame,5e3)
        % memory leakage ? 
        self.resetVideoReader();
      end
    end
    
    function bool = hasFrame(self);      

      if self.currentTime>self.timerange(2) 
        bool = false;
      else
        
        if ~self.visreader
          bool = hasFrame(self.videoReader);
        else
          bool = ~isDone(self.videoReader);
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
  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ClASSIFICATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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

    
    function updateTrackFishIds(self,trackIndex,newFishIds,counterFieldName)
    % updates all ids since the last crossing

      oldFishIds = cat(2,self.tracks(trackIndex).fishId);
      
      if any(isnan(trackIndex)) || any(isnan(newFishIds)) 
        %should not happen
        warning('NAN detected!! This should not happen.')
        return
      end
      
      if any(newFishIds~=oldFishIds)     

        if ~all(sort(newFishIds)==sort(oldFishIds)) || length(newFishIds)==1
          % not a permutation. 
          return;
        end

        f2t = self.fishId2TrackId;
        pos = self.pos;

       
        for i = 1:length(trackIndex)
          if ~isnan(newFishIds(i))
            oldFishId = oldFishIds(i);
            newFishId = newFishIds(i);

            if oldFishId~=newFishId
              verbose('switching fish %d->%d\r',oldFishId,newFishId)
              self.tracks(trackIndex(i)).fishId = newFishId;
              
              % update 
              t = self.currentFrame;

              startframe = self.tracks(trackIndex(i)).(counterFieldName);
              
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

    
    function distangleMismatchedTracks(self,currentTrackIndices,type, updateif)

      if strcmp(type,'current')
        idFieldName = 'crossedCurrentTrackIds';
        counterFieldName = 'frameOfCurrentCrossing';
      else
        idFieldName =  'crossedLastTrackIds';
        counterFieldName = 'frameOfLastCrossing';
      end
      
      % handle all simultanously
      crossedTracks = [];
      for i = 1:length(currentTrackIndices)
        crossedTracks = union(crossedTracks,self.tracks(currentTrackIndices(i)).(idFieldName));
      end
      if isempty(crossedTracks)
        return;
      end

      trackIds = cat(2,self.tracks.id);      
      
      % look if the tracks still exist
      trackIndices= find(ismember(trackIds,crossedTracks));
      crossedTracks = trackIds(trackIndices);

      % take the current fishId of the existing tracks (should be the actual one)
      %frameIdx = self.tracks(currentTrackIndex).(counterFieldName);
      [~,crossedFishIds] = ismember(crossedTracks,self.fishId2TrackId(self.currentFrame,:));
            
      if any(~crossedFishIds)
        %seems that some tracks got deleted...
        trackIndices(~crossedFishIds)=[];
        crossedFishIds(~crossedFishIds) = [];
      end
      trackIndices = trackIndices(:)'; % make sure it is a rwo vector
       
      % test assignments
      feat = self.getFeatureDataFromTracks(trackIndices);
      [assignedFishIds prob] = self.fishClassifier.predictPermutedAssignments(feat,crossedFishIds);
            
      if min(prob)>self.opts.classifier.reassignProbThres
        idx = ~isnan(assignedFishIds);
        self.updateTrackFishIds(trackIndices(idx),assignedFishIds(idx),counterFieldName);
      end
      for i = 1:length(currentTrackIndices) 
        self.tracks(currentTrackIndices(i)).(idFieldName) = []; % show that it was handled
      end
      if updateif
        self.fishClassifierUpdate(currentTrackIndices);
      end

      if self.displayif>3
        figure(3);
        clf;
        plot(squeeze(self.pos(1,:,:))')
        drawnow;
      end
    end
    
      
    function fishClassifierUpdate(self,trackIndices)
    % updates the collected features according to the current fishIDs
      
      fishIds = cat(2,self.tracks(trackIndices).fishId);
      feat = self.getFeatureDataFromTracks(trackIndices);

      % this update might update difference fishIDs. .. maybe handle this an reassign. 
      [assignedFishIds prob] = batchUpdate(self.fishClassifier,fishIds,feat); % force?
      if min(prob)>self.opts.classifier.reassignProbThres
        idx = ~isnan(assignedFishIds);
        self.updateTrackFishIds(trackIndices(idx),assignedFishIds(idx),'frameOfLastCrossing');
      end
    end
    


    function updateClassifier(self)
    % updates the identity classifier and tracks the identity to
    % trackID assisgments
      
      if isempty(self.tracks) % no tracks
        return;
      end
      
     % save current track->fishID
      t = self.currentFrame;
      if size(self.fishId2TrackId,1)< t
        self.fishId2TrackId = cat(1,self.fishId2TrackId,nan(250,self.nfish));
      end
      
      fishIds = cat(2,self.tracks.fishId);      
      trackIds = cat(2,self.tracks.id);      
      self.fishId2TrackId(t,fishIds) = trackIds;
        
      % save track positions
      if isempty(self.pos)
        self.pos = nan(2,self.nfish,250);
      end
      
      if size(self.pos,3)<t
        self.pos = cat(3,self.pos,nan(2,self.nfish,250));
      end
      self.pos(1:2,fishIds,t) = cat(1,self.tracks.centroid)';

      % check whether exist crossing or time to update
      crossif = any(self.crossings(:));
      if ~crossif && length(self.tracks)==self.nfish
        self.uniqueFishFrames = self.uniqueFishFrames + 1;
      else
        self.uniqueFishFrames = 0;
        if ~isInit(self.fishClassifier) % we do not want a mixture at the beginning
          self.resetBatchIdx(1:length(self.tracks));
        end
      end
      nfeat = cat(2,self.tracks.nextBatchIdx)-1;
      
      if ~isInit(self.fishClassifier)
        if min(nfeat)> self.opts.classifier.minBatchN  && self.uniqueFishFrames>self.opts.classifier.nFramesForInit
          % classifier not yet initialized but a lot of unique frames. Force to init! Might get some crossings
          % as well but do not care. The  classifiers will get updated later
          
          feat = self.getFeatureDataFromTracks(1:length(self.tracks));
          [~,sidx] = sort(fishIds);
          feat = feat(sidx);
          batchsample = cell(1,self.nfish);
          [batchsample{:}] = deal([]);
          batchsample(1:length(feat)) = feat;
          
          self.fishClassifier.init(batchsample);
          self.resetBatchIdx(1:length(self.tracks));
          
          % reset all potentials previous crossings
          [self.tracks.frameOfLastCrossing] = deal(self.currentFrame);
          [self.tracks.frameOfCurrentCrossing] = deal(self.currentFrame);
          [self.tracks.crossedLastTrackIds] = deal([]);
          [self.tracks.crossedCurrentTrackIds] = deal([]);
          
          self.uniqueFishFrames  = 0;
        else
          return
        end
      end
      
      
      if self.uniqueFishFrames-1 > max(self.opts.classifier.nFramesForUniqueUpdate,self.opts.classifier.nFramesAfterCrossing)
        verbose('Perform uniqie frames update.\r')
        
        trackIndices = 1:self.nfish;
        self.fishClassifierUpdate(trackIndices);

        self.uniqueFishFrames = 0; % reset
        self.resetBatchIdx(trackIndices);
      end

      
      %% distangle wrongly fishIDs after crossings
      indicesCurrent = [];
      indicesLast = [];
      for trackIndex = 1:length(self.tracks)
       
        track = self.tracks(trackIndex);

        if ~any(self.crossings(trackIndex,:))
          % no crossing. But maybe the threshold for handlines the tracks after crossing is reached

          longAfterCrossing =  self.currentFrame - track.frameOfCurrentCrossing > self.opts.classifier.nFramesAfterCrossing;
          if ~isempty(track.crossedCurrentTrackIds)  &&  longAfterCrossing
            indicesCurrent = [indicesCurrent, trackIndex];
          end
        else   % currently some crossing !

          % indices the last crossing
          indicesLast = [indicesLast,trackIndex];
        end
      end

      if ~isempty(indicesLast)
        self.distangleMismatchedTracks(indicesLast,'last',false);
        % was some crossing. So better reset 
        self.resetBatchIdx(indicesLast);
      else % only indices current if there is otherwise no crossing (to avoid update errors)
        if ~isempty(indicesCurrent)
          self.distangleMismatchedTracks(indicesCurrent,'current',false);
          self.resetBatchIdx(indicesCurrent);
        end
      end
      

      
      %% unique update
      for trackIndex = 1:length(self.tracks)
        track = self.tracks(trackIndex);
        if ~any(self.crossings(trackIndex,:))
          longAfterCrossing =  self.currentFrame - track.frameOfCurrentCrossing > self.opts.classifier.nFramesAfterCrossing;
          if longAfterCrossing &&  (track.nextBatchIdx > self.opts.classifier.nFramesForSingleUpdate)
            % crossing already handled. FishIds *should* be correct
            verbose('Perform one fish update\r');
            self.fishClassifierUpdate(trackIndex);
            self.resetBatchIdx(trackIndex);
          end
        end
      end
      
      
      
      if self.displayif>2
        if any(self.crossings(:))
          figure(2)
          self.fishClassifier.plotMeans();
        end
      end
      

    end

    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  % Track related functions
  %---------------------------------------------------
 
  
    function initializeTracks(self)
    % create an empty array of tracks
      self.tracks = struct(...
        'id', {}, ...
        'fishId', {}, ...
        'bbox', {}, ...
        'centroid',{},...
        'features',{},...
        'segment', {},...
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
        % no prediction outside of the video:
        predictedCentroid(1) = max(min(predictedCentroid(1),size(self.frame,2)),1);
        predictedCentroid(2) = max(min(predictedCentroid(2),size(self.frame,1)),1);

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
      
      highCost = 50;
      somewhatCostly = 10;
      
      if nDetections>0
        %mdist = min(pdist(double(centroids)));
        
        % Compute the cost of assigning each detection to each track.
        for i = 1:nTracks
          cost = {};
          % center position 
          center = tracks(i).centroid; % kalman predict/pointtracker
          distmat = pdist2(double(center),double(centroids),'Euclidean');
          highidx = distmat>self.maxVelocity*self.dt; 
          distmat(highidx) = highCost*distmat(highidx);
          cost{end+1} = distmat;
          
          % cost of 
          n = zeros(1,length(segments));
          for jj = 1:length(segments)
            memb = ismember(tracks(i).segment.PixelIdxList(:)',segments(jj).PixelIdxList(:)');
            n(jj) = sum(memb)/length(memb);
          end
          cost{end+1} = 1-n + (n==0)*highCost;
          

          if ~all(isnan(tracks(i).classProb))
            % id distance 
            cost{end+1} = pdist2(tracks(i).classProb,clprob,'correlation');
            cost{end}(isnan(cost{end})) = 1;
            
            % make sure that at least one fishclass is reasonable probable
            msk = max(clprob,[],2)< self.opts.tracks.probThresForFish;
            cost{end}(msk) = highCost;
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

            % TRACKS ?
            mt = min(medtau*nTracks,self.currentFrame*nTracks);
            self.medianCost(k) = self.medianCost(k)*(mt-1)/mt  + 1/mt*median(cost{k}(:));
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
      
        % is there is currently a crossing then prob of non-assignement should be scaled down 
        costOfNonAssignment = self.opts.tracks.costOfNonAssignment;
        if ~isempty(self.tracks)
          frameOfCross = cat(1,self.tracks.frameOfCurrentCrossing);
          if ~isempty(frameOfCross) && self.currentFrame >  self.opts.classifier.minCrossingFrameDist;
            currrentlyCrossing = self.currentFrame - frameOfCross<=1;
            if any(currrentlyCrossing)
              costOfNonAssignment = costOfNonAssignment/(sum(currrentlyCrossing)+1);
            end
          end
        end
        
        % determine cost of each assigment to tracks
        fcost = computeCostMatrix(self);
        
        % Solve the assignment problem.
        [self.assignments, self.unassignedTracks, self.unassignedDetections] = ...
            assignDetectionsToTracks(fcost, costOfNonAssignment);
        
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
          classprob = nan(1,self.nfish);
        end
        
        % Correct the estimate of the object's location
        % using the new detection.
        if self.opts.tracks.kalmanFilterPredcition
          correct(tracks(trackIdx).predictor, [centroid]);
        end
        
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
        tracks(trackIdx).segment = segm;
        tracks(trackIdx).features = feat;

        if tracks(trackIdx).segment.bendingStdValue<self.opts.classifier.bendingThres
          tracks(trackIdx).batchFeatures(tracks(trackIdx).nextBatchIdx,:)  = idfeat;
          tracks(trackIdx).nextBatchIdx = min(tracks(trackIdx).nextBatchIdx+1,self.opts.classifier.maxFramesPerBatch);
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
        
        if self.opts.tracks.kalmanFilterPredcition
          pred = self.newPredictor(self.centroids(i,:));
        else
          pred = [];
        end
        
        
        %save the features to later update the batchClassifier
        idfeat = self.idfeatures(i,:);
        bfeat = zeros(self.opts.classifier.maxFramesPerBatch,size(idfeat,2));
        bfeat(1,:) = idfeat;
        if ~all(isnan(self.classProb))
          clprob = self.classProb(i,:);
        else
          clprob = nan(1,self.nfish);
        end
        
        % Create a new track.
        newTrack = struct(...
          'id',        self.nextId,       ...
          'fishId',   newfishid,...
          'bbox',      self.bboxes(i,:),  ...
          'centroid',  self.centroids(i,:),  ...
          'features',  self.features(i,:),...
          'segment',   self.segments(i),  ... 
          'classProb',clprob,...
          'batchFeatures',bfeat,          ...
          'nextBatchIdx', 1,              ...
          'predictor', pred,              ...
          'frameOfCurrentCrossing', self.currentFrame, ...
          'crossedCurrentTrackIds',[],...
          'frameOfLastCrossing', self.currentFrame, ...
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
      bbox = double(cat(1,self.tracks.bbox));
      centroids = cat(1,self.tracks.centroid);
      dist = squareform(pdist(centroids));
      bbox34 = bbox(:,3:4);
      msk = dist<max(bbox34(:));

      bBoxXOverlap = bsxfun(@ge,bbox(:,1) + bbox(:,3), bbox(:,1)') & bsxfun(@lt,bbox(:,1), bbox(:,1)');
      bBoxYOverlap = bsxfun(@ge,bbox(:,2) + bbox(:,4), bbox(:,2)') & bsxfun(@lt,bbox(:,2), bbox(:,2)');
      bBoxOverlap = bBoxXOverlap & bBoxYOverlap;
      
      for i = 1:length(self.tracks)
        track = self.tracks(i);

        if sum(msk(i,:))>1  
          
          if track.consecutiveInvisibleCount>0 || track.assignmentCost>self.opts.classifier.crossCostThres ...
            || sum(bBoxOverlap(i,:))>1
            self.crossings(i,:) = msk(i,:);
          end
          
        end

      end
      
      %make symmetric CATUTION MIGHT NOT BE SYMMETRIC!
      self.crossings = self.crossings | self.crossings';
      self.crossings(find(eye(length(self.crossings)))) = any(self.crossings,1);
    
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
            % do nothing
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
      
      if self.opts.tracks.kalmanFilterPredcition
        self.predictNewLocationsOfTracks();
      end
      
      self.detectionToTrackAssignment();
      
      self.updateTracks();
      %self.deleteLostTracks(); % ? MAYBE DO NOT DELETE
      
      if length(self.tracks)<self.nfish
        self.createNewTracks();
      end

      self.updateTrackCrossings();
    end
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % Results & Plotting helper functions
    %--------------------
    
    function generateResults(self,savedTracks)
      
      if isempty(savedTracks)
        return;
      end
      self.res.tracks = savedTracks;

      if size(savedTracks,2)<self.nfish
        self.res.tracks(end+1,1:self.nfish) =savedTracks(1);
        self.res.tracks(end,:) = [];
      end
      
      nFrames = self.currentFrame;

      self.res.pos = permute(self.pos(:,:,1:nFrames),[3,1,2]);      
      f2t = self.fishId2TrackId(1:nFrames,:);

      msk = reshape(arrayfun(@(x)~isempty(x.id),savedTracks),size(savedTracks));
      
      trackIdMat = zeros(nFrames,self.nfish);
      trackIdMat(msk) = cat(1,savedTracks.id);
      %fishIdMat = zeros(nFrames,self.nfish);

      % HOW CAN THIS DONE A BIT MORE EFFICIENTLY?
      order = 1:self.nfish;
      for i = 1:nFrames
        [~,loc] = ismember(trackIdMat(i,:),f2t(i,:));
        loc(~loc) = order(~loc);
        self.res.tracks(i,:) =self.res.tracks(i,loc);
      end

      for j = 1:self.nfish
        [self.res.tracks(:,j).fishId] = deal(j);
      end

    end
    
    
    function trackinfo = getCurrentTracks(self)
    % gets the track info 
      
      trackinfo = repmat(struct,size(self.tracks));

      % need id 
      [trackinfo(:).id] = deal(self.tracks.id);
      [trackinfo(:).t] = deal(self.currentTime - self.dt);

      % maybe a bit too slow?
      for i_f = 1:length(self.saveFields)
        f = self.saveFields{i_f};
        if any(f=='.')
          if isempty(self.tracks) 
            eval(sprintf('[trackinfo(:).%s] = deal([]);',f(1:find(f=='.',1,'first')-1)));
          else
            for i = 1:length(self.tracks)
              eval(sprintf('trackinfo(i).%s = self.tracks(i).%s;',f,f));
            end
          end
        else
          [trackinfo(:).(f)] = deal(self.tracks.(f));
        end
      end
    end
    
    
    function displayTrackingResults(self)

      cjet = jet(self.nextId);

      minVisibleCount = 2;
      tracks = self.tracks;
      uframe = gather(self.uframe);
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
          cl(isnan(cl)) = [];
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
            howmany = 10;
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
        self.stepVideoWriter(uframe);
          
      end
      
    end

  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  PUBLIC METHODS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contructor
    %--------------------
    function self = FishTracker(vid,varargin) 
    % constructor
      self = self@handle();

      if ~exist('emclustering') && exist('helper')
        addpath('helper');
      end
      
      % some additional defaults [can be overwritten by eg 'detector.thres' name]
      self.opts(1).detector(1).thres = 'auto';
      self.opts.detector.dtau = 0; % set to some value ONLY if the frame rate is very high!
      self.opts.detector.mtau= 1e3;
      self.opts.detector.inverse= 0;
      self.opts.detector.rgbchannel= [];
      self.opts.detector.nAutoFrames = 9;
      self.opts.detector.excludeBorderPercentForAutoThres = 0.05;
      
      % blob anaylser
      self.opts(1).blob(1).overlapthres= 0.8; % just for init str
      self.opts.blob.minArea = 100;
      self.opts.blob.colorfeature = true;
      self.opts.blob.interpif = 1;
      
      % classifier 
      self.opts(1).classifier.crossCostThres = 5;
      self.opts(1).classifier.reassignProbThres = 0.1;
      self.opts(1).classifier.maxFramesPerBatch = 150; 
      self.opts(1).classifier.minBatchN = 20; 
      self.opts(1).classifier.nFramesForInit = 30; % unique frames needed for init classifier
      self.opts(1).classifier.minCrossingFrameDist = 5;
      self.opts(1).classifier.nFramesAfterCrossing = 25;
      self.opts(1).classifier.nFramesForUniqueUpdate = 50; % all simultanously...
      self.opts(1).classifier.nFramesForSingleUpdate = 100; % single. should be larger than unique...
      self.opts(1).classifier.bendingThres = 1.5;
      % tracks
      self.opts(1).tracks.medtau = 100;
      self.opts(1).tracks.costOfNonAssignment =  numel(self.medianCost) + sum(self.scalecost);   
      self.opts(1).tracks.probThresForFish = 0.2;
      self.opts(1).tracks.displayEveryNFrame = 8;
      self.opts(1).tracks.kalmanFilterPredcition = 0; % better without

      %lost tracks
      self.opts(1).tracks.invisibleForTooLong = 5;
      self.opts(1).tracks.ageThreshold = 5;
      
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
    
      if 1/self.dt<25 && self.detector.dtau  % high framrate
        warning(['DTAU setting might only be suitable for high framne rate. \nSuggesting to set dtau to zero to get ' ...
                 'superior detections.']);
      end
    end
 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Main tracking loop
    %--------------------
  
    function track(self,continueif)
    % Detect objects, and track them across video frames.

      if isempty(self.res) || ~exist('continueif','var') || ~continueif 

        self.initTracking();
        verbose('Start tracking of %d fish  in the time range [%1.0f %1.0f]...',...
                self.nfish,self.timerange(1),self.timerange(2));
      else
        verbose('Continue tracking ...');
      end

      while  hasFrame(self)

        self.readFrame();

        self.detectObjects();
        
        self.handleTracks();

        self.updateClassifier();

        savedTracks(self.currentFrame,1:length(self.tracks)) = self.getCurrentTracks();

        % SOMETIMES SOME TRACKS ARE EMPTY. WHY ?
        DEBUG = 0;
        if DEBUG
          if any(cellfun('isempty',{self.tracks.id}))
            warning('Detected empty tracks')
            keyboard
          end
        end
        
        
        if self.displayif
          if ~mod(self.currentFrame-1,self.opts.tracks.displayEveryNFrame)
            self.displayTrackingResults();
          end
        else
          t = datevec(seconds(self.currentTime));
          if ~mod(self.currentFrame,10)
            verbose('Currently at time %1.0fh %1.0fm %1.1fs                      \r',t(4),t(5),t(6));
          end
        end
      end
      %% CAUTION handle last crossing 
      
      % make correct output structure
      self.generateResults(savedTracks);
      
      self.closeVideoObjects();
    end
    
    
    function [combinedFT varargout] = trackParallel(self,tOverlap,minDurationPerWorker)
    %  FTNEW = FT.TRACKPARALLEL() will track the results in parallel and creates a new FT object with tracks combined
    % 
    %  FT.TRACKPARALLEL(..,TOVERLAP,MINDURATIONPERWORKER) sepecifies the overlap between workers in
    %  seconds and the minmal video time per worker in seconds
    % 
    % CAUTION: if not calculated in parallel (w.g. for too short data) the returned handle object
    % will be the IDENTICAL handle (only a reference to the same data).  Only in case of parallel
    % processing the returned object will have a NEW handle (and thus reference new data).

      if ~exist('minDurationPerWorker','var')
        minDurationPerWorker = 300; % 5 minutes;
      end
      
      if ~exist('tOverlap','var') || isempty(tOverlap)
        % allow enough time for the classifer to init in the overlapping period
        tOverlap = 5*self.dt*max([self.opts.classifier.nFramesForInit,self.opts.classifier.nFramesForUniqueUpdate]); 
        tOverlap = max(tOverlap,self.opts.tracks.medtau*self.dt);
      end
      
      assert(minDurationPerWorker>tOverlap);
      totalDuration = diff(self.timerange);
      
      if totalDuration<2*tOverlap || totalDuration < 2*(minDurationPerWorker-tOverlap)
        % just on a single computer;
        self.track();
        combinedFT = self; % SHOULD IMPLEMENT COPY OBJECT HERE !!!
        varargout{1} = [];
        return
      end
      pool = gcp;
      
      maxDuration = min(diff(self.timerange)/pool.NumWorkers,1.2e4*self.dt);
      nparts = floor((diff(self.timerange)/maxDuration));
      
      sec = linspace(self.timerange(1),self.timerange(2),nparts+1);
      secEnd = sec(2:end);
      secStart = max(sec(1:end-1)-tOverlap,sec(1));
      timeRanges = [secStart(:),secEnd(:)];
      
      nRanges = size(timeRanges,1);
      
      % copy objects
      verbose('Parallel tracking from %1.1f to %1.1fh on %d workers',...
              self.timerange(1)/3600,self.timerange(2)/3600,pool.NumWorkers);
      verbose('%d processes with %1.1fm duration and %1.1fsec overlap.',...
              nRanges,max(diff(timeRanges,[],2))/60,tOverlap);

      % CAUTION THIS YIELDS A REFERENCE COPY ONLY !!! HOWEVER, RUNNING THIS WITH PARFOR LOOP WILL CHANGE
      % EACH HANDLE TO A VALUE BASED OBJECT ON EACH WORKER SO WE DO NOT NEED A FORMAL COPY METHOD. SEEMS TO
      % BE A HACK IN MATLAB.
      FTs = repmat(self,[1,nRanges]);
      res = {};
      parfor i = 1:nRanges
        ft = FTs(i);
        ft.timerange = timeRanges(i,:);
        ft.currentTime = ft.timerange(1);
        ft.resetVideoReader();
        ft.setCurrentTime(ft.timerange(1)); % to make sure 
        ft.displayif = 0;
        ft.track();
        res{i} = ft;
      end
      combinedFT = combine(res{:});
      if nargout>1
        varargout{1} = res;
      end
    end
    
  
  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Combination of two (or more) tracking objects with overlap
    %--------------------
    function combinedObj = combine(self,varargin)

      % make sure that the objects do not refer to the same data (i.e. have different handle
      % references)
      if any(self==cat(1,varargin{:}))
        error('cannot combine handles with identical data references');
      end
        
      
      tbinsize = max(self.fishlength/self.maxVelocity,2*self.dt);

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
        assert(strcmp(combinedObj.videoReader.Name,obj.videoReader.Name));
        assert(combinedObj.nfish==obj.nfish);

        res = getTrackingResults(obj);

        if obj.timerange(1)> combinedObj.timerange(2)
          % no overlap
          warning('no overlap. Fish IDs might get mixed up!!');
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
            warning(['no valid indices found for overlap. Fish IDs might get mixed up!!. Take all ' ...
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
          [X,Y] = ndgrid(tidxCombined,1:self.nfish*dim);
          overlappedCombinedPosInterp = reshape(accumarray([X(:),Y(:)],overlappedCombinedPos(:),[],@mean),[],obj.nfish);
          
          tidx = min(max(floor((tObj(overlapedIdx) - tstart)/tbinsize)+1,1),tidxCombined(end));
          [X,Y] = ndgrid(tidx,1:self.nfish*dim);
          overlappedPosInterp = reshape(accumarray([X(:),Y(:)],overlappedPos(:),[tidxCombined(end),self.nfish*dim],@mean),[],obj.nfish);

          % now the two position vectors can be compared. 
          cost = pdist2(overlappedCombinedPosInterp',overlappedPosInterp','correlation');

          % use the hungarian matching assignments from the vision toolbox
          assignments = assignDetectionsToTracks(cost,1e4);
          [~,sidx] = sort(assignments(:,1),'ascend');
          assignments = assignments(sidx,:); % make sure that the (combined) tracks are in squential order
          % append;
          combinedRes.pos = cat(1,combinedRes.pos,res.pos(overlapedIdx(end)+1:end,:,assignments(:,2)));
          combinedRes.tracks = cat(1,combinedRes.tracks,res.tracks(overlapedIdx(end)+1:end,assignments(:,2)));
        end

        combinedObj.res = combinedRes;
        combinedObj.timerange= [combinedObj.timerange(1),obj.timerange(2)];
        combinedObj.currentFrame = 1;
        combinedObj.currentTime = self.timerange(1);
        combinedObj.pos = []; % not valid 
        combinedObj.fishId2TrackId = []; % not valid 
      end
        
    end
    
    
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Results & plotting 
    %--------------------

    function res = getTrackingResults(self)
    % returns the current results
      if isempty(self.res)
        error('No results available. First track()...');
      else
        res = self.res;
      
        % delete beyond border pixels
        posx = squeeze(self.res.pos(:,1,:));
        posy = squeeze(self.res.pos(:,2,:));

        posx(posx>size(self.frame,2) | posx<1) = NaN;
        posy(posy>size(self.frame,1) | posy<1) = NaN;

        res.pos(:,1,:) = posx;
        res.pos(:,2,:) = posy;

        % cut those time points with too high assignemnt cost
        % WHY DOES THIS HAPPEN?? THERE SHOULDN'T BE ANY EMPTY VALUES!!
        emptymsk = find(cellfun('isempty',{res.tracks.assignmentCost}));
        [res.tracks(emptymsk).assignmentCost] = deal(Inf);

        assignmentCost = reshape(cat(1,res.tracks.assignmentCost),size(res.tracks));
        assignmentThres = quantile(assignmentCost(:),0.999);
        idx = find(assignmentCost > assignmentThres);

        for i = 1:size(res.pos,2)
          posi = res.pos(:,i,:);
          posi(idx) = NaN;
          res.pos(:,i,:) = posi;
        end
      end
    end
    
    function plotVideo(self,timerange,delay)

      res = self.getTrackingResults();      
      
      if ~exist('timerange','var')
        timerange = [res.tracks(1).t,res.tracks(end).t];
      end
      
      if ~exist('delay','var')
        delay = 0;
      end
      
      self.videoPlayer = self.newVideoPlayer();
      self.resetVideoReader();
      tr = self.timerange;
      self.timerange = timerange;
      self.setCurrentTime(timerange(1));
      self.currentFrame = 0;
      t_tracks =cat(1,res.tracks(:,1).t);
      tidx = find(t_tracks>=timerange(1) & t_tracks<timerange(2)); 
      t_tracks = t_tracks(tidx);
      selected_tracks = res.tracks(tidx,:);   
      cols = uint8(255*jet(self.nfish));
      s = 0;
      while self.hasFrame() & s<size(selected_tracks,1)
        self.readFrame();
        s = s+1;
        
        t = self.currentTime-self.dt;
        
        assert(abs(t_tracks(s)-t)<self.dt);

        uframe = self.uframe;
        tracks = selected_tracks(s,:);

        % Get bounding boxes.
        bboxes = cat(1, tracks.bbox);
          
        % Get ids.
        ids = int32([tracks(:).fishId]);
        labels = cellstr(int2str(ids'));
        clabels = cols(ids,:);

        highcost = isinf(cat(1,tracks.assignmentCost));
        clabels(highcost,:) = 127; % grey
        % Draw the objects on the frame.
        uframe = insertObjectAnnotation(uframe, 'rectangle', bboxes, labels,'Color',clabels);
        

        clprob = cat(1,tracks.classProb);
        clprob(isnan(clprob)) = 0;
        dia = self.fishlength/self.nfish;
        for i_cl = 1:self.nfish
          for j = 1:length(tracks)
            pos = [bboxes(j,1)+dia*(i_cl-1)+dia/2,bboxes(j,2)+ dia/2 ...
                   + bboxes(j,4),dia/2];
            uframe = insertShape(uframe, 'FilledCircle',pos, 'color', cols(i_cl,:), 'opacity',clprob(j,i_cl));
          end
        end
        self.stepVideoPlayer(uframe);
        d = seconds(t);
        d.Format = 'hh:mm:ss';
        verbose('CurrentTime %s\r',char(d))
        pause(delay);
        
      end
      
      self.timerange = tr;
      self.closeVideoObjects();
      
      
    end
    
    
    function plotTrace(self,varargin)
      self.plotByType('TRACE',varargin{:});
    end
    
    function plotProbMap(self,varargin)
      self.plotByType('PROBMAP',varargin{:});
    end
    
    function plotDomains(self,varargin)
      self.plotByType('DOMAINS',varargin{:});
    end

    function plotVelocity(self,varargin)
      self.plotByType('VELOCITY',varargin{:});
    end

    function plot(self,varargin)

      clf;
      subplot(2,2,1)
      self.plotTrace(varargin{:});

      subplot(2,2,2)
      self.plotVelocity(varargin{:});
      
      subplot(2,2,3)
      self.plotProbMap(varargin{:});
      
      subplot(2,2,4)
      self.plotDomains(varargin{:});

    end

    function plotByType(self,plottype,fishIds,plotTimeRange)
    % PLOTTRACE(FISHIDS) plots the tracking results
    % for the given FISHIDS

      if ~exist('plottype','var') || isempty(plottype)
        plottype = 'TRACE';
      end
      
      if isempty(self.res)
        warning('No results available. First track()...');
        return;
      end
    
      if ~exist('fishIds','var') || isempty(fishIds)
        fishIds = 1:self.nfish;
      end

      if ~exist('plotTimeRange','var') || isempty(plotTimeRange)
        plotTimeRange = self.timerange;
      end

      cla;
      res = self.getTrackingResults();
      
      t = cat(1,res.tracks(:,1).t);
      plotidx = t>=plotTimeRange(1) & t<plotTimeRange(2);
      
      centroid = reshape(cat(1,res.tracks.centroid),size(res.tracks,1),self.nfish,2);
      centroidx = centroid(plotidx,fishIds,1);
      centroidy = centroid(plotidx,fishIds,2);
      
      posx = squeeze(res.pos(plotidx,1,fishIds));
      posy = squeeze(res.pos(plotidx,2,fishIds));
      
      
      switch upper(plottype)

        case 'TRACE'
          plot(posx,posy);
          title('Traces');
          xlabel('X-Position [px]')
          ylabel('Y-Position [px]')
          
          msk = isnan(posx);
          hold on;
          plot(centroidx(msk),centroidy(msk),'.r','linewidth',1);
        
        case 'VELOCITY'
      

          plot(diff(posx)/self.dt,diff(posy)/self.dt,'.');
          xlabel('X-Velocity [px/sec]')
          ylabel('Y-Velocity [px/sec]')
          
          title('Velocity');
        
        case {'PROBMAP','DOMAINS'}


          dxy = 10;
          szFrame = size(self.frame);
          sz = ceil(szFrame/dxy);
          sz = sz(1:2);
          
          dposx = min(max(floor(posx/dxy)+1,1),sz(2));
          dposy = min(max(floor(posy/dxy)+1,1),sz(1));
          P = zeros([sz,length(fishIds)]);
          for i = 1:size(posx,2)
            msk = ~(isnan(dposy(:,i)) | isnan(dposx(:,i)));
            P(:,:,i) = accumarray([dposy(msk,i),dposx(msk,i)],1,sz);
            P(:,:,i) = P(:,:,i)/sum(sum(P(:,:,i)));
          end
          P(1,1,:) = 0;
          
          if strcmp(upper(plottype),'PROBMAP')
            Z = sum(P,3);
            Z = imfilter(Z,fspecial('disk',3),'same');
            imagesc(1:szFrame(2),1:szFrame(1),Z,[0,quantile(P(:),0.99)]);
            title('Overall probability');
            xlabel('X-Position [px]')
            ylabel('Y-Position [px]')
            axis xy;
          else

            [~,cl] = max(P,[],3);
            cl = cl + (sum(P,3)~=0);
            imagesc(1:szFrame(2),1:szFrame(1),cl);
            xlabel('X-Position [px]')
            ylabel('Y-Position [px]')
            title('Domain of fish')
            axis xy;
          end
          
        otherwise
          error('do not know plot type');
      end
    end
   
  end
  
  
end


      

    