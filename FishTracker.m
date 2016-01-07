classdef FishTracker < handle;
  

  
  properties 


    timerange= [];    

    opts = struct([]);
    saveFields = {'centroid','classProb','bbox','assignmentCost', ...
                  'consequtiveInvisibleCount', ...
                  'segment.Orientation','centerLine', ...
                  'thickness','segment.MinorAxisLength',...
                  'segment.MajorAxisLength','segment.reversed',...
                  'stmInfo'};%,...              'segment.FishFeatureC','segment.FishFeature','segment.FishFeatureCRemap'};
    stimulusPresenter = [];


  end

  properties (SetAccess = private)

    maxVelocity = [];
    displayif = 3;

    useMex = 1;
    useOpenCV = 1;
    useScaledFormat = 0;
    useKNN = 0;
    stmif = 0;
    
    costinfo= {'Location','Overlap','CenterLine','Classifier','Size','Area','BoundingBox'};
    scalecost = [10,3,3,2,1,2,1];


    nfish = [];
    fishlength = [];
    fishwidth = [];

    videoFile = '';

    uniqueFishFrames = 0;
    currentFrame = 0;
    currentCrossBoxLengthScale = 1; 


    meanAssignmentCost =1;
    tracks = [];
    cost = [];
    res = [];
    savedTracks = [];
    segments = [];
    
    videoHandler = [];
    videoPlayer = [];
    videoWriter = [];
    
    fishClassifier = [];

    writefile = '';
  end
  
  properties (SetAccess = private, GetAccess=private);

    saveFieldSegIf = [];
    saveFieldSub = {};

    nextId = 1; % ID of the next track
    isInitClassifier = [];

    centerLine = [];
    thickness = [];
    bboxes = [];
    centroids = [];
    idfeatures = [];
    features = [];

    
    classProb = [];
    classProbNoise = [];
    crossings = [];

    fishId2TrackId = [];
    meanCost = ones(4,1);
    pos = [];
    
    assignments = [];
    unassignedTracks = [];
    unassignedDetections = [];
    featurecosttypes = {'Area','Size','BoundingBox'};


  end
  

  methods (Access=private)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % Gereral & init functions   
    %---------------------------------------------------

    function closeVideoObjects(self)

      if ~isempty(self.videoWriter)
        % close writer
        if ~hasOpenCV() || ~self.useOpenCV
          relaese(self.videoWriter);
        else
          v = self.videoWriter;
          clear v;
        end
        self.videoWriter = [];
        self.writefile = '';
      end
      
    end
    
    function stepVideoWriter(self,frame);
      if ~isempty(self.videoWriter)
        if hasOpenCV() && self.useOpenCV
          self.videoWriter.write(frame);
        else
          self.videoWriter.step(frame);
        end
      end
    end
    
    
    function [writer] = newVideoWriter(self,vidname)
      verbose('Init VideoWriter "%s"',vidname);
      if hasOpenCV() && self.useOpenCV
        writer = cv.VideoWriter(vidname,self.videoHandler.frameSize([2,1]),'FPS',self.videoHandler.frameRate,'fourcc','X264');
      else
        writer = vision.VideoFileWriter(vidname,'FrameRate',self.videoReader.frameRate);
      end
    end
    
    
    function [player] = newVideoPlayer(self,vidname)
      player = vision.VideoPlayer();
      %player = FishVideoPlayer();
    end
    
    
    function [handler, timerange] = newVideoHandler(self,vidname,timerange,opts)
    % note: does not set the current time. 
      if ~exist('timerange','var')
        timerange = [];
      end
      if ~exist('opts','var')
        opts = self.opts;
      end
      if self.useMex && exist(['FishVideoCapture_.' mexext])
        handler = FishVideoHandlerMex(vidname,timerange,self.useKNN,opts);
      elseif hasOpenCV() && self.useOpenCV
        handler = FishVideoHandler(vidname,timerange,self.useKNN,opts);
      else
        handler = FishVideoHandlerMatlab(vidname,timerange,self.useKNN,opts);
      end
      timerange = handler.timeRange;

    end

    
    function setDisplayType(self)

      
      if self.displayif>3
        self.fishClassifier.plotif = 1;
      end

      if self.displayif>6
        self.videoHandler.plotting(true);
      else
        self.videoHandler.plotting(false);
      end
      
      if self.displayif==3
        self.opts.tracks.displayEveryNFrame = 20;
      end

      if self.displayif==2
        self.opts.tracks.displayEveryNFrame = 50;
      end

      if self.displayif==1
        self.opts.tracks.displayEveryNFrame = 100;
      end
      
      if ~self.displayif
        self.fishClassifier.plotif = 0;
        self.videoHandler.plotting(false);
      end
      
      if self.displayif> 5
        self.opts.tracks.displayEveryNFrame = 1;
      end
      
      if self.displayif && isempty(self.videoPlayer)
        self.videoPlayer = self.newVideoPlayer();
      end
      
    end
    
    
    function predictor = newPredictor(self,centroid);
    % Create a Kalman filter object.
      predictor = configureKalmanFilter('ConstantVelocity', centroid, [20, 20], [20, 20], 5);
    end

    
    
    
    function self=setupSystemObjects(self,vid)
    % Initialize Video I/O
    % Create objects for reading a video from a file, drawing the tracked
    % objects in each frame, and playing the video.

      
    %% Create a video file reader.
      self.videoHandler = [];
      [self.videoHandler,self.timerange] = self.newVideoHandler(vid,self.timerange);
      self.videoFile = vid;

      if self.stmif
        if isempty(self.stimulusPresenter) 
          self.stimulusPresenter = FishStimulusPresenter();
        end
        self.stimulusPresenter.init();
      end
      
      
      %% look for objects 
      if isempty(self.nfish) || isempty(self.fishlength) || isempty(self.fishwidth) || self.useScaledFormat  
        [nfish,fishSize] = self.findObjectSizes();

      
        if (nfish>25 || ~nfish) && isempty(self.nfish)% too many... something wrong
          warning('The fish size and number cannot be determined')
          if self.displayif
            self.nfish = chooseNFish(vid,1); % only if interactively
            fishSize = [100,20]; % wild guess;
            close(gcf);
          else
            error('Please manual provide fishlength, fishwidth and nfish');
          end
        end
      end  
      
      if isempty(self.opts.fishlength) % otherwise already set by hand
        self.fishlength = fishSize(1);
      end
      if isempty(self.opts.fishwidth) 
        self.fishwidth = fishSize(2); 
      end
      if isempty(self.opts.nfish) 
        self.nfish = nfish;
      end
      assert(self.fishlength>self.fishwidth);

      % overwrite the given optiond
      self.fishlength = self.fishlength;
    
      % set all options
      self.setOpts();
    end

    
    function initTracking(self)
    % MAYBE RESET THE CLASSIFIER ? 

      self.videoHandler.timeRange = self.timerange;                       
      self.timerange = self.videoHandler.timeRange; % to get the boundaries right;
      self.videoHandler.setCurrentTime(self.timerange(1));
      self.videoHandler.fishlength = self.fishlength;
      self.videoHandler.fishwidth = self.fishwidth;
      self.setOpts();
      self.videoHandler.initialize();
      
      
      if isscalar(self.writefile) && self.writefile
        [a,b,c] = fileparts(self.videoFile);
        self.writefile = [a '/' b '_trackingVideo' c];
      end

      if ~isempty(self.writefile) && self.displayif
        self.videoWriter = self.newVideoWriter(self.writefile);
      else
        if ~isempty(self.writefile)
          warning('set displayif>0 for writing a video file!');
        end
        self.videoWriter = [];
      end

     
      
      %% get new fish classifier 
      self.fishClassifier = newFishClassifier(self,self.opts.classifier);
      self.isInitClassifier = isInit(self.fishClassifier);
      
      self.setDisplayType();
      if self.displayif && ~isOpen(self.videoPlayer)
        self.videoPlayer.show();
      end

      self.fishId2TrackId = nan(250,self.nfish);
      self.pos = [];
      self.res = [];
      self.tracks = [];
      self.savedTracks = [];
      
      self.initializeTracks(); % Create an empty array of tracks.
      self.nextId = 1;
      self.uniqueFishFrames = 0;
      self.meanCost(:) = 1;
      self.meanAssignmentCost = 1;

      self.segments = [];
      self.bboxes = [];
      self.centroids = [];
      self.idfeatures = [];
      self.centerLine = [];
      self.thickness = [];
      self.classProb = [];
      self.classProbNoise = [];
      self.crossings = [];
      self.saveFieldSegIf = [];
      self.saveFieldSub = {};
      
      if length(self.scalecost)< length(self.costinfo)
        verbose('Enlarging scalecost with 1s.')
        scale = self.scalecost;
        self.scalecost = ones(1,length(self.costinfo));
        self.scalecost(1:length(scale)) =  scale;
      end
      
      self.meanCost = ones(1,length(self.scalecost));
      
      self.assignments = [];
      self.unassignedTracks = [];
      self.unassignedDetections = [];
      self.cost = [];
      self.currentFrame = 0;
      
      self.currentCrossBoxLengthScale = self.opts.tracks.crossBoxLengthScalePreInit;
      
      if isempty(self.maxVelocity)
        self.maxVelocity = self.fishlength*80; % Ucrit for sustained swimming is 20 BodyLength/s (see Plaut 2000). Thus set
                                               % twice for short accelaration bursts
      end
    
    end
    
    
    function [nObjects,objectSize] =findObjectSizes(self,minAxisWidth)

      if nargin==1
        minAxisWidth = 4;
      end
      
      verbose('Detect approx fish sizes (minAxisWidth=%g)...',minAxisWidth);

      self.videoHandler.reset();
      self.videoHandler.originalif = true;
      n = min(self.videoHandler.history,floor(self.videoHandler.timeRange(2)*self.videoHandler.frameRate));
      self.videoHandler.initialize(0);
      s = 0;
      for i = 1:n

        [segm] = self.videoHandler.step();
        fprintf('%1.1f%%\r',i/n*100); % some output

        if isempty(segm)
          continue;
        end
        s = s+1;
        idx = [segm.MinorAxisLength]>minAxisWidth;
        count(s) = length(segm(idx));
        width(s) = median([segm(idx).MinorAxisLength]);
        height(s) = median([segm(idx).MajorAxisLength]);
      end

      if ~s
        error('Something wrong with the bolbAnalyzer ... cannot find objects.');
      end
      
      nObjects = median(count);
      objectSize = ceil([quantile(height,0.9),quantile(width,0.9)]); % fish are bending
      verbose('Detected %d fish of size %1.0fx%1.0f pixels.',nObjects, objectSize);
      
      if self.displayif>2
        figure;
        bwmsk = self.videoHandler.getCurrentBWImg();
        imagesc(bwmsk);
        title(sprintf('Detected %d fish of size %1.0fx%1.0f pixels.',nObjects, objectSize));
      end

      if self.useScaledFormat
        % get good color conversion
        bwmsks = {};
        cframes = {};
        for i =1:10 % based on not just one frame to get better estimates
          [segm] = self.videoHandler.step();
          bwmsks{i} =self.videoHandler.getCurrentBWImg();
          cframes{i} = im2double(self.videoHandler.getCurrentFrame());
        end
        
        [scale,delta] = getColorConversion(bwmsks,cframes);
          
        if ~isempty(scale)
          self.videoHandler.setToScaledFormat(scale,delta);
          self.videoHandler.resetBkg(); 
          self.videoHandler.computeSegments = false;
          % re-generate some background
          verbose('Set to scaled format. Regenerate background..')
          for i =1:(n/2)
            [~,frame] = self.videoHandler.step();    
            fprintf('%1.1f%%\r',i/n*2*100); % some output
          end
          self.videoHandler.computeSegments = true;
        else
          self.videoHandler.setToRGBFormat();
          self.videoHandler.frameFormat = ['GRAY' self.videoHandler.frameFormat(end)];
        end
      end

      self.videoHandler.originalif = false;
      self.videoHandler.reset();
    end
    
    
    function overlap = bBoxOverlap(self,bbtracks,bbsegs,margin)
      
      overlap = zeros(size(bbtracks,1),size(bbsegs,1));
      for i = 1:size(bbtracks,1)
        
        % better judge overlap from bounding box 
        bbtrack = double(bbtracks(i,:));
        bbtrack(1:2) = bbtrack(1:2) - margin/2;
        bbtrack(3:4) = bbtrack(3:4) + margin;
        
        mnx = min(bbsegs(:,3),bbtrack(3));
        mny = min(bbsegs(:,4),bbtrack(4));
      
        bBoxXOverlap = min(min(max(bbtrack(1) + bbtrack(3) - bbsegs(:,1),0), ...
                               max(bbsegs(:,1) + bbsegs(:,3) - bbtrack(1),0)),mnx);
        bBoxYOverlap = min(min(max(bbtrack(2) + bbtrack(4) - bbsegs(:,2),0), ...
                               max(bbsegs(:,2) + bbsegs(:,4) - bbtrack(2),0)),mny);
        
        overlap(i,:) = bBoxXOverlap.*bBoxYOverlap./mnx./mny;
      end
    end
    
    
    function detectObjects(self)
    % calls the detectors

      segm = self.segments;
      
      
      if ~isempty(segm)

        % check for minimal distance. 
        self.centroids = cat(1,segm.Centroid);
        self.bboxes = int32(cat(1,segm.BoundingBox));
        
        
        overlap = self.bBoxOverlap(self.bboxes,self.bboxes,2);
        overlap(1:length(segm)+1:end) = 0;
        [i,j] =  find(overlap);
        
        while ~isempty(i) && length(self.segments)>self.nfish

          if segm(i(1)).Area > segm(j(1)).Area
            deli = j(1);
          else
            deli = i(1);
          end
          
          if self.displayif >4
            figure(23)
            imagesc(self.videoHandler.getCurrentFrame());
            hold on;
            plot(self.centroids(:,1),self.centroids(:,2),'x');
            rectangle('position',self.bboxes(i(1),:),'edgecolor','k')
            rectangle('position',self.bboxes(j(1),:),'edgecolor','k')
            rectangle('position',self.bboxes(deli,:),'edgecolor','r')
            drawnow;
          end
          
          segm(deli) = [];
          self.segments(deli) = [];
          self.centroids(deli,:) = [];
          self.bboxes(deli,:)= [];
          
          % recompute
          overlap = self.bBoxOverlap(self.bboxes,self.bboxes,self.fishwidth);
          overlap(1:length(segm)+1:end) = 0;
          [i,j] =  find(overlap);
        end
                                
        
        if self.useMex % there is no centerline information without
                       % the mex version...
          self.centerLine  = cat(3,segm.CenterLine);
          self.thickness  = cat(1,segm.Thickness);

          if size(self.centerLine,3)<length(segm)
            % center Line not computed for some of the detections
            self.centerLine = nan(self.videoHandler.nprobe,2,length(segm));
            self.thickness = nan(length(segm),self.videoHandler.nprobe);
            for i = 1:length(segm)
              if ~isempty(segm(i).CenterLine)
                self.centerLine(:,:,i) = segm(i).CenterLine;
                self.thickness(i,:) = segm(i).Thickness;
              end
            end
          end
          
        else
          self.centerLine = [];
          self.thickness = [];
        end
        


        if ~self.useMex
          for i = 1:length(segm)
            % better take MSER regions if existent
            if ~isempty(segm(i).MSERregions)
              self.centroids(i,:) = double(segm(i).MSERregionsOffset + segm(i).MSERregions(1).Location);
            end
          end
        end
        
        

        % MAYBE TAKE DCT2 !???
        idfeat = permute(cat(4,segm.FishFeature),[4,1,2,3]);
        self.idfeatures =  idfeat;
        
        for ct = self.featurecosttypes
          self.features.(ct{1}) = cat(1,segm.(ct{1}));
        end
        
        
        self.classProb = predict(self.fishClassifier,self.idfeatures(:,:));
        %self.classProbNoise = cat(1,segm.MinorAxisLength)./cat(1,segm.MajorAxisLength);
        self.classProbNoise = cat(1,segm.bendingStdValue);
      else
        self.centroids = [];
        self.bboxes = [];
        self.idfeatures = [];
        self.features = [];
        self.classProb = [];
        self.centerLine = [];
        self.thickness = [];
        self.classProb = [];
        self.classProbNoise = [];
      end
      
      

    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ClASSIFICATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function cla = newFishClassifier(self,opts);
    % Create Classifier objects
      featdim = [self.videoHandler.getFeatureSize()];
      cla = FishBatchClassifier(self.nfish,featdim);
      
      % set proporties
      for f = fieldnames(opts)'
        if isprop(cla,f{1})
          cla.(f{1}) = opts.(f{1});
        end
      end

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

    
    function nc = connectedComponents(self,trackIndices,assignedFishIds)
      oldFishIds = [self.tracks(trackIndices).fishId];
      A = zeros(self.nfish);
      A(sub2ind(size(A),oldFishIds,assignedFishIds)) = 1;
      [~,~,nc] = networkComponents(A | A');
      nc(cellfun('length',nc)==1) = []; % delete self-components
    
    end
      
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
      
      % check whether it is only swapping
      nc = self.connectedComponents(trackIndices,assignedFishIds);

      % build Gaussian kernel
      sigma = self.tracks(1).classProbHistory.tau; % same for dist/clprob based
      dt = 1; % in frames
      l = 0:dt:sigma*8;
      l = [l(end:-1:2) l];
      G = (exp(-(l.^2/2/sigma^2))/sqrt(2*pi)/sigma*dt)';
      g2 = ceil(length(G)/2);

      %% switch the tracks for each connected component
      f2t = self.fishId2TrackId;
      pos = self.pos;

      for i = 1:length(nc)
        msk = ismember(oldFishIds,nc{i});

        if crossingflag
          tmin = calcDistanceBasedSwitchPoint(trackIndices(msk),assignedFishIds(msk));
        else
          tmin = calcClassProbBasedSwitchPoint(trackIndices(msk),assignedFishIds(msk));
        end
        
        idx = find(msk);
        t = tmin:self.currentFrame;        
        % swap
        for j = idx(:)'
          self.fishId2TrackId(t,assignedFishIds(j)) = f2t(t,oldFishIds(j));
          self.pos(:,assignedFishIds(j),t) = pos(:,oldFishIds(j),t);
          self.tracks(trackIndices(j)).fishId = assignedFishIds(j);
          self.resetBatchIdx(trackIndices(j)); 
        end
      end
      
      
      function tmin = calcDistanceBasedSwitchPoint(inTrackIndices,inAssignedFishIds)
      % to be called from switchFish: calculated the end point based on the distance of
      % the tracks. It is assumed that the given indices are only ONE permutation
        localOldFishIds = [self.tracks(inTrackIndices).fishId];

        tcurrent = self.currentFrame;
        % should be the same crossing, since we expect just one fully connected component
        if isempty([self.tracks(inTrackIndices).crossedTrackIds])
          % not currently crossing. Thus: compute from current frame
          tend = tcurrent;
          tstart = max([self.tracks(inTrackIndices).lastFrameOfCrossing]);
        else
          % during the crossing
          tstart = min([self.tracks(inTrackIndices).firstFrameOfCrossing]);
          tend =  max([self.tracks(inTrackIndices).lastFrameOfCrossing]);
        end
        
        % check min distance within crossing
        dist = squeeze(sum((self.pos(:,localOldFishIds,tstart:tend) - self.pos(:,inAssignedFishIds,tstart:tend)).^2,1))';
        dist = [repmat(dist(1,:),g2,1);dist;repmat(dist(end,:),g2,1)];
        cdist = conv2(dist,G);
        cdist = cdist(2*g2:end-2*g2+1,:);
        mdist = mean(cdist,2); %just one component and symmetric distance. 
        [~,tminsub] = min(mdist);
        tmin = tstart + tminsub-1;

        verbose('N(dist) = %d',tmin-tcurrent);
      
      end % nested function
        
      function tmin = calcClassProbBasedSwitchPoint(inTrackIndices,inAssignedFishIds)
      % get time from classification accuracy 
        
        localOldFishIds = [self.tracks(inTrackIndices).fishId];
        
        % at most from last crossing time
        tstart = min([self.tracks(inTrackIndices).lastFrameOfCrossing]); 
        tcurrent = self.currentFrame;
        
        tspan  = tcurrent-tstart; % from currentframe!!
        pdiff = [];
        for i = 1:length(inTrackIndices)
          trIdx = inTrackIndices(i);
          [p,w] = self.tracks(trIdx).classProbHistory.getData(tspan);
          if isempty(pdiff)
            pdiff = zeros(length(w),length(inTrackIndices));
          end
          pdiff(:,i) = p(:,inAssignedFishIds(i)) -  p(:,localOldFishIds(i));
          pdiff(w<self.tracks(trIdx).classProbHistory.reasonableThres,i) = 0;
        end
        
        % need to ensure that the last diff values are not noisy and below
        % zeros. Otherwise the convolution boundary will fail..
        mpdiff = mean(pdiff,2);
        if mpdiff(end)<0
          idx = find(mpdiff>0,1,'last');
          mpdiff(idx:end) = mpdiff(idx);
        end
        
        mpdiff = [mpdiff(1)*ones(g2,1);mpdiff;mpdiff(end)*ones(g2,1)];
        cpdiff = conv2(mpdiff,G);
        cpdiff = cpdiff(2*g2:end-2*g2+1);
        
        tpmin = find(cpdiff(end:-1:1)<=nanmean(cpdiff),1,'first');
        tmin = tcurrent-tpmin+1; 

        verbose('N(clprob) = %d',tmin-tcurrent);

        
      end % nested function

      
    end
    
    
    function resetBatchIdx(self,trackidx)
      [self.tracks(trackidx).nextBatchIdx] = deal(1);
    end


    
    function changed = fishClassifierUpdate(self,trackIndices,updateif,forceif)
    % updates and tests the collected features according to the current fishIDs. It does NOT
    % force and update, update is only done if the batch is likly to come from the true
    % fish. In case a valid permutation is detected and reassignProb is higher than a
    % threshold, fishIDs are exchanged.
    %
    % resets the batchIdxs after update

      fishIds = cat(2,self.tracks(trackIndices).fishId);
      changed = 0;
      % test on the whole set of possible fishIds
      [assignedFishIds prob minsteps probdiag] = self.predictFish(trackIndices,...
                                                        1:self.nfish,self.opts.classifier.nFramesForSingleUpdate);
      same = assignedFishIds==fishIds; 

      if updateif && all(same) || self.enoughEvidenceForForcedUpdate(prob,minsteps,probdiag) 
        % do not update if classes might be mixed up (wait for more data)
        % change classifier
        feat = self.getFeatureDataFromTracks(trackIndices);
        [assignedFishIds prob] = self.fishClassifier.batchUpdate(fishIds,feat,1); % force
        self.resetBatchIdx(trackIndices);
        changed  = 1;
        return
      end
      
      if length(same)>1 && ~all(same) &&  (forceif || self.enoughEvidenceForReassignment(prob,minsteps,probdiag))

        verbose('mixed and min prob is %1.2f',min(prob));
        if all(ismember(assignedFishIds,fishIds))
          % only permute if true permutation  
          verbose('Switching fish with fishClassifierUpdate')
          self.switchFish(trackIndices,assignedFishIds,self.opts.tracks.purelyDistanceBasedSwitching);
          self.resetBatchIdx(trackIndices);
          changed = 1;
        else
          verbose('like to switch fish, but could not find a permutation.');
        end
      elseif length(same)>1 && ~all(same) 
        verbose('considered not mixed: min prob is %1.2f\r',min(prob));
      end
      
    end

    
    function updatePos(self)
    % save current track->fishID
      t = self.currentFrame;
      fishIds = cat(1,self.tracks.fishId);
      tracks = self.tracks(~isnan(fishIds));
      
      if size(self.fishId2TrackId,1)< t
        self.fishId2TrackId = cat(1,self.fishId2TrackId,nan(250,self.nfish));
      end
      % save track positions
      if isempty(self.pos)
        self.pos = nan(2,self.nfish,250);
      end
      
      if size(self.pos,3)<t
        self.pos = cat(3,self.pos,nan(2,self.nfish,250));
      end
      
      if isempty(tracks) % no tracks
        return;
      end
      
      fishIds = cat(2,tracks.fishId);      
      trackIds = cat(2,tracks.id);      
      self.fishId2TrackId(t,fishIds) = trackIds;
      
      self.pos(1:2,fishIds,t) = cat(1,tracks.centroid)';
      
     
    end
    
    
    function success = initClassifier(self)
    % checks condition and inits classifier (if conds are met) and resets relevant counters

      success = false;
      thres = 5;
      if length(self.tracks)<self.nfish 
        if  self.currentFrame > thres*self.opts.classifier.nFramesForInit
          verbose(['nfish setting might be wrong or some fish are lost\r'])
        end
        return
      end

      if self.uniqueFishFrames>self.opts.classifier.nFramesForInit  || ...
          self.currentFrame> thres*self.opts.classifier.nFramesForInit % cannot wait for ages..
        
        fishIds = cat(2,self.tracks.fishId);      

        feat = self.getFeatureDataFromTracks(1:length(self.tracks));
        [~,sidx] = sort(fishIds);
        feat = feat(sidx);
        batchsample = cell(1,self.nfish);
        [batchsample{:}] = deal([]);
        batchsample(1:length(feat)) = feat;
        
        self.fishClassifier.init(batchsample);
        self.currentCrossBoxLengthScale = self.opts.tracks.crossBoxLengthScale; % update
        self.resetBatchIdx(1:length(self.tracks));
        
        % reset all potentials previous crossings
        [self.tracks.lastFrameOfCrossing] = deal(self.currentFrame);
        [self.tracks.firstFrameOfCrossing] = deal(self.currentFrame);
        [self.tracks.crossedTrackIds] = deal([]);
        
        self.uniqueFishFrames  = 0;
        success = true;
      
      elseif self.currentFrame < thres/2*self.opts.classifier.nFramesForInit % cannot wait for ages..
        if ~isempty(self.crossings)
          self.resetBatchIdx(1:length(self.tracks));% we do not want a mixture at the beginning
        end
      end
      
      self.isInitClassifier = success;
      
    end
    
    function  [assignedFishIds, prob, minsteps, probdiag] = predictFish(self,trackIndices,fishIdSet,maxSteps);
    % predicts the fish according to the values saved in the classProbHistory up to the lastFrameOfCrossing
      
      %assumedFishIds = [self.tracks(trackIndices).fishId];
      % feat = self.getFeatureDataFromTracks(trackIndices);
      % [assignedFishIds1,prob1] = self.fishClassifier.predictPermutedAssignments(feat,assumedFishIds);
    
      % we could also use the class probhistory
      C = zeros(length(trackIndices),length(fishIdSet));
      minsteps = maxSteps;   

      assignedFishIds = nan(size(trackIndices));
      prob = nan(size(assignedFishIds));
      probdiag = nan(size(assignedFishIds));

      
      for i = 1:length(trackIndices)
        track = self.tracks(trackIndices(i));
        nsteps = max(self.currentFrame-track.lastFrameOfCrossing,1);
        [p nsteps] = track.classProbHistory.mean(nsteps);
        minsteps = min(nsteps,minsteps);
        
        if minsteps==0
          return;
        end
        
        C(i,:) = 1-p(fishIdSet);
      end
      assignments = assignDetectionsToTracks(C',2);
      
      assignedFishIds(assignments(:,2)) = fishIdSet(assignments(:,1)); % fishIDs have 1:nfish order
      for i = 1:size(assignments,1)
        prob(assignments(i,2)) = 1-C(assignments(i,2),assignments(i,1));
        probdiag(i) = 1-C(i,i);
      end
    end

    
    function bool = enoughEvidenceForReassignment(self,prob,minsteps,probdiag)
      if any(isnan(prob))
        bool = false;
      else
        bool = min(prob)>self.opts.classifier.reassignProbThres && minsteps>=self.opts.classifier.minBatchN;
      end
    end
    
    function bool = enoughEvidenceForForcedUpdate(self,prob,minsteps,probdiag)
      if any(isnan(prob)) || length(prob)~=self.nfish
        bool = false;
      else
        bool = minsteps >= self.opts.classifier.nFramesAfterCrossing;
        bool = bool && max(prob)<=self.opts.classifier.reassignProbThres;
      end
      
      
    end
    
    function bool = enoughEvidenceForBeingHandled(self,prob,minsteps,probdiag)
      if any(isnan(prob))
        bool = false;
      else
        bool = (min(prob)>self.opts.classifier.reassignProbThres ...
                && minsteps>=self.opts.classifier.nFramesAfterCrossing) ... 
               || minsteps>=self.opts.classifier.nFramesForUniqueUpdate; % to avoid very long crossings
      end
    end
    
    
    function handlePreviouslyCrossedTracks(self)
    % expect that the trackIndices are not currently crossing (they
    % can have just reached a newCrossing though).
      
      trackIndices = find(self.testBeyondBound() & ~self.testHandled());

      if isempty(trackIndices)
        return;
      end
      
      if self.displayif>6
        self.plotCrossings_(1,trackIndices);
      end

      toTrackIndices = trackIndices;
      tranges = zeros(length(trackIndices),2);
      trackIds = cat(2,self.tracks.id);      
      fishIds = cat(2,self.tracks.fishId);      
      
      
      crossedTrackIdStrs = arrayfun(@(x)num2str(sort(x.crossedTrackIds)), self.tracks(trackIndices),'uni',0);
      [u,idxct,idxu] = unique(crossedTrackIdStrs);
      
      
      % we should distinguish the following cases
      % 1)  only one trackindex. (numel(trackIndices)==1)
      %    in this case, we should check :
      %      a) whether the track is very likely a certain
      %      fishID. If yes, we can erase it from the other tracks
      %      crossed track fields,
      %      otherwise, we should wait for the other tracks to
      %      finish. (for instance by just increasing the
      %      lastCrossed by one). [If this track does come into a
      %      new crossing, then all crossed tracks will also enter this new
      %      crossing and the tracks will be handled later]
      % 2)  if there are more than one trackindex. We test them in
      %     a permuted way, allowing for all fish currently in the
      %     crossing. If it is not clear which fish left the crossing,
      %     we just leave it and wait for the other to finish. If
      %      instead we find a true permutation of the fish that are
      %     leaving, we delete them out of the crossed fields of all
      %     others. 
      % 3)  if there is only one fish in the crossed id field, then,
      %     if it is the same trackid, it's OK. just delete it. If it
      %     is another trackid, we have a problem. One should force a
      %     all fish update (maybe we could just put ALL fish into the
      %     crossID field and let it handle the next time.)
      
      for i = 1:length(u)
        crossedTrackIds =  self.tracks(trackIndices(idxct(i))).crossedTrackIds;
        thisTrackIndices = trackIndices(idxu==i);

        
        tstart = min([self.tracks(thisTrackIndices).firstFrameOfCrossing]);
        tend = max([self.tracks(thisTrackIndices).lastFrameOfCrossing]);
        
        [~,crossedFishIds] = find(ismember(self.fishId2TrackId(tstart:tend,:),crossedTrackIds));
        crossedFishIds = unique(crossedFishIds)';

        assumedFishIds = fishIds(thisTrackIndices);
        if ~all(ismember(assumedFishIds,crossedFishIds))
          %shoudl no happen because fishIDs swapping should be blocked
          %during crossing
          keyboard;
        end


        % we force the permutation to be in the valiud fish only (otherwise too many errors for many fish)
        [assignedFishIds prob minsteps probdiag] = self.predictFish(thisTrackIndices,crossedFishIds,...
                                                          self.opts.classifier.nFramesForUniqueUpdate); 
        
        if  all(ismember(assignedFishIds,assumedFishIds)) 
          % we have a valid permutation and enough confidence (or the very same ordering)
          nc = self.connectedComponents(thisTrackIndices,assignedFishIds);
          
          for i_nc = 1:length(nc)
            msk = ismember(assumedFishIds,nc{i_nc});
            
            if self.enoughEvidenceForReassignment(prob(msk),minsteps,probdiag(msk))
              % permutation : switch tracks and delete from others. 
              verbose('valid permutation.. switch...')
              self.switchFish(thisTrackIndices(msk),assignedFishIds(msk),true); % set the crossing flag
              deleteCrossedTrackIds(thisTrackIndices(msk));
            end
          end
          equalIndices = find(assignedFishIds == assumedFishIds);
          for i = equalIndices(:)'
            if self.enoughEvidenceForReassignment(prob(i),minsteps,probdiag(i)) ...
                || length(self.tracks(thisTrackIndices(i)).crossedTrackIds)==1
              deleteCrossedTrackIds(thisTrackIndices(i));
            end
          end
% $$$           
% $$$           if isempty(nc) && self.enoughEvidenceForBeingHandled(prob,minsteps,probdiag)
% $$$             % signal that things are handled if enough confidence
% $$$             deleteCrossedTrackIds(thisTrackIndices);
% $$$           end
% $$$           
        else
          % now the assignedFishIds are outside of the assumedFishIds but still insight the possible crossed. That
          % means that there is at least another fish involved in the crossing. (see above the asssumedFishIds are
          % necessaryly in the crossed, so we do not need to account for that problem)
          
          % wait. Will trigger next time
          for i = thisTrackIndices(:)'
            self.tracks(i).lastFrameOfCrossing = self.tracks(i).lastFrameOfCrossing+ 1;
          end
          verbose('wait for others to exit crossing\r')
        end

      end
    
    
    
      if self.displayif>6
        self.plotCrossings_(2,trackIndices);
      end

      
      function deleteCrossedTrackIds(localTrackIndices)
        [self.tracks(localTrackIndices).crossedTrackIds] = deal([]);
        thisTrackIds = [self.tracks(localTrackIndices).id];

        %delete these ids from all others. They should not currently cross. 
        for i = 1:length(self.tracks)
          self.tracks(i).crossedTrackIds = setdiff(self.tracks(i).crossedTrackIds,thisTrackIds);
          if isempty(self.tracks(i).crossedTrackIds)
            self.tracks(i).crossedTrackIds = []; % some weired things happening with the size of the empty set in matlab

          end
        end

      end % nested function
      
    end
    
    
    
    function plotCrossings_(self,section,updateIdx)


      figure(2);
      if section==1
        clf;
      end
      
      cols = jet(self.nfish);
      [r1,r2] = getsubplotnumber(length(updateIdx));
      
      s = 0;
      for trackIndex = updateIdx
        s = s+1;
        a = subplot(r1,r2,s);
        set(a,'ColorOrder',cols,'ColorOrderIndex',1);
        hold on;
        
        t1 = self.tracks(trackIndex).firstFrameOfCrossing;
        t2 = self.tracks(trackIndex).lastFrameOfCrossing;
        t0 = max(t1-30,1);
        t = self.currentFrame;
        
        if section==1
          ids = self.tracks(trackIndex).crossedTrackIds;
          
          %first plot all traces
          plot(squeeze(self.pos(1,:,t0:t))',squeeze(self.pos(2,:,t0:t))','--','Linewidth',1);
          hold on; 
          
          % plot those that are in the current crossing
          fids = find(ismember(self.fishId2TrackId(t2,:),ids));
          set(a,'ColorOrder',cols(fids,:));set(a,'ColorOrderIndex',1)
          plot(squeeze(self.pos(1,fids,t1:t2))',squeeze(self.pos(2,fids,t1:t2))','o','Linewidth',2);
          
        else
          %section 2
          fid = self.tracks(trackIndex).fishId;
          plot(squeeze(self.pos(1,fid,t0:t))',squeeze(self.pos(2,fid,t0:t))','-','Linewidth',2,'color',cols(fid,:));
          frameSize = self.videoHandler.frameSize;
          axis ij;
          xlim([0,frameSize(2)])
          ylim([0,frameSize(1)]);
        

        end
        
      end
      if section==2
        drawnow;
      end
    end
    
    function handled = testHandled(self);
      handled = false(1,length(self.tracks));
      for i = 1:length(self.tracks)
        handled(i) = isempty(self.tracks(i).crossedTrackIds);
      end
      
      %handled = arrayfun(@(x)isempty(x.crossedTrackIds),self.tracks);
    end
    
    function  hitBound = testBeyondBound(self);
      % hit bound will be true even if already handled!!
      hitBound = self.currentFrame-[self.tracks.lastFrameOfCrossing] > self.opts.classifier.nFramesAfterCrossing;
    end
      
    function newCrossing = testNewCrossing(self)
      newCrossing  = ismember(1:length(self.tracks), cat(2,self.crossings{:}));
    end
    
    function [bool,pdiff] = testTrackMisalignment(self,trackIndices);
      fishIds = cat(1,self.tracks(trackIndices).fishId);
      clp = cat(1,self.tracks(trackIndices).clpMovAvg);
      idx = s2i(size(clp),[(1:size(clp,1))',fishIds]);
      prob_correct = mean(clp(idx));
      clp(idx) = 0;
      prob_notcorrect = mean(max(clp,[],2)); % might be Nan
      pdiff = prob_notcorrect - prob_correct ;
      bool = pdiff > self.opts.classifier.reassignProbThres/2; % will be zero if pdiff should be NaN
    end
    
      
    function updateClassifier(self)
    % updates the identity classifier and tracks the identity to
    % trackID assisgments
      
      self.crossings = self.detectTrackCrossings(); % only detects crossings

      % update unique frames
      crossif = ~isempty(self.crossings);
      if ~crossif && length(self.tracks)==self.nfish
        self.uniqueFishFrames = self.uniqueFishFrames + 1;
      else
        self.uniqueFishFrames = 0;
      end
      
      if ~self.isInitClassifier
        if ~self.initClassifier()
          return
        end
      end


      %% handle previous crossings
      self.handlePreviouslyCrossedTracks();

      handled = self.testHandled();
      allhandled = all(handled);


      %% update the current crossings
      self.updateTrackCrossings(); % updates the self.tracks fields

      
      %% test whether tracks might be misaligned
      if sum(handled)>1
        handledIndices = find(handled);
        [misif,pdiff] = self.testTrackMisalignment(handledIndices);
        if misif 
          if pdiff < self.opts.classifier.reassignProbThres
            verbose('Probably misaligned (%1.2f): better test again and switch..\r',pdiff);
            forceif = 0;
          else
            verbose('Surely misaligned (%1.2f): try to force switch..',pdiff);
            forceif = 1;
          end
          if self.fishClassifierUpdate(handledIndices,0,forceif); % use function for testing/switching 
            self.uniqueFishFrames = 0;
          end
        end
      end
      
      
        
      %% unique fish update
      if (self.uniqueFishFrames-1 > self.opts.classifier.nFramesForUniqueUpdate)  && allhandled
        if self.fishClassifierUpdate(1:length(self.tracks),1,0);
          verbose('Performed unique frames update.\r')
        end
        self.uniqueFishFrames = 0; % anyway reset;
      end

      
      %% single fish update
      enoughData =  [self.tracks.nextBatchIdx] > self.opts.classifier.nFramesForSingleUpdate;
      singleUpdateIdx = find(handled & enoughData);
      if ~isempty(singleUpdateIdx)
        if ~self.fishClassifierUpdate(singleUpdateIdx,1,0) 
          if  allhandled
            self.fishClassifierUpdate(1:length(self.tracks),0,0); % something might be wrong. Try all
          end
          self.resetBatchIdx(singleUpdateIdx); % anyway reset          
        else
          verbose('Performed single fish update.\r');          
        end
      end

    end

    
    function cleanupCrossings(self)
      %switch if necessary
      trackIndices = 1:length(self.tracks);
      [assignedFishIds,prob,~,probdiag] = self.predictFish(trackIndices,1:self.nfish,Inf);
      if self.enoughEvidenceForReassignment(prob,Inf,probdiag)
        self.switchFish(trackIndices,assignedFishIds,self.opts.tracks.purelyDistanceBasedSwitching);
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
        'velocity',{},...
        'centerLine',{},...
        'thickness',{},...
        'segment', {},...
        'features',{},...
        'classProb',{},...
        'classProbHistory',{},...
        'clpMovAvg',{},...
        'batchFeatures',{},...
        'nextBatchIdx',{},...
        'predictor', {}, ...
        'firstFrameOfCrossing', {}, ...
        'lastFrameOfCrossing', {}, ...
        'crossedTrackIds',{},...
        'age', {}, ...
        'totalVisibleCount', {}, ...
        'consequtiveInvisibleCount', {},...
        'assignmentCost', {},...
        'nIdFeaturesLeftOut', {},...
        'stmInfo',{});
    end
    
    
    function predictNewLocationsOfTracks(self)

      szFrame = self.videoHandler.frameSize;
      for i = 1:length(self.tracks)
        bbox = self.tracks(i).bbox;
        
        % Predict the current location of the track.
        predictedCentroid = predict(self.tracks(i).predictor);
        
        % Shift the bounding box so that its center is at
        % the predicted location.
        predictedCentroid(1:2) = int32(predictedCentroid(1:2)) - bbox(3:4) / 2;
        % no prediction outside of the video:

        predictedCentroid(1) = max(min(predictedCentroid(1),szFrame(2)),1);
        predictedCentroid(2) = max(min(predictedCentroid(2),szFrame(1)),1);

        self.tracks(i).bbox = [predictedCentroid(1:2), bbox(3:4)];
        self.tracks(i).centroid = predictedCentroid;
      end
    end
      
    
    function fcost = computeCostMatrix(self)
      
      nTracks = length(self.tracks);
      nDetections = size(self.centroids, 1);
      fcost = zeros(nTracks, nDetections,length(self.costinfo));
      highCost = 5;
      somewhatCostly = 2;
      
      if nDetections>0
        %mdist = min(pdist(double(centroids)));
        maxDist = self.maxVelocity/self.videoHandler.frameRate; % dist per frame
% $$$         minDist = 2;
% $$$         meanDist = sqrt(sum(cat(1,self.tracks.velocity).^2,2)); % track.velocity is in per frame
% $$$         %nvis = max(cat(1,self.tracks.consequtiveInvisibleCount)-self.opts.tracks.invisibleForTooLong,0);
% $$$         nvis = cat(1,self.tracks.consequtiveInvisibleCount);
% $$$         thresDist = min(max(meanDist,minDist),maxDist).*(1+nvis);
        thresDist = maxDist;
        
        % Compute the cost of assigning each detection to each track.
        for k = 1:length(self.costinfo)
          switch lower(self.costinfo{k})

            case {'location','centroid'}

              %% center position 
              centers = cat(1,self.tracks.centroid);
              dst = pdist2Euclidean(centers,self.centroids);
              %dst = bsxfun(@rdivide,dst,thresDist).^2;
              dst = (dst/thresDist).^2;

              dst(dst>5) = 5;
              
              fcost(:,:,k) = dst;

             
            case {'centerline'}
              if self.useMex
                
                %dst = zeros(nTracks,nDetections);
                %dstrev = zeros(nTracks,nDetections);
                %centerLineRev = self.centerLine(end:-1:1,:,:);
                %for i = 1:nTracks
                %  dst(i,:) = mean(sqrt(sum(bsxfun(@minus,self.centerLine,self.tracks(i).centerLine).^2,2)),1);
                %  dstrev(i,:) = mean(sqrt(sum(bsxfun(@minus,centerLineRev,self.tracks(i).centerLine).^2,2)),1);
                %end
                %dst = min(dst,dstrev);

                %% use the mex-file instead of the above
                dst = pdist2CenterLine(cat(3,self.tracks.centerLine),self.centerLine);
                % CAUTION: DOES NAN ALWAY "WORK" IN C ??!??  seems to work on my system tough... 
                
                nidx = any(isnan(dst),2);
                %dst = bsxfun(@rdivide,dst,thresDist).^2; % better square?
                dst = (dst/thresDist).^2; % better square?

                dst(dst>5) = 5;
              
                % makes not sense to compare the distance to others if one is nan
                dst(nidx,:) = 1;
                
                fcost(:,:,k) = dst;
                
              else
                fcost(:,:,k) = 1;
              end

              
            case {'overlap'}

              bbsegs = cat(1,self.segments.BoundingBox);
              
              for i = 1:nTracks
                %% cost of overlap
                %overlap = zeros(1,length(segments));
                %for jj = 1:length(segments)
                % memb = ismember(tracks(i).segment.PixelIdxList(:)',segments(jj).PixelIdxList(:)');
                % overlap(jj) = sum(memb)/length(memb);
                %end
                
                % better judge overlap from bounding box 
                bbtrack = double(self.tracks(i).bbox);
                
                mnx = min(bbsegs(:,3),bbtrack(3));
                mny = min(bbsegs(:,4),bbtrack(4));
                
                bBoxXOverlap = min(min(max(bbtrack(1) + bbtrack(3) - bbsegs(:,1),0), ...
                                       max(bbsegs(:,1) + bbsegs(:,3) - bbtrack(1),0)),mnx);
                bBoxYOverlap = min(min(max(bbtrack(2) + bbtrack(4) - bbsegs(:,2),0), ...
                                       max(bbsegs(:,2) + bbsegs(:,4) - bbtrack(2),0)),mny);
                
                overlap = bBoxXOverlap.*bBoxYOverlap./mnx./mny;
                
                fcost(i,:,k) = 1-overlap ;%+ (overlap==0)*somewhatCostly;
              end
              
            case {'classprob','classification','classifier'}

              if self.isInitClassifier
                clProb = cat(1,self.tracks(:).classProb);
                dst = pdist2Euclidean(clProb,self.classProb);%correlation??!??
                dst(isnan(dst)) = 1; 
                
                % prob for fish
                msk = max(self.classProb,[],2)< self.opts.tracks.probThresForFish;
                if ~any(isnan(msk))
                  dst(:,msk) = somewhatCostly;
                end
                fcost(:,:,k) = dst; 
              else
                fcost(:,:,k) = 1;
              end
              
            otherwise
              % assume feature cost (might result in error if not..)
              detectedFeatures = self.features.(self.costinfo{k});
              feat = cat(1,self.tracks.features);
              trackFeatures = cat(1,feat.(self.costinfo{k}));
              if size(detectedFeatures,2)>1
                fcost(:,:,k) = pdist2Euclidean(trackFeatures,detectedFeatures);
              else
                fcost(:,:,k) = abs(bsxfun(@minus,trackFeatures,detectedFeatures'));
              end
          end % switch

        end % costtype k loop
      
      end
    end
    
    function plotWrongAssignments(self,idx)
      clf;
      imagesc(self.videoHandler.getCurrentFrame());
      hold on;

      detectionIdx = self.assignments(:, 2);
      trackIdx = self.assignments(:, 1);

      centroids = self.centroids(detectionIdx, :);
      latestcentroids = cat(1,self.tracks.centroid);
      plot(latestcentroids(:,1),latestcentroids(:,2),'<k');
      plot(centroids(:,1),centroids(:,2),'or');
      plot(self.centroids(:,1),self.centroids(:,2),'xm');

      plot(latestcentroids(idx,1),latestcentroids(idx,2),'.w');


      bboxes = self.bboxes(detectionIdx, :);
      latestbboxes = cat(1,self.tracks.bbox);

      for i = 1:size(latestbboxes,1)
        rectangle('position',latestbboxes(i,:),'edgecolor','k');
      end
      for i = 1:size(bboxes,1)
        rectangle('position',bboxes(i,:),'edgecolor','r');
      end
    end
    
      
    function detectionToTrackAssignment(self)

      nDetections = size(self.centroids, 1);

      self.cost = [];
      self.assignments= [];
      self.unassignedTracks = [];
      self.unassignedDetections = [1:nDetections];


      if nDetections && length(self.tracks)
        
        scalecost = self.scalecost/sum(self.scalecost);
        scale = scalecost./self.meanCost;

        % is there is currently a crossing then prob of non-assignement should be scaled down 
        costOfNonAssignment =  self.meanAssignmentCost*self.opts.tracks.costOfNonAssignment;

        if self.opts.tracks.adjustCostDuringCrossing
          if length(self.tracks)==self.nfish
            currrentlyCrossing = sum(self.currentFrame - [self.tracks.lastFrameOfCrossing] ...
                                     < self.opts.classifier.nFramesAfterCrossing);
            if any(currrentlyCrossing)
              costOfNonAssignment = costOfNonAssignment/min(currrentlyCrossing/2+1,3);
            end
          end
        end
        
        
        % determine cost of each assigment to tracks
        fcost = self.computeCostMatrix();

        % scale the cost
        sfcost = sum(bsxfun(@times,fcost,shiftdim(scale(:),-2)),3);

        % reduce the cost if invisible --> leads to oscillations...
        %sfcost1 = bsxfun(@rdivide,sfcost,(nvis+1)); % do not save this cost...
        
        % Solve the assignment problem. First for the visible
        nvis = cat(1,self.tracks.consequtiveInvisibleCount);
        validIdx = find(~nvis);

        [self.assignments, self.unassignedTracks, self.unassignedDetections] = ...
            assignDetectionsToTracks(sfcost(validIdx,:),  costOfNonAssignment);

        self.unassignedTracks = validIdx(self.unassignedTracks);
        self.assignments(:,1) = validIdx(self.assignments(:,1)); 

        invisibleIdx = find(nvis);       
        if ~isempty(self.unassignedDetections) && ~isempty(invisibleIdx)
          % Solve the assignment problem for the invisible tracks
          sfcost1 = sfcost(invisibleIdx,self.unassignedDetections);
          %sfcost1 = bsxfun(@rdivide,sfcost1,(nvis(invisibleIdx)+1)); 
          mnvis = max(nvis);
          costOfNonAssignment1 = costOfNonAssignment*(1+mnvis);


          [assignments, unassignedTracks, unassignedDetections] = ...
              assignDetectionsToTracks(sfcost1,  costOfNonAssignment1);

          self.assignments = [self.assignments;...
                              [invisibleIdx(assignments(:,1)),...
                              self.unassignedDetections(assignments(:,2))]];
          self.unassignedDetections = self.unassignedDetections(unassignedDetections);
          self.unassignedTracks = [self.unassignedTracks;invisibleIdx(unassignedTracks)];
        end
        

        % check max velocity of the assignments
        if ~isempty(self.assignments)
          detectionIdx = self.assignments(:, 2);
          centroids = self.centroids(detectionIdx, :);
          trackIdx = self.assignments(:, 1);
          latestCentroids = cat(1,self.tracks(trackIdx).centroid);
          velocity = sqrt(sum((centroids - latestCentroids).^2,2))./(nvis(trackIdx) +1);
          
          idx = find(velocity>self.maxVelocity/self.videoHandler.frameRate);

          if ~isempty(idx)
            self.unassignedTracks = [self.unassignedTracks; trackIdx(idx)];
            self.unassignedDetections = [self.unassignedDetections; detectionIdx(idx)];
            self.assignments(idx,:) = [];
            verbose('%d assignment(s) exceeded velocity.\r',length(idx));
          end
        end
        
        
        outcost = nan(1,size(sfcost,2));
        if ~isempty(self.assignments)
          for i = 1:size(self.assignments,1)
            outcost(self.assignments(i,2)) = sfcost(self.assignments(i,1),self.assignments(i,2));
          end
          for i = 1:length(self.unassignedDetections)
            outcost(self.unassignedDetections(i)) = min(sfcost(:,self.unassignedDetections(i)));
          end
        end

        self.cost = outcost;

        % save the individual mean cost
        mt = min(self.opts.tracks.costtau,self.currentFrame);
        self.meanCost = self.meanCost*(mt-1)/mt  + 1/mt*mean(reshape(fcost,[],size(fcost,3)),1);

        % save the mean assignment costs
        if ~isempty(self.assignments)
          costs = self.cost(self.assignments(:,2));
          if ~isempty(costs)
            self.meanAssignmentCost = self.meanAssignmentCost*(mt-1)/mt ...
                + 1/mt*mean(costs);
          end
        end
        
        if self.displayif>4 && ~isempty(self.assignments) 
          figure(19);
          trackIdx = self.assignments(:, 1);
          detectionIdx = self.assignments(:, 2);
          self.plotWrongAssignments(1:length(trackIdx));
          drawnow;
          
        end
        
        
        if self.displayif>4 && self.currentFrame>4
          figure(5);
          clf;
          [r1,r2] = getsubplotnumber(length(self.segments));
          for i = 1:length(self.segments)
            subplot(r1,r2,i);
            imagesc(self.segments(i).FilledImage);
            if any(i==self.unassignedDetections)
              title(sprintf('unassigned [%1.3f]',outcost(i)));
            else
              title(sprintf('trackID %d [%1.3f]',self.assignments(self.assignments(:,2)==i,1),outcost(i)));
            end
            
          end
          drawnow;
          
          if self.displayif>5
            figure(6);
            subplot(self.nfish+2,1,1);
            hold on;   
            plot(self.currentFrame,self.meanAssignmentCost','xk');
            title('Assignment costs');
            
            subplot(self.nfish+2,1,2);
            hold on;
            set(gca,'colororderindex',1);
            plot(self.currentFrame,self.meanCost','x');
            title('mean cost');
            
            
            plotcost = nan(self.nfish,size(fcost,3));
            ssfcost = bsxfun(@times,fcost,shiftdim(scale(:),-2));
            for i = 1:size(self.assignments,1)
              plotcost(self.assignments(i,1),:) = ssfcost(self.assignments(i,1),self.assignments(i,2),:);
            end
            

          
            col = 'rgbykmc';
            for i = 1:size(plotcost,1)
              subplot(self.nfish+2,1,2+i);
              hold on;
              for j = 1:size(plotcost,2)
                h(j) = plot(self.currentFrame,plotcost(i,j),['x',col(j)]);
              end
              idx = find(i==self.assignments(:,1));
              if ~isempty(idx)
                plot(self.currentFrame,sfcost(self.assignments(idx,1),self.assignments(idx,2)),'ok');
              end
              ylim([0,self.meanAssignmentCost*5])
            end
    % $$$           if self.currentFrame==5
    % $$$             legend(h,self.costinfo);
    % $$$           end
            

          end
          drawnow;
        end

        % compute the class prob switching matrix
        if length(self.unassignedDetections)>2
          verbose('there were %d unassigned detections\r',length(self.unassignedDetections))
        end

        
      end
      
    end
    
    
    
    function updateTracks(self)

      %% update assigned
      assignments = self.assignments;
      numAssignedTracks = size(assignments, 1);
      classopts = self.opts.classifier;
      trackopts = self.opts.tracks;
      
      for i = 1:numAssignedTracks
        trackIdx = assignments(i, 1);
        detectionIdx = assignments(i, 2);
        centroid = self.centroids(detectionIdx, :);

        %% Correct the estimate of the object's location  using the new detection.
        if trackopts.kalmanFilterPredcition
          correct(self.tracks(trackIdx).predictor, [centroid]);
        end

        tau = trackopts.tauVelocity;
        self.tracks(trackIdx).velocity = self.tracks(trackIdx).velocity*(1-1/tau) ...
            + 1/tau*(self.tracks(trackIdx).centroid-centroid);
        
        %% detect whether the fish feature might be reversed
        reversed = 0;
        if self.useMex
          thickness = self.thickness(detectionIdx,:);
          centerLine = self.centerLine(:,:,detectionIdx);
          n2 = ceil(size(centerLine,1)/2);
          centeredCenterLine = bsxfun(@minus,centerLine,centerLine(n2,:));
          vel = self.tracks(trackIdx).velocity; % already updated velocity 
          % check whether centerLIne is in the right direction
          p = centeredCenterLine * vel';
          if sum(p(1:n2))>sum(p(n2:end))
            %reverse
            reversed = 1;
            centerLine = centerLine(end:-1:1,:);
            thickness = thickness(end:-1:1);

            if self.opts.detector.fixedSize %% only for saveing needed
              self.segments(detectionIdx).FilledImageFixedSizeRotated = ...
                  self.segments(detectionIdx).FilledImageFixedSizeRotated(:,end:-1:1,:);
            end
          
          end
        else
          centerLine = [];
          thickness = [];
        end
        self.tracks(trackIdx).centerLine = centerLine;
        self.tracks(trackIdx).thickness = thickness;
        self.segments(detectionIdx).reversed = reversed;

        %% update classifier 
        if ~isempty(self.classProb)
          classprob = self.classProb(detectionIdx,:);
          probNoise = self.classProbNoise(detectionIdx);
        else
          classprob = nan(1,self.nfish);
          probNoise = NaN;
        end


        if reversed  && classopts.discardPotentialReversed
          % not update.. fishfeature reversed (cannot be changed outside of mex)
          probNoise = NaN;
          classprob = nan(1,self.nfish);
        else
          self.tracks(trackIdx).classProb =  classprob; % only update tracks if not reversed
        end
        % always update history, though (weight will be zero for NaN)
        reasonable = self.tracks(trackIdx).classProbHistory.update(classprob, probNoise);
        self.tracks(trackIdx).clpMovAvg = self.tracks(trackIdx).classProbHistory.movavg;         
        
        %% Replace predicted bounding box with detected bounding box.
        self.tracks(trackIdx).bbox = self.bboxes(detectionIdx, :);
        
        %% Update track's age.
        self.tracks(trackIdx).age = self.tracks(trackIdx).age + 1;

        %% segment/features
        self.tracks(trackIdx).segment = self.segments(detectionIdx);
        for ct = self.featurecosttypes
          self.tracks(trackIdx).features.(ct{1}) = self.features.(ct{1})(detectionIdx,:);
        end

        %% save IDfeatures
        thres = classopts.crossCostThres*self.meanAssignmentCost;
        if self.cost(detectionIdx) <= thres && (reasonable  || ~self.isInitClassifier)

          self.tracks(trackIdx).batchFeatures(self.tracks(trackIdx).nextBatchIdx,:)  = self.idfeatures(detectionIdx,:);
          self.tracks(trackIdx).nextBatchIdx = ...
              min(self.tracks(trackIdx).nextBatchIdx+1,classopts.maxFramesPerBatch);
        else
          %[~,w] =self.tracks(trackIdx).classProbHistory.getData(1);
          %verbose('Above thres : w=%1.2f, cost=%1.2f (thres == %1.2f)',w,self.cost(detectionIdx),thres);
          self.tracks(trackIdx).nIdFeaturesLeftOut = self.tracks(trackIdx).nIdFeaturesLeftOut +1;
        end
        
        self.tracks(trackIdx).centroid = centroid;
        self.tracks(trackIdx).assignmentCost =  self.cost(detectionIdx);
        
        %% Update visibility.
        self.tracks(trackIdx).totalVisibleCount =  self.tracks(trackIdx).totalVisibleCount + 1;
        self.tracks(trackIdx).consequtiveInvisibleCount = 0;
      end
      

      
      
      %% update unassigned
      for i = 1:length(self.unassignedTracks)
        ind = self.unassignedTracks(i);
        self.tracks(ind).age = self.tracks(ind).age + 1;
        self.tracks(ind).consequtiveInvisibleCount =  self.tracks(ind).consequtiveInvisibleCount + 1;
      
        % do not update the tracks with nan (last classprob will remain)
        % only put into history:
        self.tracks(ind).classProbHistory.update(nan(1,self.nfish),NaN);

      end
    
    end

    
    function deleteLostTracks(self)
      if isempty(self.tracks) || ~self.isInitClassifier
        return;
      end
      
      
      % Compute the fraction of the track's age for which it was visible.
      ages = [self.tracks(:).age];
      totalVisibleCounts = [self.tracks(:).totalVisibleCount];
      visibility = totalVisibleCounts ./ ages;
      
      % Find the indices of 'lost' tracks.
      lostInds = (ages < self.opts.tracks.ageThreshold & visibility < 0.6);
      lostInds =  lostInds | [self.tracks(:).consequtiveInvisibleCount] >= self.opts.tracks.invisibleForTooLong;

      if any(lostInds)

        % Delete lost tracks.
        lostTrackIds = [self.tracks(lostInds).id];
        fishIds = [self.tracks.fishId];
        if any(~isnan(fishIds(lostInds)))
          warning('want to delete a true track');
          keyboard;
        end
          
        self.tracks = self.tracks(~lostInds);

        %delete the track ids from the crossings
        for i = 1:length(self.tracks)
          self.tracks(i).crossedTrackIds = setdiff(self.tracks(i).crossedTrackIds,lostTrackIds);
% $$$           if ~isempty(setdiff(self.tracks(i).crossedTrackIds,lostTrackIds))
% $$$             self.tracks(i).crossedTrackIds = []; % we have to delete all the
% $$$                                                  % crossings. otherwise mixup happens
% $$$           end
        end
      
        verbose('Deleted %d tracks',sum(lostInds));
      end
      
    end
    
    
    function createNewTracks(self)
    % create new tracks for unassigndetections (if less the number of fish)

      availablefishids = setdiff(1:self.nfish,cat(2,self.tracks.fishId));
      if  ~self.opts.tracks.withTrackDeletion && ~length(availablefishids)
        return;
      end
      
      
      if self.isInitClassifier
        % already have fish.  take only one detection with least
        % assignment cost to the other classes
% $$$         cost = self.cost;
% $$$         [mcost,idx] = min(self.cost(self.unassignedDetections));
% $$$         
% $$$         if mcost < self.opts.tracks.costOfNonAssignment
% $$$           unassignedDetections = self.unassignedDetections(idx);
% $$$         else
% $$$           return; % rather do nothing
% $$$         end
        unassignedDetections = self.unassignedDetections;
      else
        %take all at the beginning

        if ~length(availablefishids)
          return
        end
        unassignedDetections = self.unassignedDetections;
      end

      if isempty(unassignedDetections)
        return;
      end
      
      s =0;
      for i = unassignedDetections(:)'
        s = s+1;
        if s<=length(availablefishids)
          newfishid = availablefishids(s);
        else
          newfishid = NaN; 
          % just take the first detections. Maybe sort which to take?
          if ~self.opts.tracks.withTrackDeletion
            break
          end
        end
        
        if self.opts.tracks.kalmanFilterPredcition
          pred = self.newPredictor(self.centroids(i,:));
        else
          pred = [];
        end
        
        
        %save the features to later update the batchClassifier
        idfeat = self.idfeatures(i,:);
        bfeat = zeros(self.opts.classifier.maxFramesPerBatch,size(idfeat,2),'single');
        bfeat(1,:) = idfeat;
        if ~all(isnan(self.classProb))
          clprob = self.classProb(i,:);
          clprobnoise = self.classProbNoise(i);
        else
          clprob = nan(1,self.nfish);
          clprobnoise = nan;
        end
        
        if self.useMex
          centerLine = self.centerLine(:,:,i);
          thickness = self.thickness(i,:);
        else
          centerLine = [];
          thickness = [];
        end

        feat = [];
        for ct = self.featurecosttypes
          feat.(ct{1}) = self.features.(ct{1})(i,:);
        end
        
        
        % Create a new track.
        newTrack = struct(...
          'id',        self.nextId,       ...
          'fishId',   newfishid,...
          'bbox',      self.bboxes(i,:),  ...
          'centroid',  self.centroids(i,:),  ...
          'velocity', [0,0], ...
          'centerLine',  centerLine,...
          'thickness',  thickness,...
          'segment',   self.segments(i),  ... 
          'features',   feat,  ... 
          'classProb',clprob,...
          'classProbHistory',FishClassProbHistory(self.nfish),...
          'clpMovAvg',clprob,...
          'batchFeatures',bfeat,          ...
          'nextBatchIdx', 1,              ...
          'predictor', pred,              ...
          'firstFrameOfCrossing', self.currentFrame, ...
          'lastFrameOfCrossing', self.currentFrame, ...
          'crossedTrackIds',[],...
          'age',       1,                 ...
          'totalVisibleCount', 1,         ...
          'consequtiveInvisibleCount', 0, ...
          'assignmentCost',Inf, ...
          'nIdFeaturesLeftOut',0,...
          'stmInfo',[]);

        % update the classprobhistory and set the parameters
        %newTrack.classProbHistory.lambda = self.fishwidth/self.fishlength;
        newTrack.classProbHistory.lambda = self.opts.classifier.bendingThres;
        newTrack.classProbHistory.update(clprob,clprobnoise);
        
        % Add it to the array of tracks.
        self.tracks(end + 1) = newTrack;

        
        % Increment the next id.
        self.nextId = self.nextId + 1;
        

      end
    end
    

    function crossings = detectTrackCrossings(self)
      
      crossings = cell(1,0);
      if length(self.tracks)<2
        return
      end

      crossmat = false(length(self.tracks),length(self.tracks));
      bbox = double(cat(1,self.tracks.bbox));
      centroids = cat(1,self.tracks.centroid);
      msk = pdist2Euclidean(centroids,centroids)<self.currentCrossBoxLengthScale*self.fishlength;

      bBoxXOverlap = bsxfun(@ge,bbox(:,1) + bbox(:,3), bbox(:,1)') & bsxfun(@lt,bbox(:,1), bbox(:,1)');
      bBoxYOverlap = bsxfun(@ge,bbox(:,2) + bbox(:,4), bbox(:,2)') & bsxfun(@lt,bbox(:,2), bbox(:,2)');
      bBoxOverlap = bBoxXOverlap & bBoxYOverlap;

      costThres = self.opts.classifier.crossCostThres*self.meanAssignmentCost;
      for i = 1:length(self.tracks)
        
        if sum(msk(i,:))>1  || sum(bBoxOverlap(i,:))>1  
          if self.tracks(i).consequtiveInvisibleCount>0 || self.tracks(i).assignmentCost>costThres 
            crossmat(i,:) =  msk(i,:)  | bBoxOverlap(i,:);
          end
        end
      end
      
      %make symmetric CAUTION: MIGHT NOT BE SYMMETRIC!
      if any(crossmat(:))
        crossmat = crossmat | crossmat';
        [~,~,crossings] = networkComponents(crossmat);
        crossings(cellfun('length',crossings)==1) = []; % delete self-components
      end % else just empty (see above)
    end
    
    
    function updateTrackCrossings(self)
    % detects and updates the crossing of the track. We want to have the
    % tracks time of joining the cross and the time of the track leaving. 

      trackIds = cat(2,self.tracks.id);

      for i_cross = 1:length(self.crossings) % now cell array of trackIndices that cross
      
        crossedIndices = self.crossings{i_cross}(:)';
        newCrossedTracks = [];        
        allCrossedTracks = [];
        s = 0;
        while isempty(allCrossedTracks) || length(allCrossedTracks)~=length(newCrossedTracks) || ...
              any(sort(allCrossedTracks)~=sort(newCrossedTracks))

          newCrossedTracks = trackIds(crossedIndices);                
          allCrossedTracks = union(newCrossedTracks,cat(2,self.tracks(crossedIndices).crossedTrackIds));
          allCrossedTracks = allCrossedTracks(:)';
          crossedIndices = find(ismember(trackIds,allCrossedTracks));
          s = s+1;
          
        end
        
          
        for trackIndex = crossedIndices(:)'
          % handle each trackIndex one by one

          if isempty(self.tracks(trackIndex).crossedTrackIds)
            %new crossing. set the first frame
            self.tracks(trackIndex).firstFrameOfCrossing = self.currentFrame;            
          end
          % alway make sure that all tracks in a crossing have the
          % same Id field to avoid mixing up later.

          self.tracks(trackIndex).crossedTrackIds = allCrossedTracks;
          self.tracks(trackIndex).lastFrameOfCrossing = self.currentFrame;      

          self.resetBatchIdx(trackIndex); % we reset the batch to avoid the possible mixing of features
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

      self.updatePos();
      
      self.updateClassifier();

      self.createNewTracks();

      
      if self.opts.tracks.withTrackDeletion % CAUTION: some strange osciallations can happen...
        self.deleteLostTracks(); 
      end

    end
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % Results & Plotting helper functions
    %--------------------
    
    function savedTracks = appendSavedTracks(self,savedTracks)
      
      if self.opts.tracks.withTrackDeletion
        warning('DO NOT CONVERT SAVEDTRACKS')
        return
      end
      
      s2mat = strucarr2strucmat(savedTracks);
      if isempty(self.savedTracks)
        self.savedTracks = s2mat;
        savedTracks(:) =  [];
        return
      end
      for f = fieldnames(s2mat)'
        d = length(size(s2mat.(f{1})));% at least 3
        self.savedTracks.(f{1}) = cat(d,self.savedTracks.(f{1}),s2mat.(f{1}));
      end
      savedTracks(:) =  [];
    end
    
      
    
    function generateResults(self)

      if isempty(self.savedTracks.id) 
        verbose('WARNING: cannot generate results!')
        verbose('WARNING: not all fish detected. Maybe adjust "nfish" setting.');
        return;
      end

% $$$       if size(self.savedTracks,2)<self.nfish
% $$$         self.res.tracks(end+1,1:self.nfish) =savedTracks(1);
% $$$         self.res.tracks(end,:) = [];
% $$$       end
      
      nFrames = self.currentFrame;

      self.res.pos = permute(self.pos(:,:,1:nFrames),[3,1,2]);      

      %apply a slight smoothing to the pos
      %sz = size(self.res.pos);
      %nconv = 3;
      %self.res.pos = reshape(conv2(self.res.pos(:,:),ones(nconv,1)/nconv,'same'),sz);

      f2t = self.fishId2TrackId(1:nFrames,:)';
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
      idx = s2i(size(Loc),[Loc(:),fridx(:)]);
      for f = fieldnames(self.savedTracks)'
        if isempty(self.savedTracks.(f{1}))
          continue;
        end
        sz = size(self.savedTracks.(f{1}));
        d = length(sz); % at least 3
        trackdat = permute(self.savedTracks.(f{1}),[d,2,1,3:d-1]);
        fishdat = reshape(trackdat(idx,:),[self.nfish,nFrames,sz(2),sz(1),sz(3:d-1)]);
        self.res.tracks.(f{1}) = permute(fishdat,[2,1,3:d+1]);
      end
      
      fishid =  (1:self.nfish)' * ones(1,nFrames);
      fishid(msk) = NaN;
      self.res.tracks.fishId = fishid';

    end
    
    
    function trackinfo = getCurrentTracks(self)
    % gets the track info 
      
      trackinfo = repmat(struct,[1,self.nfish]);

      % only take those that are sure for now
      %fishIds = cat(self.tracks.fishId);
      %tracks = self.tracks(~isnan(fishIds));
      tracks = self.tracks;

      % need id 
      n =length(tracks);
      [trackinfo(1:n).id] = deal(tracks.id);
      t = self.videoHandler.currentTime;
      [trackinfo(1:n).t] = deal(t);

      if isempty(self.saveFieldSegIf)
        self.saveFieldSegIf = zeros(length(self.saveFields),1);
        for i_f = 1:length(self.saveFields)
          f = self.saveFields{i_f};
          idx = find((f=='.'));
          if ~isempty(idx)
            self.saveFieldSegIf(i_f) = idx;

            if (length(idx)>1)
              error('Cannot save 2nd level features');
            end

          end
          
          if self.saveFieldSegIf(i_f) 
            fsub = f;
            fsub(idx) = '_';
            self.saveFieldSub{i_f} = fsub;
            if ~strcmp(fsub(1:idx(1)-1),'segment') 
              error('expected "segment"');
            end
          end
        end
      end
      
      if any(self.saveFieldSegIf)
        segs = cat(1,tracks.segment);
      end
      
      for i_f = 1:length(self.saveFields)
        f = self.saveFields{i_f};
        if self.saveFieldSegIf(i_f)

          if isempty(tracks) 
            eval(sprintf('[trackinfo(:).%s] = deal([]);',self.saveFieldSub{i_f}));
          else
            [trackinfo(1:n).(self.saveFieldSub{i_f})] = deal(segs.(f(self.saveFieldSegIf(i_f)+1:end)));
          end
        else
          [trackinfo(1:n).(f)] = deal(tracks.(f));
        end
      end
      
    end
    
    
    function displayTrackingResults(self)
      cjet = jet(self.nextId);

      minVisibleCount = 2;
      tracks = self.tracks;
      uframe = getCurrentFrame(self.videoHandler);

      if ~isa(uframe,'uint8')
        uframe = uint8(uframe*255);
      end
      
      if size(uframe,3)==1
        uframe = repmat(uframe,[1,1,3]);
      end
      

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
          predictedTrackInds = [reliableTracks(:).consequtiveInvisibleCount] > 0;
          isPredicted = cell(size(labels));
          isPredicted(predictedTrackInds) = {' predicted'};

          for i = 1:length(reliableTracks)
            if ~predictedTrackInds(i) & reliableTracks(i).assignmentCost>1
              isPredicted{i} = sprintf('  %1.0d',round(reliableTracks(i).assignmentCost));
            end
          end
          

          labels = strcat(labels, isPredicted);
          cols = jet(self.nfish);
          cols_grey = [0.5,0.5,0.5;cols];
          
          
          if length(reliableTracks)==self.nfish && self.isInitClassifier

            ids = [reliableTracks(:).fishId];
            pids = arrayfun(@(x)x.classProbHistory.leakyID(),reliableTracks);
            clabels = uint8(cols(ids,:)*255);
            pclabels = uint8(cols(pids,:)*255);
            % Draw the objects on the frame.
            uframe = insertShape(uframe, 'FilledRectangle', bboxes,'Color',pclabels,'Opacity',0.2);
            uframe = insertObjectAnnotation(uframe, 'rectangle', bboxes, labels,'Color',clabels);

          else
            % Draw the objects on the frame inm grey 
            uframe = insertObjectAnnotation(uframe, 'rectangle', bboxes, labels,'Color', uint8([0.5,0.5,0.5]*255));
          end
          
          center = cat(1, reliableTracks.centroid);
          center = cat(2,center,max(self.fishlength/2,self.fishwidth)*ones(size(center,1),1));
          
          crossedTrackIdStrs = arrayfun(@(x)num2str(sort(x.crossedTrackIds)), reliableTracks,'uni',0);
          [u,idxct,idxu] = unique(crossedTrackIdStrs);
          for i =1:length(u)
            crossId = reliableTracks(idxct(i)).crossedTrackIds;
            [~,cross] =  ismember(crossId,[reliableTracks.id]);
            cross(~cross) = [];
            
            uframe = insertObjectAnnotation(uframe, 'circle', center(cross,:), 'Crossing!','Color', uint8(cols(i,:)*255));
          end
          
        
          %% insert some  features
          if self.displayif>1
            for i = 1:length(reliableTracks)
              
              id = reliableTracks(i).fishId;
              if isnan(id)
                id = 1;
              end
              
              segm = reliableTracks(i).segment;
              
              
              if isfield(segm,'MSERregions') && ~isempty(segm.MSERregions)
                for j = 1:length(segm.MSERregions)
                  pos = segm.MSERregionsOffset + segm.MSERregions(j).Location;
                  uframe = insertMarker(uframe, pos, '*', 'color', uint8(cols(id,:)*255), 'size', 4);
                end
              end
              pos = reliableTracks(i).centroid;
              uframe = insertMarker(uframe, pos, 'o', 'color', uint8(cols(id,:)*255), 'size', 5);
              
              if ~isempty(reliableTracks(i).centerLine) && ~any(isnan(reliableTracks(i).centerLine(:)))
                uframe = insertMarker(uframe, reliableTracks(i).centerLine, 's', 'color', uint8([1,1,1]*255), 'size', 2);
              end
              
            end
          end
          
          
            
          if self.displayif>2
            %% insert more markers
            if length(self.tracks)==self.nfish
              howmany = 25;
              idx = max(self.currentFrame-howmany,1):self.currentFrame;
              trackpos = self.pos(:,:,idx);
              f2t = self.fishId2TrackId(idx,:);
              delidx = find(any(any(isnan(trackpos),1),2));
              trackpos(:,:,delidx) = [];
              f2t(delidx,:) = [];
              cli = zeros(length(idx),self.nfish);
              for iii = 1:length(self.tracks)
                [classProbs, w]= self.tracks(iii).classProbHistory.getData(length(idx));
                [~,classIdx] = max(classProbs,[],2);
                msk = any(isnan(classProbs),2) | w < self.tracks(iii).classProbHistory.reasonableThres;
                classIdx(msk) = 0;
                cli(end-size(classIdx)+1:end,iii) = classIdx;
              end
              cli(delidx,:) = [];
              trackIds = [self.tracks.id];
              [~,f2i] = ismember(f2t,trackIds);
              cli = cat(2,ones(size(cli,1),1),cli+1);
              if ~isempty(trackpos)
                for ii = 1:self.nfish
                  idx1 = f2i(:,ii)+1;
                  inds = s2i(size(cli),[(1:size(cli,1))',idx1]);
                  inds2 = cli(inds);
                  pos1 = squeeze(trackpos(:,ii,~isnan(inds2)))';
                  cols1 = uint8(cols_grey(inds2(~isnan(inds2)),:)*255);
                  cols2 = uint8(cols(ii,:)*255);
                  uframe = insertMarker(uframe, pos1, 'o', 'color', cols2 , 'size', 3);
                  uframe = insertMarker(uframe, pos1, 'x', 'color', cols1 , 'size', 2);
                end
              end
            end
          end
          
        end
      end
      
      self.videoPlayer.step(uframe);
      self.stepVideoWriter(uframe);
      
    end

  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  PUBLIC METHODS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Constructor
    %--------------------
    function self = FishTracker(vid,varargin) 
    % FISHTRACKER(VID,...) starts a FishTracker objects on a given video file
    % VID. If VID==[] a uidialog is opened to select the video file. Run
    % 
    % Options can the set with FISHTRACKER(VID,OPTNAME1,OPTVALUE1,...) or
    % FISHTRACKER(VID,OPTS) where OPTS is a structure with fields identical
    % to the names of the options.
    %  
    % To see an help for the possible options run (without arguments):
    %  >> FishTracker
    %
    % Example:
    %  >> ft = FishTracker([],'nfish',3);
    %  >> ft.track([0,20]); 
    %  >> ft.plot();

      

      def.vid = [];
      doc.vid = 'Video file name';
      
      %% properties
      def.opts.displayif = 3;
      doc.displayif = {'0: no display, set to higher values ' ...
                        '(2-7) for more plots'};
      
      def.opts.nfish = [];
      doc.nfish = {'Number of fish. Needs to be fixed ' ...
                   '(ATTEMPT to be estimated if empty)'};
      
      def.opts.fishlength = [];
      doc.fishlength = {'Approx length of fish in pixel (will' ...
                        'be estimated if empty)'};
      def.opts.fishwidth = [];
      doc.fishwidth = {'Approx width of fish in pixel (will' ...
                        'be estimated if empty)'};

      def.opts.maxVelocity = [];
      doc.maxVelocity = {'Maximal velocity in px/sec [will be ' ...
                          'estimated if empty]'};
      
      def.opts.stimulusPresenter = [];
      doc.stimulusPresenter = {'Stimulus presenter object ','[for STMIF=true]'};
      
      def.opts.stmif = false;
      doc.stmif = 'Use online stimulation';

      def.opts.useMex = true;
      doc.useMex = {'Use the C++ version FishVideoHandler','(FAST)'};

      def.opts.useOpenCV = true;
      doc.useOpenCV = 'uses mexopencv library (FAST)';

      def.opts.useKNN = false;
      doc.useKNN = {'Use the KNN instead of thresholding ' ...
                    'for background segmentation (SLOW)'};

      def.opts.useScaledFormat = false;
      doc.useScaledFormat = 'Use adaptive scaled gray format';

      
      %% additional options      
      % detector
      def.opts.detector(1).history = 250;  %[nframes]
      doc.detector(1).history = 'Background update time [nFrames]';
      
      def.opts.detector.nskip = 5; 
      doc.detector.nskip = 'Skip frmaes for thresh. (useKNN=0)';

      def.opts.detector.inverted = false;  
      doc.detector.inverted = {'Set 1 for IR videos (white fish on ' ...
                          'dark background)'};
      
      def.opts.detector.adjustThresScale = 0.91;   
      doc.detector.adjustThresScale = {'0..1 : reduce to avoid many wrong ' ...
                          ' detections in case of the','thresholder'};

      def.opts.detector.fixedSize = 0;  
      doc.detector.fixedSize = 'Set 1 for saving the fixedSize image';

      % reader
      def.opts.reader(1).resizeif = false;
      doc.reader(1).resizeif = {'Whether to resize the frame before','tracking'};

      def.opts.reader.resizescale = 0.75; 
      doc.reader.resizescale = {'Fraction to resizing the frame for ' ...
                          'RESIZEIF = true'};
      
      
      % blob anaylser
      def.opts(1).blob(1).computeMSERthres= 2; % just for init str
      doc.blob(1).computeMSERthres =  {'When to compute MSER (SLOW; only ' ...
                          'for useMex=0)'};

      def.opts.blob.interpif = 1;
      doc.blob.interpif = {'Whether to interpolate while','debending'};
      
      def.opts.blob.colorfeature = false; 
      doc.blob.colorfeature = {'Use color information for fish','feature'};

      def.opts.blob.difffeature = true; 
      doc.blob.difffeature = {'Use background subtracted gray' ...
                          'images for fish feature'};

      % classification
      def.opts(1).classifier.crossCostThres = 2.0; 
      doc.classifier.crossCostThres = {'candidates for crossings: scales ' ...
                          'mean assignment cost'};

      def.opts.classifier.bendingThres = 4; % fishwidth/fishheight
      doc.classifier.bendingThres = {'For exclusion of features for ' ...
                          'classfication'};

      def.opts.classifier.reassignProbThres = 0.45; 
      doc.classifier.reassignProbThres = {'minimal probability for ' ...
                          'reassignment/crossing exits'};
      
 
      def.opts.classifier.discardPotentialReversed = true;
      doc.classifier.discardPotentialReversed = {'Discard potential reversed during','classification (better ?)'};

      def.opts.classifier.npca = 40; 
      doc.classifier.npca = 'Number of PCA components';

      def.opts.classifier.nlfd = 10; 
      doc.classifier.nlfd = {'Number of LFD components.','Set to 0 to turn off.'};

      def.opts.classifier.outliersif = false; 
      doc.classifier.outliersif = 'Exclude outliers before update';

      def.opts.classifier.minBatchN = 8; 
      doc.classifier.minBatchN = 'minimum sample size for batch update'; 

      def.opts.classifier.nFramesForInit = 200; 
      doc.classifier.nFramesForInit = {'Number of unique frames for ' ...
                          'initialize the classifier'};

      def.opts.classifier.nFramesAfterCrossing = 16;
      doc.classifier.nFramesAfterCrossing = {'When to check for ' ...
                          'permutations after crossings'};
      
      def.opts.classifier.nFramesForUniqueUpdate = 70;
      doc.classifier.nFramesForUniqueUpdate = {'Unique frames for update all' ...
                          ' simultaneously'};
  
      def.opts.classifier.nFramesForSingleUpdate = 250; 
      doc.classifier.nFramesForSingleUpdate = {'single fish update. Should be ' ...
                          'larger than unique.'};

      def.opts.classifier.maxFramesPerBatch = 280;
      doc.classifier.maxFramesPerBatch = {'Maximal frames for saving ' ...
                          'the fish feature'};

      % tracks
      def.opts(1).tracks.costtau = 500;
      doc.tracks.costtau = {'Time constant for computing the mean ' ...
                          'cost [nFrames]'};

      def.opts.tracks.tauVelocity = 5; 
      doc.tracks.tauVelocity = {'Time constant to compute the ' ...
                          'velocity [nFrames]'};
      
      def.opts.tracks.probThresForFish = 0.25;
      doc.tracks.probThresForFish = {'Classification probability to ' ...
                          'assume a fish feature'};
      
      def.opts.tracks.crossBoxLengthScale = 1; 
      doc.tracks.crossBoxLengthScale = {'How many times the bbox is regarded' ...
                          'as a crossing'};

      def.opts.tracks.crossBoxLengthScalePreInit = 0.75; 
      doc.tracks.crossBoxLengthScalePreInit = {'Before classifier init ' ...
                          '(reduced somewhat)'};

      def.opts.tracks.purelyDistanceBasedSwitching = false;
      doc.tracks.purelyDistanceBasedSwitching = {'Switching strictly','on min. distance.'};
      
      def.opts.tracks.displayEveryNFrame = 10;
      doc.tracks.displayEveryNFrame = {'How often to update the track ' ...
                          'display window (SLOW)'};
      
      def.opts.tracks.kalmanFilterPredcition = false; 
      doc.tracks.kalmanFilterPredcition = {'Whether to use Kalman filter ' ...
                          '(DEPRECIATED)'};
      
      def.opts.tracks.costOfNonAssignment =  2; 
      doc.tracks.costOfNonAssignment = {'Higher values force (possible  ' ...
                          'spurious) assignments', '[re. to  mean cost]'};

      def.opts.tracks.adjustCostDuringCrossing = true; 
      doc.tracks.adjustCostDuringCrossing = {'Whether to reduced non-assignment' ...
                          ' cost during crossings'};
      

      STRICT = 1; % be not too strict. Some settings of the videoHandler are not
                  % explicitly given here
      VERBOSELEVEL = 0;
      ALLFIELDNAMES = 1;
      
      parseInputs;
      if HELP;  return;end

      opts.classifier.invisibleForTooLong = 5; % regulates the uniqueUpdate

      
      %% other developmental parameters
      %lost tracks
      opts.tracks.invisibleForTooLong = 15; % on track basis. Only for deletion
      opts.tracks.ageThreshold = 10;
      opts.tracks.withTrackDeletion = false; % BUG !!! TURN OFF. maybe needed later 
      
      
      if ~exist('chooseNFish') && exist('helper','dir')
        addpath('helper');
      end
      
      if isempty(vid)
        vid = getVideoFile();
      end

      if isempty(opts.nfish) 
        chooseNFish(vid,1);
      end

      % construct object
      self = self@handle();  

      self.opts = opts;
      self.setProps();
      self.setupSystemObjects(vid); % will call setOpts
      
    end

    function setProps(self,opts)

      if nargin<2
        opts = self.opts;
      end
      
      for f = fieldnames(self.opts)'
        if ~any(strcmp(f{1},{'tracks','classifier','blob','detector','reader'}))  && isprop(self,f{1})
          self.(f{1}) = opts.(f{1});
        end
      end
    end
    
    
    function setOpts(self,opts)
      
      if nargin<2
        opts = self.opts;
      end
      
      for f = fieldnames(opts)'
        if any(strcmp(f{1},{'tracks','classifier','blob','detector','reader'}))
          self.opts(1).(f{1}) = opts.(f{1});
        end
      end
      self.videoHandler.setOpts(self.opts);
    end
    
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Main tracking loop
    %--------------------
  
    function track(self,trange,saveif,writefile)
    % TRACK(TIMERANGE). Detect objects, and track them across video frames. CAUTION:
    % Old tracking data will be overwritten
      
      if exist('writefile','var')
        self.writefile = writefile;
      end

      if ~exist('saveif','var') || isempty(saveif)
        saveif = false;
      end
      
      if exist('trange','var') 
        self.timerange = trange;
      else
        self.timerange = []; % take all
      end

      if saveif
        [fname,test] = self.getDefaultFileName();
        if test
          warning(sprintf('Filename %s will be overwritten !',fname));
        end
      end
      
      self.initTracking();

      if ~self.stmif
        verbose('Start tracking of %d fish  in the time range [%1.0f %1.0f]...',...
                self.nfish,self.timerange(1),self.timerange(2));
      else
        verbose('Start tracking of %d fish with stimulation', self.nfish);
      end
        
      localTimeReference = tic;
      localTime = toc(localTimeReference);
      s = 0;
      while  hasFrame(self.videoHandler) && ...
          (~self.displayif || (self.displayif && isOpen(self.videoPlayer))) ...
          && (~self.stmif || ~self.stimulusPresenter.isFinished(localTime))

        self.currentFrame = self.currentFrame + 1;
        s = s+1;

        self.segments = self.videoHandler.step();
        
        self.detectObjects();
        self.handleTracks();

        localTime = toc(localTimeReference);        
        if self.stmif && ~isempty(self.stimulusPresenter) 
          self.tracks = self.stimulusPresenter.step(self.tracks,self.videoHandler.frameSize,localTime);
        end


        if self.opts.tracks.withTrackDeletion
          savedTracks(1:length(self.tracks),s) = self.getCurrentTracks();
        else
          savedTracks(1:self.nfish,s) = self.getCurrentTracks();
        end
        
        
        if ~mod(s,100)
          savedTracks = self.appendSavedTracks(savedTracks);
          s = size(savedTracks,2);
        
          if ~self.displayif && ~self.stmif &&  ~mod(self.currentFrame,25000) && saveif
            % occasionally save for long videos 
            verbose('Save tracking process to disk..');
            self.save();
          end
        end
        
        
        if self.displayif
          if ~mod(self.currentFrame-1,self.opts.tracks.displayEveryNFrame)
            self.displayTrackingResults();
          end
        end
        
        
        if ~mod(self.currentFrame,40) && ~isempty(self.tracks)
          t = datevec(seconds(self.videoHandler.getCurrentTime()));
          verbose(['Currently at time %1.0fh %1.0fm %1.1fs                           \r'],t(4),t(5),t(6));
        end
      end

      self.cleanupCrossings();

      % make correct output structure
      if exist('savedTracks','var')
        self.appendSavedTracks(savedTracks);
        self.generateResults();
      end
      self.closeVideoObjects();
    
      if saveif
        % occasionally save for long videos 
        verbose('Save tracking to disk...');
        self.save();
      end

    end
    
    
    function [combinedFT varargout] = trackParallel(self,inTimeRange,tOverlap,minDurationPerWorker)
    %  FTNEW = FT.TRACKPARALLEL() will track the results in parallel and creates a
    %  new FT object with tracks combined
    % 
    %  FT.TRACKPARALLEL(..,TOVERLAP,MINDURATIONPERWORKER) sepecifies the overlap
    %  between workers in seconds and the minmal video time per worker in
    %  seconds
    % 
    % CAUTION: if not calculated in parallel (w.g. for too short data) the
    % returned handle object will be the IDENTICAL handle (only a reference to
    % the same data).  Only in case of parallel processing the returned object
    % will have a NEW handle (and thus reference new data).

      if ~exist('inTimeRange','var')
        inTimeRange = [];
      end
      self.videoHandler.timeRange = inTimeRange;
      self.timerange = self.videoHandler.timeRange;
      self.videoHandler.reset();
      
      dt = 1/self.videoHandler.frameRate;
      if ~exist('minDurationPerWorker','var')
        minDurationPerWorker = 300; % 5 minutes;
      end
      
      if ~exist('tOverlap','var') || isempty(tOverlap)
        % allow enough time for the classifer to init in the overlapping period
        tOverlap = 5*dt*max([self.opts.classifier.nFramesForInit,self.opts.classifier.nFramesForUniqueUpdate]); 
        tOverlap = max(tOverlap,self.opts.tracks.costtau*dt);
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
      pool = parpool(2);
      
      maxDuration = min(diff(self.timerange)/pool.NumWorkers,3600);
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
        ft.setupSystemObjects(ft.videoFile); % redo the videoHandler;
        ft.displayif = 0;
        ft.track(timeRanges(i,:));
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
        
      dt = 1/self.videoHandler.frameRate;
      tbinsize = max(self.fishlength/self.maxVelocity,2*dt);

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
        combinedObj.currentTime = self.timerange(1);
        combinedObj.pos = []; % not valid 
        combinedObj.fishId2TrackId = []; % not valid 
      end
        
    end
    
    
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Results & plotting 
    %--------------------

    function setDisplay(self,number)
    % FT.SETDISPLAY(VALUE) sets the amount of plotting during the tracking of the fish. Use 0
    % for no plot, 1..7 will display more and more plots (slowing the tracking time down
    % considerable). Larger values than 3 are only reasonable during debugging. the
    % VALUE can also be set as option "displayif" when creating the FISHTRACKER
    % object. 
      self.displayif = number;
    end
      
    
    function res = getTrackingResults(self,forceif)
    % returns the current results
      if isempty(self.res) || (exist('forceif','var') && forceif)
        self.generateResults(); % maybe not done yet
      end
      
      if isempty(self.res)
        error('No results available. First track()...');
      else
        res = self.res;
      
        % delete beyond border pixels
        posx = squeeze(self.res.pos(:,1,:));
        posy = squeeze(self.res.pos(:,2,:));

        sz = self.videoHandler.frameSize;
        posx(posx>sz(2) | posx<1) = NaN;
        posy(posy>sz(1) | posy<1) = NaN;

        res.pos(:,1,:) = posx;
        res.pos(:,2,:) = posy;

        assignmentCost = res.tracks.assignmentCost;
 
        if all(isnan(assignmentCost(:)))
          error('somehow only empty results');
        end

        assignmentCost(isnan(assignmentCost)) = Inf;
        assignmentThres = quantile(assignmentCost(:),0.999);
        idx = find(assignmentCost > assignmentThres);
        for i = 1:size(res.pos,2)
          posi = res.pos(:,i,:);
          posi(idx) = NaN;
          res.pos(:,i,:) = posi;
        end
      end
    end
    
    
    function smoothBlibs(self)
      
      if isempty(self.res) || (exist('forceif','var') && forceif)
        self.generateResults(); % maybe not done yet
      end
      
      if isempty(self.res)
        error('No results available. First track()...');
      end
      
      res = self.res;
    end
    
      
    function checkVideoHandler(self)
    % checks wether the videoHandler is still OK
      
      
      if ~iscell(self.videoFile) && ~exist(self.videoFile)
        verbose('Cannot find videofile: %s',self.videoFile);
        self.videoFile = getVideoFile(self.videoFile);
        self.videoHandler = [];
        self.videoHandler = self.newVideoHandler(self.videoFile,self.timerange);
      end
      

    end

    
    function playVideo(self,timerange,writefile)

      res = self.getTrackingResults();

      if hasOpenCV() && self.useOpenCV
        vr = @FishVideoReader;
      else
        vr = @FishVideoReaderMatlab;
      end
      
      videoReader = vr(self.videoFile,self.timerange);

      if isinf(videoReader.duration) 
        verbose('Cannot find videofile: %s',self.videoFile);
        self.videoFile = getVideoFile(self.videoFile);
        videoReader = vr(self.videoFile,self.timerange);
      end
      
      if ~exist('timerange','var')
        t = res.tracks.t(:,1);
        timerange = t([1,end]);
      end
      
      if ~exist('writefile','var') || isempty(writefile)
        writefile = '';
      end
      if isscalar(writefile) && writefile
        [a,b,c] = fileparts(self.videoFile);
        writefile = [a '/' b '_playVideo' c];
      end
      videoWriter = [];
      if ~isempty(writefile) 
        videoWriter = self.newVideoWriter();
      end

      if isempty(self.videoPlayer)
        self.videoPlayer = self.newVideoPlayer();
      end
      if  ~isOpen(self.videoPlayer)
        self.videoPlayer.show();
      end

      currentTime = videoReader.currentTime;      
      videoReader.timeRange = timerange;
      videoReader.reset();
      videoReader.setCurrentTime(timerange(1)); 

      t_tracks =res.tracks.t(:,1);
      tidx = find(t_tracks>=timerange(1) & t_tracks<timerange(2)); 
      t_tracks = t_tracks(tidx);
      cols = uint8(255*jet(self.nfish));
      s = 0;
      while videoReader.hasFrame() && s<length(tidx) && isOpen(self.videoPlayer)
        uframe = videoReader.readFrameFormat('RGBU');
        s = s+1;
        
        t = videoReader.currentTime;
        %assert(abs(t_tracks(s)-t)<=1/videoReader.frameRate);

        % Get bounding boxes.
        bboxes = shiftdim(res.tracks.bbox(tidx(s),:,:),1);
          
        % Get ids.

        ids = res.tracks.fishId(tidx(s),:);
        foundidx = ~isnan(ids);
        ids = int32(ids(foundidx));
        labels = cellstr(int2str(ids'));
        clabels = cols(ids,:);

        highcost = isinf(res.tracks.assignmentCost(tidx(s),:));
        highcost = highcost(foundidx);
        clabels(highcost,:) = 127; % grey
        % Draw the objects on the frame.
        uframe = insertObjectAnnotation(uframe, 'rectangle', bboxes(foundidx,:), labels,'Color',clabels);
        

        clprob = shiftdim(res.tracks.classProb(tidx(s),:,:),1);
        clprob(isnan(clprob)) = 0;
        dia = self.fishlength/self.nfish;
        for i_cl = 1:self.nfish
          for j = 1:length(foundidx)
            if ~foundidx(j)
              continue;
            end
            pos = [bboxes(j,1)+dia*(i_cl-1)+dia/2,bboxes(j,2)+ dia/2 ...
                   + bboxes(j,4),dia/2];
            uframe = insertShape(uframe, 'FilledCircle',pos, 'color', cols(i_cl,:), 'opacity',clprob(j,i_cl));
          end
        end

        self.videoPlayer.step(uframe);
        if ~isempty(videoWriter)
          videoWriter.write(uframe);
        end
        
        d = seconds(t);
        d.Format = 'hh:mm:ss';
        verbose('CurrentTime %s\r',char(d))
      end
      fprintf('\n')
      clear videoWriter videoReader
    end
    
    function [savemat,exists] = getDefaultFileName(self)
      if iscell(self.videoFile) && ~isempty(self.videoFile{2})
        vid = self.videoFile{2};
      elseif ischar(self.videoFile)
        vid =self.videoFile;
      else
        vid = '';
      end
      
      [savepath,savename] = fileparts(vid);
      if isempty(savename)
        savename = 'FishTrackerSave';
      end
      
      savemat = [savepath filesep savename '.mat'];
      exists = ~~exist(savemat,'file');
        
    end
    
    
    function  save(self,savename,savepath,vname)
    % SELF.SAVE([SAVENAME,SAVEPATH,VNAME]) saves the object

      [defname,exists] = self.getDefaultFileName();
      if ~exist('savename','var') || isempty(savename)
        [~,savename] = fileparts(defname);
      end
      if ~exist('savepath','var') || isempty(savepath)
        [savepath] = fileparts(defname);
      end
      
      if ~exist('vname','var') || isempty(vname)
        vname = 'ft';
      end
      eval([vname '=self;']);
      fname = [savepath filesep savename '.mat'];
      save([savepath filesep savename],vname,'-v7.3');

      verbose('saved variable %s to %s',vname,fname);
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
    
    function plotClassProb(self,plotTimeRange,fishIds);
      
         
      if ~exist('fishIds','var') || isempty(fishIds)
        fishIds = 1:self.nfish;
      end

      if ~exist('plotTimeRange','var') || isempty(plotTimeRange)
        plotTimeRange = self.timerange;
      end

      clf;
      res = self.getTrackingResults();
      t = res.tracks.t(:,1);
      plotidx = t>=plotTimeRange(1) & t<plotTimeRange(2);
      
      if ~isfield(res.tracks,'classProb') || isempty(res.tracks.classProb)
        error('No classprob data');
      end

      prob = res.tracks.classProb(plotidx,fishIds,fishIds);
      probid = nanmean(prob(:,1:length(fishIds)+1:end),2);
      if self.nfish<=4
        p = perms(1:self.nfish);
        probmax = probid;
        for i = 1:size(p)
          idx = sub2ind([self.nfish,self.nfish],(1:self.nfish)',p(i,:)');
          probmax = nanmax(probmax,nanmean(prob(:,idx),2));
        end
      else
        probmax = nanmax(prob(:,:,:),[],3);
      end
      
      nconv = 50;
      probid = conv(probid,ones(nconv,1)/nconv,'same');
      probmax = conv(probmax,ones(nconv,1)/nconv,'same');
            
      tt = t(plotidx);

      
      a(1) = subplot(3,1,1);
      plot([probid,probmax]);

      a(2) = subplot(3,1,2);
      seq = res.tracks.consequtiveInvisibleCount(plotidx,fishIds);
      plot(seq);

      a(3) = subplot(3,1,3);
      x = res.pos(plotidx,:,fishIds);
      v2 = squeeze(sum(abs(diff(x,1,1)),2)); 
      plot(v2);
      linkaxes(a,'x');
      
      
    end
    
    
    function plotCenterLine(self,fishIds,plotTimeRange)
      
        
      if ~exist('fishIds','var') || isempty(fishIds)
        fishIds = 1:self.nfish;
      end

      if ~exist('plotTimeRange','var') || isempty(plotTimeRange)
        plotTimeRange = self.timerange;
      end


      res = self.getTrackingResults();
      t = res.tracks.t(:,1);
      plotidx = t>=plotTimeRange(1) & t<plotTimeRange(2);
      
      if  ~isfield(res.tracks,'centerLine') || isempty(res.tracks.centerLine)
        error('No centerLine data');
      end

      clx = permute(res.tracks.centerLine(plotidx,fishIds,1,:),[4,1,2,3]);
      cly = permute(res.tracks.centerLine(plotidx,fishIds,2,:),[4,1,2,3]);
      clx = convn(clx,ones(2,1)/2,'valid');
      cly = convn(cly,ones(2,1)/2,'valid');
      lap = abs(nanmean(diff(clx,2,1))) + abs(nanmean(diff(cly,2,1),1));
      mx = nanmean(clx(:,:,:),1);
      my = nanmean(cly(:,:,:),1);

% $$$       % find sudden changes
% $$$       clf;
% $$$       %tt = t(plotidx); 
% $$$       tt = find(plotidx);% t is wrong in the txt file...
% $$$       v2 = squeeze(abs(diff(mx,1,2)) + abs(diff(my,1,2))); 
% $$$       a(1) = subplot(3,1,1);
% $$$       plot(tt(1:end-1),v2);
% $$$       title('Difference position')
% $$$       a(2) = subplot(3,1,2);
% $$$       if isfield(res.tracks,'consequtiveInvisibleCount')
% $$$         seq = res.tracks.consequtiveInvisibleCount(plotidx,fishIds);
% $$$         plot(tt,seq);
% $$$         title('Conseq. invisible counts');
% $$$       end
% $$$       
% $$$       a(3) = subplot(3,1,3);
% $$$       prob = res.tracks.classProb(plotidx,fishIds,fishIds);
% $$$       plot(tt,prob(:,1:length(fishIds)+1:end));
% $$$       title('Class prob');
% $$$       linkaxes(a,'x');
      
      
      figure;
      cmap = jet(self.nfish);
      szFrame = self.videoHandler.frameSize;

      cla;

      for i = 1:length(fishIds)
        col = cmap(fishIds(i),:);
        
        plot(clx(:,:,i),cly(:,:,i),'color',col,'linewidth',2);
        hold on;        
        plot(clx(1,:,i),cly(1,:,i),'o','color',col,'linewidth',1);

        idx = lap(1,:,i) > 1.5;
        %scatter(mx(1,idx,i),my(1,idx,i),50,'r','o','filled'); 
        plot(clx(:,idx,i),cly(:,idx,i),'color','r','linewidth',1);
      end
      colormap(gray)
    
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

    
    function plotByType(self,plottype,plotTimeRange,fishIds)
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
      t = res.tracks.t(:,1);
      plotidx = t>=plotTimeRange(1) & t<plotTimeRange(2);
      
      centroid = res.tracks.centroid;
      centroidx = centroid(plotidx,fishIds,1);
      centroidy = centroid(plotidx,fishIds,2);
      
      posx = squeeze(res.pos(plotidx,1,fishIds));
      posy = squeeze(res.pos(plotidx,2,fishIds));
      posx = conv2(posx,ones(5,1)/5,'same');
      posy = conv2(posy,ones(5,1)/5,'same');

      dt = 1/self.videoHandler.frameRate;
      
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
      

          plot(diff(posx)/dt,diff(posy)/dt,'.');
          xlabel('X-Velocity [px/sec]')
          ylabel('Y-Velocity [px/sec]')
          
          title('Velocity');
        
        case {'PROBMAP','DOMAINS'}


          dxy = 10;
          szFrame = self.videoHandler.frameSize;
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
            if size(P,3)==3
              Z = P;
            else
              Z = sum(P,3);
            end
            Z = imfilter(Z,fspecial('disk',3),'same');
            if size(P,3)==3
              image(1:szFrame(2),1:szFrame(1),Z/quantile(P(:),0.99));
            else
              imagesc(1:szFrame(2),1:szFrame(1),Z,[0,quantile(P(:),0.99)]);
            end
            
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


      

    