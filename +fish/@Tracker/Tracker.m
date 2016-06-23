classdef Tracker < handle;


  properties 
    timerange= [];    
    opts = struct([]);
  end

  properties (SetAccess = private)


    videoFile = '';
    writefile = '';
    currentFrame = 0;

    nfish = [];
    fishlength = [];
    fishwidth = [];

    maxVelocity = [];        
    avgVelocity = [];
    
    stmif = 0;
    useMex = 1;
    useOpenCV = 1;
    useScaledFormat = 0;
    useKNN = 0;
    
    verbosity = 1;
    displayif = 3;    

    
    videoHandler = [];
    videoPlayer = [];
    videoWriter = [];
    stimulusPresenter = [];
    
    fishClassifier = [];
    daGraph = [];
    tracks = [];
    res = [];

  end
  
  properties (SetAccess = private, GetAccess=private);

    saveFieldSegIf = [];
    saveFieldsIn = struct([]);
    saveFieldsOut = struct([]);

    costinfo= {'Location','Overlap','CenterLine','Classifier','Size','Area','BoundingBox'}; % bounding box ? ()
    scalecost = [];

    saveFields = {};
    
    nextId = 1; % ID of the next track
    isInitClassifier = [];

    centerLine = [];
    thickness = [];
    bboxes = [];
    centroids = [];
    idfeatures = [];
    features = [];
    cost = [];
    segments = [];
    
    classProb = [];
    classProbNoise = [];
    crossings = [];

    fishId2TrackId = [];
    meanCost = ones(4,1);
    pos = [];
    tabs = [];
    uniqueFishFrames = 0;
    
    assignments = [];
    unassignedTracks = [];
    unassignedDetections = [];
    featurecosttypes = {'Area','Size'};%,'BoundingBox'};

    leakyAvgFrame = [];
    savedTracks = [];
    savedTracksFull = [];

    currentCrossBoxLengthScale = 1; 
    
    % constants
    nFramesExtendMemory = 250;
    nFramesAppendSavedTracks = 100;
    nStepsSaveProgress = 250; % times nFramesAppendSavedTracks
    nFramesVerbose = 40;
    
    % autoinit
    nFramesForInit  = [];
    nFramesAfterCrossing = [];
    nFramesForUniqueUpdate = [];
    clpMovAvgTau = [];
    nFramesForSingleUpdate = [];
    maxFramesPerBatch = [];
    minBatchN = [];
    
    meanAssignmentCost =1;
    maxAssignmentCost =1;
    maxClassificationProb =0;
    
    avgTimeScale = [];
    
    timeStamp = [];
    lastTimeStamp = [];
    
  end

  
  methods(Static)
  
    function obj = loadobj(S)

      deleted_fields = {'saveFieldSub'};

      if isstruct(S)
        % compatibility with old FishTracker files
        opts = S.opts;
        for f = fieldnames(S.opts)'
          if ~isstruct(S.opts.(f{1})) && isfield(S,f{1})
            opts.(f{1}) = S.(f{1});
          end
        end

        opts.stmif = 0; % set to zero
        if iscell(S.videoFile)
          fish.helper.verbose('Realtime-stimulus experiment. Load grabbed video file.')
          S.videoFile = S.videoFile{2};
        end

        if ~exist(S.videoFile,'file') 
          [a,b,c] = fileparts(S.videoFile);
          if exist([b,c],'file') 
            S.videoFile = [b c];
          else
            fish.helper.verbose(['WARNING: could not find video file "%s". Please ' 'select!!'],[b c]);
            S.videoFile = '';
          end
        end

        obj = fish.Tracker(S.videoFile,'STRICT',0,opts);
        for f = fieldnames(S)'
          if isprop(obj,f{1}) 
            if ~any(strcmp(f{1},{'videoHandler','videoWriter','videoPlayer','videoFile'}))
              obj.(f{1}) = S.(f{1});
            end
          else
            if ~any(strcmp(f{1},deleted_fields))
              fish.helper.verbose('Cannot initialize field %s',f{1});
            end
          end
          
        end
      else
        obj = S;
      end
    end
  
    
    %% test
    
    function [success,t_elapsed,varargout] = runTest(tmax,opts,pathToVideo,pathToMat,plotif)
    %[SUCCESS,T_ELAPSED] = RUNTEST(TMAX,OPTS,PATHTOVIDEO,PATHTOMAT,PLOTIF) makes a validation versus
    %the idTracker [1] with the video provided by idTracker. OPTS is a struct
    %with fields given to initialize fish.Tracker. PATHTOMAT is the trajectory
    %file from idTracker tracking the same video. PATHTOVIDEO is the video file
    %from idTracker (5 Zebrafish). TMAX sets the trange of the tracking
    %comparison (in seconds). Set TMAX to [] to test the whole
    %7.4min video. T_ELAPSED is the required tracking time in
    %seconds. if PLOTIF==0 no output is given. PLOTIF==1 plots a result
    %plot and PLOTIF>1 displays the tracks currently tracked.
    %
    % [..,FTPOS,IDPOS] = RUNTEST returns additionally the position results.
    % [..,FT] = RUNTEST returns additionally the fish.Tracker object.
    %  
    % Example: 
    % >> fish.Tracker.runTest(50,[],[],[],2); % with plotting
    %
    % [1] Perez-Escudero et al Nature Methods 2014

      if nargin<1 || isempty(tmax)
        trange = [];
      else
        trange = [0,tmax];
      end
      
      if nargin<2 || isempty(opts)
        args = {};
      else
        args{1} = opts;
      end

      if nargin<3 || isempty(pathToVideo)
        if isunix()
          pathToVideo = '~/Videos/5Zebrafish_nocover_22min.avi';   
        end
        
        if ~exist('pathToVideo','var') || ~exist(pathToVideo,'file')
          pathToVideo = [fish.helper.getFTRoot() 'data' filesep '5Zebrafish_nocover_22min.avi'];   
        end
        if ~exist(pathToVideo,'file')
          error(['Please set PATHTOVIDEO to dowloaded video file available at "http://' ...
                 'www.cajal.csic.es/files/gpolavieja/5Zebrafish_nocover_22min.avi" ' ...
                 'or put the downloaded avi-file into the "+fish/../data" directory.']);
        end
      end
      
      if nargin<4 || isempty(pathToMat)
        pathToMat = [fish.helper.getFTRoot() 'data' filesep 'trajectories.mat']; 
      end

      if nargin<5
        plotif = 1;
      end
      
      MAXDISTANCE = 0;
      
      % handle ID
      id = load(pathToMat);
      idres.pos = permute(id.trajectories,[1,3,2]);

      % run benchmark
      fish.helper.verbose('Starting run test.');
      ft = fish.Tracker(pathToVideo,args{:});
      ft.setDisplay(max(plotif-1,0));
      ft.addSaveFields('firstFrameOfCrossing', 'lastFrameOfCrossing');

      tic;
      ft.track(trange);
      t_elapsed=toc;

      
      nfish = ft.nfish;
      ftres = ft.getTrackingResults();
      ftresnan = ft.getTrackingResults(1);

      dist = zeros(nfish);
      offs = size(idres.pos,1) - size(ftres.pos,1);
      
      
      for i = 1:nfish
        for j = 1:nfish
          dist(i,j) = nanmean(sqrt(sum((ftres.pos(:,:,j) - idres.pos(1:end-offs,:,i)).^2,2)));
        end
      end
      assignments = fish.helper.assignDetectionsToTracks(dist,1e3);
      
      ftpos = ftres.pos;
      ftposnan = ftresnan.pos;
      idpos = idres.pos(1:end-offs,:,assignments(assignments(:,1),2));
 
      t = ftres.tracks.t(:,1);
      d = sqrt(sum((ftpos(:,:,:) - idpos(:,:,:)).^2,2))/ft.fishlength;
      
      success = mean(nanmax(d,[],3)>1)<0.05;

      varargout{1} = ftpos;
      varargout{2} = idpos;
      varargout{3} = ft;
      
      if plotif
        figure;
        
        r1 = 4;
        r2 = 1;
        s = 0;
        a = [];
        

        %distances
        s = s+1;
        a(end+1) = subplot(r1,r2,s,'align');
        

        if MAXDISTANCE
          plot(t,nanmax(d,[],3));
        else
          fish.helper.errorbarpatch(t,nanmean(d,3),fish.helper.stderr(d,3));
        end
        
        hold on;
        plot(t([1,end]),[1,1],'r:');
        
        set(a(s),'fontsize',8);
        if MAXDISTANCE
          ylabel(sprintf('Max distance\n[Fish length]'),'fontsize',10);
        else
          ylabel(sprintf('Avg. distance\n[Fish length]'),'fontsize',10);
        end
        
        xlim(t([1,end]));
        box off;
        
        %detections
        s = s+1;
        a(end+1) = subplot(r1,r2,s,'align');
        nconv = 75;
        ftnan = conv(sum(isnan(ftposnan(:,1,:)),3),ones(nconv,1)/nconv,'same');
        idnan = conv(sum(isnan(idpos(:,1,:)),3),ones(nconv,1)/nconv,'same');
        plot(t,[ftnan,idnan]*100/ft.nfish);
        xlim(t([1,end]))
        set(a(s),'fontsize',8);
        
        ylabel(sprintf('Avg. lost\n tracks [%%]'),'fontsize',10);
        
        legend({'fish.Tracker','idTracker'},'location','NorthWest','fontsize',8)
        box off;
        
        
        %crossings
        s = s+1;
        a(end+1) = subplot(r1,r2,s,'align');
        
        
        cross = diff(ftres.tracks.lastFrameOfCrossing-ftres.tracks.firstFrameOfCrossing)>0;
        ccross = conv2(sum(double(cross),2),ones(nconv,1)/nconv,'same');
        plot(t(1:end-1),ccross,'k')
        set(a(s),'fontsize',8);
        ylabel(sprintf('Avg. # of \ncrossing events'),'fontsize',10);
        xlim(t([1,end]));
        box off;
        
        
        % probability
        s = s+1;
        a(end+1) = subplot(r1,r2,s,'align');
        clp = ftres.tracks.classProb;
        dclp = zeros(size(clp,1),size(clp,2));
        mclp = zeros(size(clp,1),size(clp,2));
        for i =1:size(clp,2)
          mclp(:,i) = max(clp(:,i,setdiff(1:size(clp,2),i)),[],3);
          dclp(:,i) = clp(:,i,i);
        end
        
        percent = sum(sum(mclp>dclp))/numel(mclp);
        %fprintf('%1.1f%% correct.\n',(1-percent)*100);
        
        mclp = conv2(mclp,ones(nconv,1)/nconv,'same');
        dclp = conv2(dclp,ones(nconv,1)/nconv,'same');
        
        h1 = fish.helper.errorbarpatch(t,mean(dclp,2),fish.helper.stderr(mclp,2),'color',[0,0.5,0]);
        hold on;
        h2 = fish.helper.errorbarpatch(t,mean(mclp,2),fish.helper.stderr(dclp,2),'color',[1.0000, 0.6445, 0]);
        xlim(t([1,end]));
        box off;
        ylabel(sprintf('Class\nprobability'),'fontsize',10);
        legend([h1(1),h2(1)],{'Avg.','Other'},'location','NorthWest','fontsize',8)
        set(a(end),'fontsize',8)
        
        
        xlabel('Time [sec]','fontsize',10);
        b = fish.helper.labelsubplot(gcf);
        fish.helper.shiftaxes(b,[0.02])
        drawnow;
      end

      if nargout>1
        if ~success
          fish.helper.verbose('WARNING: RunTest failed due to high inaccuracies in the tracking.');
        end
      else
        assert(success,'RunTest failed due to high inaccuracies in the tracking.');
      end
    end

    
    function [screenBoundingBox] = calibrateStimulusScreen(camIdx,screenIdx,plotif);
    % [SCREENBOUNDINGBOX] =
    % CALIBRATESTIMULUSSCREEN(CAMIDX,SCREENIDX gets the
    % SCREENBOUNDINGBOX used for stimulus presentation.
    
      opts = [];
      opts.nfish = 4;
      opts.fishlength = 100;
      opts.fishwidth = 30;
      opts.stmif = 1;
      opts.stimulus.presenter = 'fish.stimulus.PresenterCalibration';
      opts.stimulus.screen =   screenIdx;
      opts.classifier.timeToInit = Inf; 
      opts.detector.inverted = 1;
      opts.tracks.useDagResults = 0;
      
      ft = fish.Tracker({camIdx,''},opts);

      if nargin<3
        plotif = 1;
      end
      
      if plotif>1
        ft.setDisplay(1);
        ft.setDisplay('tracks',1);      
      else
        ft.setDisplay(0);
      end
      
      ft.addSaveFields('bbox');

      ft.stimulusPresenter.width = 50; % maybe needs to be adjusted
      ft.stimulusPresenter.tmax = 30;
      ft.stimulusPresenter.freq = 0.2;
      
      fish.helper.verbose(['\n\n****** \tStarting calibration. ' ...
                          'Make sure that the IR filter is NOT ' ...
                          'installed! \n\tHit Enter to proceed!\n']);
      pause;
      
      % start detecting
      ft.track(); % track the markings
      ft.stimulusPresenter.flip(); % turn stim off

      pos = ft.deleteInvisible('pos');
      bbox = ft.deleteInvisible('bbox');

      [screenBoundingBox,xyframe] = getCalibrationBox(pos,bbox,ft.videoHandler.frameSize);
      
      if plotif
        figure;
        imagesc(ft.videoHandler.getCurrentFrame());
        colormap('gray')
        hold on;
        plot(xyframe(:,:,1),xyframe(:,:,2),'linewidth',1);
        rectangle('position',screenBoundingBox,'linewidth',1, ...
                  'edgecolor','r','facecolor','none')
        title('Estimated Screen size');
      end

      fish.helper.verbose('Found Bounding Box [%d,%d,%d,%d]',round(screenBoundingBox));

    end
    
    
  end % STATIC METHODS
  

  methods (Access=private)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % external methods
    
    
    function setProps(self,opts)
      
      if nargin<2
        opts = self.opts;
      end
      
      for f = fieldnames(self.opts)'
        if ~any(strcmp(f{1},{'cost','tracks','classifier','blob','detector','reader','stimulus'}))  && isprop(self,f{1})
          self.(f{1}) = opts.(f{1});
        end
      end
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % Gereral & init functions   
    %---------------------------------------------------

    function closeVideoObjects(self)

      if ~isempty(self.videoWriter)
        % close writer
        if ~fish.helper.hasOpenCV() || ~self.useOpenCV
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
        if fish.helper.hasOpenCV() && self.useOpenCV
          self.videoWriter.write(frame);
        else
          self.videoWriter.step(frame);
        end
      end
    end
    
    
    function [writer] = newVideoWriter(self,vidname)
      fish.helper.verbose('Init VideoWriter "%s"',vidname);
      if fish.helper.hasOpenCV() && self.useOpenCV
        writer = cv.VideoWriter(vidname,self.videoHandler.frameSize([2,1]),'FPS',self.videoHandler.frameRate,'fourcc','X264');
      else
        writer = vision.VideoFileWriter(vidname,'FrameRate',self.videoHandler.frameRate);
      end
    end
    
    
    function [player] = newVideoPlayer(self,vidname)
      player = vision.VideoPlayer();
      %player = fish.core.FishVideoPlayer();
    end
    
    
    function [handler, timerange] = newVideoHandler(self,vidname,timerange,opts)
    % note: does not set the current time. 
      if ~exist('timerange','var')
        timerange = [];
      end
      if ~exist('opts','var')
        opts = self.opts;
      end
      if self.useMex && fish.core.FishVideoHandlerMex.installed()
        self.useMex = 1;
        self.useOpenCV = 1;
        handler = fish.core.FishVideoHandlerMex(vidname,timerange,self.useKNN,opts);
      elseif fish.helper.hasOpenCV() && self.useOpenCV
        self.useMex = 0;
        self.useOpenCV = 1;
        handler = fish.core.FishVideoHandler(vidname,timerange,self.useKNN,opts);
      else
        self.useMex = 0;
        self.useOpenCV = 0;
        handler = fish.core.FishVideoHandlerMatlab(vidname,timerange,self.useKNN,opts);
      end
      timerange = handler.timeRange;

    end

    
    function setDisplayType(self)
      
      
      if self.displayif >1
        self.opts.display.level = self.displayif;
      end
      dopts = self.opts.display;
      
      if isempty(self.writefile) 

% $$$         if dopts.level==4
% $$$           self.opts.display.displayEveryNFrame = 1;
% $$$         end
% $$$         
% $$$         if dopts.level==3
% $$$           self.opts.display.displayEveryNFrame = 25;
% $$$         end
% $$$         
% $$$         if dopts.level==2
% $$$           self.opts.display.displayEveryNFrame = 25;
% $$$         end
% $$$         
% $$$         if dopts.level==1
% $$$           self.opts.display.displayEveryNFrame = 100;
% $$$         end
      
      else
        self.opts.display.displayEveryNFrame = 1;
      end
      
      if self.displayif && isempty(self.videoPlayer) && dopts.tracks
        self.videoPlayer = self.newVideoPlayer();
      end
      
      %% setting
      self.fishClassifier.plotif = ~~dopts.classifier && self.displayif;
      self.videoHandler.plotting(~~dopts.videoHandler && self.displayif);
                
      if self.stmif
        self.stimulusPresenter.progressBar = self.displayif && ...
            self.opts.display.stimulusProgress;
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
      self.removeSaveFields();
      
      if self.stmif
        if isempty(self.opts.stimulus.presenter) 
          self.stimulusPresenter = fish.stimulus.Presenter(self.opts.stimulus);
        else
          if ischar(self.opts.stimulus.presenter)
            self.stimulusPresenter = eval('%s(self.opts.stimulus)',self.opts.stimulus.presenter);
          elseif isa(self.opts.stimulus.presenter,'fish.stimulus.Presenter')
            self.stimulusPresenter = self.opts.stimulus.presenter;
          else
            error('Expect fish.stimulus.Presenter object');
          end
        end
        self.stimulusPresenter.setOpts(self.opts.stimulus);
        self.stimulusPresenter.init();
        self.addSaveFields('stmInfo','fishId','predFishId');
      end
      
      %% set all options
      self.checkOpts();
      self.setOpts();
      
    end

    
    function initTracking(self)
    % Inits the objects to make it ready to start the tracking
      
      if self.opts.tracks.withTrackDeletion
        error('Track deletion & handling tracks with DAG is currently not supported.')
      end
      
      global VERBOSELEVEL;
      VERBOSELEVEL = self.verbosity;
      
      
      self.videoHandler.timeRange = self.timerange;                       
      self.timerange = self.videoHandler.timeRange; % to get the boundaries right;
      self.videoHandler.setCurrentTime(self.timerange(1));
      self.videoHandler.fishlength = self.fishlength;
      self.videoHandler.fishwidth = self.fishwidth;
      
      self.timeStamp = self.videoHandler.getCurrentTime();
      
      if isempty(self.maxVelocity)
        self.maxVelocity = 15*self.avgVelocity;
      end
      
      self.avgTimeScale = 1/(self.avgVelocity/self.videoHandler.frameRate);
      fish.helper.verbose('Set time scale to %1.2f [frame/BL]',self.avgTimeScale);
      
      % initialize graph
      self.daGraph = [];
      self.daGraph = fish.core.FishDAGraph(self.nfish,self.nfish);

      self.checkOpts();
      self.setOpts();
      self.videoHandler.initialize();
      
      if isscalar(self.writefile) && self.writefile
        [a,b,c] = fileparts(self.videoFile);
        self.writefile = [a filesep b '_trackingVideo' c];
      end

      if ~isempty(self.writefile) 
        self.videoWriter = self.newVideoWriter(self.writefile);
        self.opts.display.displayEveryNFrame = 1;
        self.opts.display.tracks = true;
        self.displayif = 1;
      else
        if ~isempty(self.writefile)
          fish.helper.verbose('WARNING: display tracks for writing a video file!');
        end
        self.videoWriter = [];
      end

      %% get new fish classifier 
      self.fishClassifier = newFishClassifier(self,self.opts.classifier);
      self.isInitClassifier = isInit(self.fishClassifier);

      
      %% set the display  (before stimPresenter.reset)
      self.setDisplayType();     
      
      
      %% init stimulus
      if self.stmif
        if isempty(self.stimulusPresenter)
          error('No stimulus Presenter given');
        end
        self.stimulusPresenter.reset();
      end
      
      
      if self.displayif && self.opts.display.tracks && ~isOpen(self.videoPlayer)
        self.videoPlayer.show();
      end

      self.fishId2TrackId = nan(250,self.nfish);
      self.pos = [];
      self.res = [];
      self.tracks = [];
      self.savedTracks = [];
      self.savedTracksFull = [];
      self.leakyAvgFrame = [];
      

      if self.opts.tracks.keepFullTrackStruc
        fish.helper.verbose('WARNING: Will keep full tracks-structure. This will cost huge amount of memory!');
      end

      self.tracks = self.initializeTracks(); % Create an empty array of tracks.
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
      self.saveFieldsIn = struct([]);
      self.saveFieldsOut = struct([]);
      
      self.costinfo = fieldnames(self.opts.cost)';
      self.scalecost = zeros(1,length(self.costinfo));
      for i =1:length(self.costinfo)
        self.scalecost(i) = self.opts.cost.(self.costinfo{i});
      end
      delidx = find(~self.scalecost);
      self.costinfo(delidx) = [];
      self.scalecost(delidx) = [];
      assert(length(self.scalecost)>0);
      
      self.meanCost = ones(1,length(self.scalecost));
      
      self.assignments = [];
      self.unassignedTracks = [];
      self.unassignedDetections = [];
      self.cost = [];
      self.currentFrame = 0;
      
      % set time scale dependent pars
      self.currentCrossBoxLengthScale = self.opts.tracks.crossBoxLengthScalePreInit;
      
      self.nFramesAfterCrossing = min(max(ceil(self.opts.classifier.timeAfterCrossing*self.avgTimeScale),1),100);
      self.nFramesForUniqueUpdate = min(max(ceil(self.opts.classifier.timeForUniqueUpdate*self.avgTimeScale),1),250);
      self.clpMovAvgTau = max(ceil(self.opts.classifier.clpMovAvgTau*self.avgTimeScale),1);

      self.minBatchN = min(max(ceil(self.nFramesAfterCrossing*0.75),4),50);  
      self.nFramesForSingleUpdate = min(3*self.nFramesForUniqueUpdate,1000); 
      self.maxFramesPerBatch = self.nFramesForSingleUpdate + 50;

      if ~isinf(self.opts.classifier.timeToInit)
        self.nFramesForInit = min(max(ceil(self.opts.classifier.timeToInit*self.avgTimeScale),1),200);
      else
        self.nFramesForInit = Inf; % disables the classifier
      end
      
    end
    
    
    
    
    
    function detectObjects(self)
    % calls the detectors

      segm = self.segments;
      
      
      if ~isempty(segm)

        % check for minimal distance. 
        self.centroids = cat(1,segm.Centroid);
        self.bboxes = double(cat(1,segm.BoundingBox));
        
        
        overlap = getBBoxOverlap(self.bboxes,self.bboxes,2);
        overlap(1:length(segm)+1:end) = 0;
        [i,j] =  find(overlap);
        
        while ~isempty(i) && length(self.segments)>self.nfish

          if segm(i(1)).Area > segm(j(1)).Area
            deli = j(1);
          else
            deli = i(1);
          end
          
          if self.displayif && self.opts.display.detectedObjects
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
          overlap = getBBoxOverlap(self.bboxes,self.bboxes,self.fishwidth);
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
      cla = fish.core.FishBatchClassifier(self.nfish,featdim);
      
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

    
    function [nc,ncidx] = connectedComponents(self,trackIndices,assignedFishIds)
      oldFishIds = [self.tracks(trackIndices).fishId];
      A = zeros(self.nfish);
      A(sub2ind(size(A),oldFishIds,assignedFishIds)) = 1;
      [~,~,nc] = fish.helper.networkComponents(A | A');
      % delete self-components
      nc(cellfun('length',nc)==1) = []; 
      
      % recode into the trackIndex
      ncidx = cell(size(nc));
      for i = 1:length(nc)
        ncidx{i} = find(ismember(oldFishIds,nc{i}));
      end
      
      % add those single components which are actually present
      idx = find(oldFishIds == assignedFishIds);
      for i = 1:length(idx)
        ncidx{end+1} = idx(i);
        nc{end+1} = oldFishIds(idx(i));
      end
      
    end
    
        
    
    function resetBatchIdx(self,trackidx)
      [self.tracks(trackidx).nextBatchIdx] = deal(1);
    end

    
    function changed = fishClassifierUpdate(self,trackIndices,updateif)
    % updates and tests the collected features according to the current fishIDs. It does NOT
    % force and update, update is only done if the batch is likly to come from the true
    % fish. In case a valid permutation is detected and reassignProb is higher than a
    % threshold, fishIDs are exchanged.
    %
    % resets the batchIdxs after update
      changed = 0;
      fishIds = cat(2,self.tracks(trackIndices).fishId);

      % test on the whole set of possible fishIds
      [assignedFishIds prob steps probdiag] = self.predictFish(trackIndices,1:self.nfish,...
                                                        self.nFramesForSingleUpdate);
      same = assignedFishIds==fishIds; 

      if updateif && (all(same) || self.enoughEvidenceForForcedUpdate(prob,steps,probdiag))
        % do not update if classes might be mixed up (wait for more data)
        % change classifier
        feat = self.getFeatureDataFromTracks(trackIndices);
        [assignedFishIds prob] = self.fishClassifier.batchUpdate(fishIds,feat,1); % force
        self.resetBatchIdx(trackIndices);
        self.uniqueFishFrames = 0;
        changed  = 1;
        return
      end
      
      %we here also check whether some unaccounted switching occurred (maybe outside of
      %a crossing)
      if length(same)>1 && ~all(same) && ~self.opts.classifier.onlyDAGMethod
        % we have some mixed classes 

        % might be some NaNs...
        msk = ~isnan(assignedFishIds);
        trackIndices = trackIndices(msk);
        assignedFishIds = assignedFishIds(msk);
        prob = prob(msk);
        probdiag = probdiag(msk);
        steps = steps(msk);
        fishIds = cat(2,self.tracks(trackIndices).fishId);
        
        [nc,ncidx] = connectedComponents(self,trackIndices,assignedFishIds);
        
        for i = 1:length(ncidx)
          if length(ncidx{i})==1
            continue;
          end
          idx = ncidx{i};
          if all(ismember(nc{i},fishIds)) ...
              && (self.enoughEvidenceForReassignment(prob(idx),steps(idx),probdiag(idx)))

            % only permute if true permutation and enough evidence
            fish.helper.verbose('Switching fish with fishClassifierUpdate [%1.2f,%d]',min(prob(idx)),min(steps(idx)))

            switchFish(self,trackIndices(idx),assignedFishIds(idx),0);
            changed = 1;
            self.uniqueFishFrames = 0;
          end
        end
      end
      
    end

    
    function updatePos(self)
    % save current track->fishID
      t = self.currentFrame;
      fishIds = cat(1,self.tracks.fishId);
      tracks = self.tracks(~isnan(fishIds));

      if size(self.fishId2TrackId,1)< t
        self.fishId2TrackId = cat(1,self.fishId2TrackId,nan(self.nFramesExtendMemory,self.nfish));
      end
      % save track positions
      if isempty(self.pos)
        self.pos = nan(2,self.nfish,self.nFramesExtendMemory);
        self.tabs = nan(self.nFramesExtendMemory,1);
      end

      if size(self.pos,3)<t
        self.pos = cat(3,self.pos,nan(2,self.nfish,self.nFramesExtendMemory));
      end

      if size(self.tabs,1)<t
        self.tabs = cat(1,self.tabs,nan(self.nFramesExtendMemory,1));
      end

      self.tabs(t,1) = self.timeStamp;
      
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
      thres = 3;
      if length(self.tracks)<self.nfish 
        if  self.currentFrame > thres*self.nFramesForInit
          fish.helper.verbose(['nfish setting might be wrong or some fish are lost\r'])
        end
        return
      end

      if self.uniqueFishFrames>self.nFramesForInit  || ...
          self.currentFrame> thres*self.nFramesForInit % cannot wait for ages..
        
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
        
        % initialize DAG
        self.daGraph.reset(self.pos(:,:,self.currentFrame),self.currentFrame);
        
        
        self.uniqueFishFrames  = 0;
        success = true;
        
      elseif self.currentFrame < thres/2*self.nFramesForInit % cannot wait for ages..
        if ~isempty(self.crossings)
          self.resetBatchIdx(1:length(self.tracks));% we do not want a mixture at the beginning
        end
      end
      
      self.isInitClassifier = success;
      
    end
    
    function  [assignedFishIds, prob, steps, probdiag] = predictFish(self,trackIndices,fishIdSet,maxSteps);
    % predicts the fish according to the running probhist mean 
      
      C = zeros(length(trackIndices),length(fishIdSet));

      assignedFishIds = nan(size(trackIndices));
      prob = nan(size(assignedFishIds));
      probdiag = nan(size(assignedFishIds));
      steps = nan(size(assignedFishIds));

      oldFishIds = cat(1,self.tracks(trackIndices).fishId);
      for i = 1:length(trackIndices)
        track = self.tracks(trackIndices(i));
        if ~isempty(track.classProbHistory)
          nsteps = min(max(self.currentFrame-track.lastFrameOfCrossing,1),maxSteps);
          [p nsteps] = track.classProbHistory.mean(nsteps);
          steps(i) = nsteps;        
        else
          p = track.clpMovAvg;
          steps(i) = maxSteps; % fake;
        end

        C(i,:) = 1-p(fishIdSet);
        probdiag(i) = p(oldFishIds(i));
      end
      C(isnan(C)) = Inf;
      
      assignments = fish.helper.assignDetectionsToTracks(C',2);
      
      assignedFishIds(assignments(:,2)) = fishIdSet(assignments(:,1)); % fishIDs have 1:nfish order
      for i = 1:size(assignments,1)
        prob(assignments(i,2)) = 1-C(assignments(i,2),assignments(i,1));
      end
      prob(isinf(prob)) = nan;
    
    end

    
    function bool = enoughEvidenceForReassignment(self,prob,steps,probdiag)
      if any(isnan(probdiag))
        bool = false;
      else
        bool = mean(prob-probdiag)>self.opts.classifier.reassignProbThres*self.maxClassificationProb && min(steps)>=self.minBatchN;

        if bool
          fish.helper.verbose('Enough evidence: %1.2f',sum(prob-probdiag));
        end
        %bool = min(prob)>self.opts.classifier.reassignProbThres && min(steps)>=self.minBatchN;
      end
    end
    
    function bool = enoughEvidenceForForcedUpdate(self,prob,steps,probdiag)
      if any(isnan(probdiag)) || length(prob)~=self.nfish
        bool = false;
      else
        bool = min(steps) >= self.nFramesAfterCrossing;
        bool = bool && max(prob)<=self.opts.classifier.forcedReassignProbThres*self.maxClassificationProb;
      end
    end
    
    function bool = enoughEvidenceForBeingHandled(self,prob,steps,probdiag)
      if any(isnan(probdiag))
        bool = false;
      else
        bool = (min(prob)>self.maxClassificationProb*self.opts.classifier.handledProbThres ...
                && min(steps)>=self.nFramesAfterCrossing) ... 
               || max(steps)>=self.nFramesForUniqueUpdate; % to avoid very long crossings
      end
    end
    
    
    function handlePreviouslyCrossedTracks(self)
    % expect that the trackIndices are not currently crossing (they
    % can have just reached a newCrossing though).
      
      trackIndices = find(self.testBeyondBound() & ~self.testHandled());

      if isempty(trackIndices)
        return;
      end
      
      if self.displayif && self.opts.display.crossings
        plotCrossings_(self,1,trackIndices);
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

        %this should no happen because fishIDs swapping should be blocked during crossing
        assert(all(ismember(assumedFishIds,crossedFishIds)));
        

        % we force the permutation to be in the valid fish only (otherwise too many errors for many fish)
        [assignedFishIds prob steps probdiag] = self.predictFish(thisTrackIndices,crossedFishIds,...
                                                          self.nFramesForUniqueUpdate); 

        if  all(ismember(assignedFishIds,assumedFishIds)) 
          % we have a valid assignment (note that it is always a permutation inside
          % those that cross because that is guaranteed in predictFish)

          
          OLD = 0;
          if OLD
            [nc ncidx] = self.connectedComponents(thisTrackIndices,assignedFishIds);
            
            for i_nc = 1:length(ncidx)
              idx = ncidx{i_nc};
              
              if length(idx)==1 
                % identity 
                if  self.enoughEvidenceForBeingHandled(prob(idx),steps(idx),probdiag(idx)) ...
                    || length(self.tracks(thisTrackIndices(idx)).crossedTrackIds)<=1
                  % delete from crossings
                  subDeleteCrossedTrackIds(thisTrackIndices(idx));
                end
              else
                % permutation
                if self.enoughEvidenceForReassignment(prob(idx),steps(idx),probdiag(idx))
                  % enough evidence switch tracks and delete from others. 
                  fish.helper.verbose('valid permutation [%1.2f,%d]: switch...',min(prob(idx)),min(steps(idx)))
                  switchFish(self,thisTrackIndices(idx),assignedFishIds(idx),true)
                  subDeleteCrossedTrackIds(thisTrackIndices(idx)); % always handle
                  
                  
                end
              end
            end
          else
            if any(assignedFishIds ~= assumedFishIds)
              
              if self.enoughEvidenceForReassignment(prob,steps,probdiag)
                % enough evidence switch tracks and delete from others. 
                fish.helper.verbose('valid permutation [%1.2f,%d]: switch...',min(prob),min(steps))

                switchFish(self,thisTrackIndices,assignedFishIds,true)
                subDeleteCrossedTrackIds(thisTrackIndices); % always handle;
                
              end
            
            else
              if  self.enoughEvidenceForBeingHandled(prob,steps,probdiag)
                subDeleteCrossedTrackIds(thisTrackIndices);
              end
              
              for idx=1:length(thisTrackIndices)
                if  length(self.tracks(thisTrackIndices(idx)).crossedTrackIds)<=1
                  % delete all from crossings
                  subDeleteCrossedTrackIds(thisTrackIndices(idx));
                  break
                end
              end
            end
          end
          
          
            
        else
          % Still inside the fishID previously involved in the crossing. 
          % We better put
          % the track back into the crossing
          if ~any(isnan(prob))
            fish.helper.verbose('Put tracks back into crossings!')
            fishIds = [self.tracks.fishId];
            self.crossings{end +1 } = find(ismember(fishIds,crossedFishIds));
          end
        end

      end
      
      
      
      if self.displayif && self.opts.display.crossings
        plotCrossings_(self,2,trackIndices);
      end

      
      function subDeleteCrossedTrackIds(localTrackIndices)
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
      
    function computeLeakyAvgFrame(self,frame);
      
      if isempty(self.leakyAvgFrame)
        self.leakyAvgFrame = im2double(frame);
      else
        p = self.opts.display.leakyAvgFrameTau;
        self.leakyAvgFrame = (1-1/p)*self.leakyAvgFrame + im2double(frame);
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
      hitBound = self.currentFrame-[self.tracks.lastFrameOfCrossing] >= self.nFramesAfterCrossing;
    end
    
    function newCrossing = testNewCrossing(self)
      newCrossing  = ismember(1:length(self.tracks), cat(2,self.crossings{:}));
    end
    
    function [bool,pdiff] = testTrackMisalignment(self,trackIndices);
      fishIds = cat(1,self.tracks(trackIndices).fishId);
      clp = cat(1,self.tracks(trackIndices).clpMovAvg);
      idx = fish.helper.s2i(size(clp),[(1:size(clp,1))',fishIds]);
      prob_correct = mean(clp(idx));
      clp(idx) = 0;
      prob_notcorrect = mean(max(clp,[],2)); % might be Nan
      pdiff = prob_notcorrect - prob_correct ;
      bool = pdiff > self.opts.classifier.reassignProbThres; % will be zero if pdiff should be NaN
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

      if ~self.opts.classifier.onlyDAGMethod
        
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
            
            fish.helper.verbose('Probably misaligned (%1.2f)..\r',pdiff);
            self.fishClassifierUpdate(handledIndices,0); % use function for testing/switching 
            
          end
        end

      else
        % only DAG
        allhandled = true ; 
        handled = true;
      end
      
          
      
      %% unique fish update
      if (self.uniqueFishFrames-1 > self.nFramesForUniqueUpdate)  && allhandled
        if self.fishClassifierUpdate(1:length(self.tracks),1);
          fish.helper.verbose('Performed unique frames update.\r')
        end
        self.uniqueFishFrames = 0; % anyway reset;
      end

      
      %% single fish update
      enoughData =  [self.tracks.nextBatchIdx] > self.nFramesForSingleUpdate;
      singleUpdateIdx = find(handled & enoughData);
      if ~isempty(singleUpdateIdx)
        if ~self.fishClassifierUpdate(singleUpdateIdx,1) 
          if  allhandled
            self.fishClassifierUpdate(1:length(self.tracks),0); % something might be wrong. Try all
          end
          self.resetBatchIdx(singleUpdateIdx); % anyway reset          
        else
          fish.helper.verbose('Performed single fish update.\r');          
        end
      end

    end

    
    function cleanup(self) 
    %switch if necessary
      
      trackIndices = 1:length(self.tracks);
      if ~isempty(trackIndices)
        self.fishClassifierUpdate(trackIndices,0);
      end

    end
    

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % Track related functions
    %---------------------------------------------------
    
    
    function tracks = initializeTracks(self)
    % create an empty array of tracks
      tracks = struct(...
        'id', {}, ...
        'fishId', {}, ...
        'predFishId', {}, ...
        'bbox', {}, ...
        'centroid',{},...
        'velocity',{},...
        'centerLine',{},...
        'thickness',{},...
        'segment', {},...
        'features',{},...
        'classProb',{},...
        'classProbW',{},...
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
        'consecutiveInvisibleCount', {},...
        'switchedFrames',{},...
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
        predictedCentroid(1:2) = predictedCentroid(1:2) - bbox(3:4) / 2;
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
        maxDist = self.maxVelocity*(self.timeStamp-self.lastTimeStamp)*self.fishlength; % dist per frame
    % $$$         minDist = 2;
    % $$$         meanDist = sqrt(sum(cat(1,self.tracks.velocity).^2,2)); % track.velocity is in per frame
    % $$$         %nvis = max(cat(1,self.tracks.consecutiveInvisibleCount)-self.opts.tracks.invisibleForTooLong,0);
    % $$$         nvis = cat(1,self.tracks.consecutiveInvisibleCount);
    % $$$         thresDist = min(max(meanDist,minDist),maxDist).*(1+nvis);
        thresDist = maxDist;

        % Compute the cost of assigning each detection to each track.
        for k = 1:length(self.costinfo)
          switch lower(self.costinfo{k})

            case {'location','centroid'}

              %% center position 
              centers = cat(1,self.tracks.centroid);
              dst = fish.helper.pdist2Euclidean(centers,self.centroids);
              %dst = bsxfun(@rdivide,dst,thresDist).^2;

              dst = (dst/thresDist);%.^2;
              %  dst(dst>highCost) = highCost; % also ^2?

              fcost(:,:,k) = dst;

              
            case {'centerline'}
              if self.useMex
                
                %% use the mex-file
                dst = fish.helper.pdist2CenterLine(cat(3,self.tracks.centerLine),self.centerLine);
                
                nidx = any(isnan(dst),2);
                %dst = bsxfun(@rdivide,dst,thresDist).^2; % better square?

                dst = (dst/thresDist);%.^2; % better square?
               %dst(dst>highCost) = highCost;
                
                
                % makes not sense to compare the distance to others if one is nan
                dst(nidx,:) = 1;
                
                fcost(:,:,k) = dst;
                
              else
                fcost(:,:,k) = 1;
              end

              
            case {'overlap'}

              bbsegs = double(cat(1,self.segments.BoundingBox));
              bbtracks = double(cat(1,self.tracks.bbox));
              overlap = getBBoxOverlap(bbtracks,bbsegs,0);
              fcost(:,:,k) = 1-overlap ;%+ (overlap==0)*somewhatCostly;
              
            case {'classprob','classification','classifier'}

              if self.isInitClassifier
                clProb = cat(1,self.tracks(:).classProb);
                dst = fish.helper.pdist2Euclidean(clProb,self.classProb);%correlation??!??
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
              fcost(:,:,k) = fish.helper.pdist2Euclidean(trackFeatures,detectedFeatures);
          end % switch

        end % costtype k loop
        
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
        costOfNonAssignment =  self.maxAssignmentCost + self.meanAssignmentCost*self.opts.tracks.costOfNonAssignment;

        % determine cost of each assigment to tracks
        fcost = self.computeCostMatrix();

        % scale the cost
        sfcost = sum(bsxfun(@times,fcost,shiftdim(scale(:),-2)),3);

        % Solve the assignment problem. First for the visible
        nvis = cat(1,self.tracks.consecutiveInvisibleCount);
        validIdx = find(~nvis);

        [self.assignments, self.unassignedTracks, self.unassignedDetections] = ...
            fish.helper.assignDetectionsToTracks(sfcost(validIdx,:),  costOfNonAssignment);

        self.unassignedTracks = validIdx(self.unassignedTracks);
        self.assignments(:,1) = validIdx(self.assignments(:,1)); 

        invisibleIdx = find(nvis);       

        if ~isempty(self.unassignedDetections) && ~isempty(invisibleIdx)
          % Solve the assignment problem for the invisible tracks
          sfcost1 = sfcost(invisibleIdx,self.unassignedDetections);

          costOfNonAssignment1 = costOfNonAssignment*(1+nvis(invisibleIdx)*self.opts.tracks.invisibleCostScale);

          if self.opts.tracks.adjustCostDuringCrossing
            crossMsk = ~self.testHandled();
            crossMsk = crossMsk(invisibleIdx);
            costOfNonAssignment1(crossMsk) = costOfNonAssignment1(crossMsk)*self.opts.tracks.crossingCostScale;
          end

          [assignments, unassignedTracks, unassignedDetections] = ...
              fish.helper.assignDetectionsToTracks(sfcost1,  costOfNonAssignment1',max(costOfNonAssignment1)*ones(1,size(sfcost1,2)));

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
          
          dt = (self.timeStamp-self.lastTimeStamp);
          idx = find(velocity>self.maxVelocity*dt*self.fishlength);

          if ~isempty(idx) && dt
            self.unassignedTracks = [self.unassignedTracks; trackIdx(idx)];
            self.unassignedDetections = [self.unassignedDetections; detectionIdx(idx)];
            self.assignments(idx,:) = [];
            fish.helper.verbose('%d assignment(s) exceeded velocity.\r',length(idx));
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
        mt = max(min(self.opts.tracks.costtau*self.avgTimeScale,self.currentFrame),1);
        self.meanCost = self.meanCost*(mt-1)/mt  + 1/mt*mean(reshape(fcost,[],size(fcost,3)),1);

        % save the mean assignment costs
        if ~isempty(self.assignments)
          costs = self.cost(self.assignments(:,2));
          if ~isempty(costs)
            self.meanAssignmentCost = self.meanAssignmentCost*(mt-1)/mt ...
                + 1/mt*mean(costs);
            self.maxAssignmentCost = self.maxAssignmentCost*(mt-1)/mt ...
                + 1/mt*max(costs);
          end
        end
        
        %Plotting
        if self.displayif && self.opts.display.assignments && ~isempty(self.assignments) 
          figure(19);
          plotWrongAssignments(self,1:length(self.assignments(:, 1)));
          drawnow;
        end
        
        if self.displayif && self.opts.display.assignments && self.currentFrame>4
          plotAssignmentDetails(self);
        end
        
        % compute the class prob switching matrix
        if length(self.unassignedDetections)>2
          fish.helper.verbose('there were %d unassigned detections\r',length(self.unassignedDetections))
        end

        
      end
      
    end
    
    
    
    function updateTracks(self)

    %% update assigned
      assignments = self.assignments;
      numAssignedTracks = size(assignments, 1);
      trackopts = self.opts.tracks;
      maxClass = 0;
      
      for i = 1:numAssignedTracks
        trackIdx = assignments(i, 1);
        detectionIdx = assignments(i, 2);
        centroid = self.centroids(detectionIdx, :);

        %% Correct the estimate of the object's location  using the new detection.
        if trackopts.kalmanFilterPredcition
          correct(self.tracks(trackIdx).predictor, [centroid]);
        end

        tau = trackopts.tauVelocity*self.avgTimeScale;
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
            
            % additional test whether the distance is also shorter
            cL = self.tracks(trackIdx).centerLine;
            flag = 0;
            if ~isempty(cL)
              cldist = (centerLine-cL).^2;
              cldistrev = (centerLine-cL(end:-1:1,:)).^2;
              if sum(cldist(:))<sum(cldistrev(:))
                flag = 1;
              end
            end
            
            if ~flag
              %reverse
              reversed = 1;
              centerLine = centerLine(end:-1:1,:);
              thickness = thickness(end:-1:1);
              
              if self.opts.detector.fixedSize %% only for saveing needed
                self.segments(detectionIdx).FilledImageFixedSizeRotated = ...
                    self.segments(detectionIdx).FilledImageFixedSizeRotated(:,end:-1:1,:);
              end
              
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


        if reversed  && self.opts.classifier.discardPotentialReversed
          % not update.. fishfeature reversed (cannot be changed outside of mex)
          probNoise = NaN;
          classprob = nan(1,self.nfish);
        else
          self.tracks(trackIdx).classProb =  classprob; % only update classprob in tracks if not reversed
        end

        if self.opts.classifier.onlyDAGMethod
          w = ~isnan(probNoise); % ignore noise
          reasonable = ~~w;
        else
          % always update history, though (weight will be zero for NaN)
          [reasonable w]  = self.tracks(trackIdx).classProbHistory.update(classprob, probNoise);
        end
        
        if w>0
          %update moving average
          tmp = (1-1/self.clpMovAvgTau)*self.tracks(trackIdx).clpMovAvg;
          self.tracks(trackIdx).clpMovAvg =  tmp  + (w/self.clpMovAvgTau) * classprob;

          maxClass = max(maxClass,max(classprob));
        end
        self.tracks(trackIdx).classProbW = w;
        
        
        
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
        thres =self.maxAssignmentCost + self.meanAssignmentCost;

        if self.cost(detectionIdx) <= thres && (reasonable  || ~self.isInitClassifier)
          batchIdx = self.tracks(trackIdx).nextBatchIdx;
          self.tracks(trackIdx).batchFeatures(batchIdx,:)  = self.idfeatures(detectionIdx,:);
          self.tracks(trackIdx).nextBatchIdx = min(batchIdx+1,self.maxFramesPerBatch);
        else
          self.tracks(trackIdx).nIdFeaturesLeftOut = self.tracks(trackIdx).nIdFeaturesLeftOut +1;
        end
        
        self.tracks(trackIdx).centroid = centroid;
        self.tracks(trackIdx).assignmentCost =  self.cost(detectionIdx);
        
        %% Update visibility.
        self.tracks(trackIdx).totalVisibleCount =  self.tracks(trackIdx).totalVisibleCount + 1;
        self.tracks(trackIdx).consecutiveInvisibleCount = 0;
      end

      %% update max classprob for reassignment scaling
      if maxClass
        mt = min(self.opts.tracks.costtau*self.avgTimeScale,self.currentFrame-self.nFramesForInit);
        self.maxClassificationProb = self.maxClassificationProb*(mt-1)/mt + 1/mt*maxClass;
      end
      
      
      %% update unassigned
      for i = 1:length(self.unassignedTracks)
        ind = self.unassignedTracks(i);
        self.tracks(ind).age = self.tracks(ind).age + 1;
        self.tracks(ind).consecutiveInvisibleCount =  self.tracks(ind).consecutiveInvisibleCount + 1;


        % do not update the tracks with nan, but use clp average instead
        self.tracks(ind).classProb = self.tracks(ind).clpMovAvg;
        
        if ~self.opts.classifier.onlyDAGMethod
          % only nan into history:
          self.tracks(ind).classProbHistory.update(nan(1,self.nfish),NaN);
        end
        
      end
      
      %% update DAG
      if self.isInitClassifier % classifier should be init
        self.daGraph.updateFromTracks(self.tracks,self.fishlength);

        
        %% predict fishIDs according to DAG
        currentPredFishIds = [self.tracks.predFishId];

        predFishIds = self.daGraph.predictFishIds();
        trackIndices = find(currentPredFishIds ~= predFishIds);
        if ~isempty(trackIndices)
          % switching!
          self.resetBatchIdx(trackIndices); 
          if self.verbosity>1
            fish.helper.verbose('Dag switched a fish...\r');
          end
        end
        for i = 1:length(trackIndices)
          self.tracks(trackIndices(i)).predFishId = predFishIds(trackIndices(i));
        end
        if self.opts.classifier.onlyDAGMethod
          for i = 1:length(self.tracks)
            self.tracks(i).fishId = predFishIds(i);
          end
        end

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
      lostInds =  lostInds | [self.tracks(:).consecutiveInvisibleCount] >= self.opts.tracks.invisibleForTooLong;

      if any(lostInds)

        % Delete lost tracks.
        lostTrackIds = [self.tracks(lostInds).id];
        fishIds = [self.tracks.fishId];
        if any(~isnan(fishIds(lostInds)))
          error('want to delete a true track');
        end
        
        self.tracks = self.tracks(~lostInds);

        % WE SHOULD PROB AVOID DELTEING TRACKS CURRENTLY CROSSING !
        %delete the track ids from the crossings
        for i = 1:length(self.tracks)
          self.tracks(i).crossedTrackIds = setdiff(self.tracks(i).crossedTrackIds,lostTrackIds);
        end
        
        fish.helper.verbose('Deleted %d tracks',sum(lostInds));
      end
      
    end
    
    
    function createNewTracks(self)
    % create new tracks for unassigndetections (if less the number of fish)

      if length(self.tracks)>=self.nfish
        return
      end
      
      availablefishids = setdiff(1:self.nfish,cat(2,self.tracks.fishId));
      if  ~self.opts.tracks.withTrackDeletion && ~length(availablefishids)
        return;
      end
      
      
      if self.isInitClassifier
        % already have fish.  take only one detection with least
        % assignment cost to the other classes
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
        bfeat = zeros(self.maxFramesPerBatch,size(idfeat,2),'single');
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
        
        if self.opts.classifier.onlyDAGMethod
          history = [];
          w = ~isnan(clprobnoise);
        else
          history = fish.core.FishClassProbHistory(self.nfish);
          % update the classprobhistory and set the parameters
          history.lambda = self.fishlength/self.fishwidth;
          [~,w] = history.update(clprob,clprobnoise);
        end
        
        % Create a new track.
        newTrack = struct(...
          'id',        self.nextId,       ...
          'fishId',   newfishid,...
          'predFishId',   newfishid,...
          'bbox',      self.bboxes(i,:),  ...
          'centroid',  self.centroids(i,:),  ...
          'velocity', [0,0], ...
          'centerLine',  centerLine,...
          'thickness',  thickness,...
          'segment',   self.segments(i),  ... 
          'features',   feat,  ... 
          'classProb',clprob,...
          'classProbW',w,...
          'classProbHistory',history,...
          'clpMovAvg',zeros(size(clprob)),...
          'batchFeatures',bfeat,          ...
          'nextBatchIdx', 1,              ...
          'predictor', pred,              ...
          'firstFrameOfCrossing', self.currentFrame, ...
          'lastFrameOfCrossing', self.currentFrame, ...
          'crossedTrackIds',[],...
          'age',       1,                 ...
          'totalVisibleCount', 1,         ...
          'consecutiveInvisibleCount', 0, ...
          'switchedFrames',0,...
          'assignmentCost',Inf, ...
          'nIdFeaturesLeftOut',0,...
          'stmInfo',[]);

        
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
      bbox = cat(1,self.tracks.bbox);
      centroids = cat(1,self.tracks.centroid);
      msk = fish.helper.pdist2Euclidean(centroids,centroids)<self.currentCrossBoxLengthScale*self.fishlength;

      bBoxXOverlap = bsxfun(@ge,bbox(:,1) + bbox(:,3), bbox(:,1)') & bsxfun(@lt,bbox(:,1), bbox(:,1)');
      bBoxYOverlap = bsxfun(@ge,bbox(:,2) + bbox(:,4), bbox(:,2)') & bsxfun(@lt,bbox(:,2), bbox(:,2)');
      bBoxOverlap = bBoxXOverlap & bBoxYOverlap;

      costThres = self.maxAssignmentCost + self.opts.classifier.crossCostThres*self.meanAssignmentCost;
      for i = 1:length(self.tracks)

        if sum(msk(i,:))>1  || sum(bBoxOverlap(i,:))>1  
          if self.tracks(i).consecutiveInvisibleCount>0 || self.tracks(i).assignmentCost>costThres 
            crossmat(i,:) =  msk(i,:)  | bBoxOverlap(i,:);
          end
        end
      end
      
      %make symmetric CAUTION: MIGHT NOT BE SYMMETRIC!
      if any(crossmat(:))
        crossmat = crossmat | crossmat';
        [~,~,crossings] = fish.helper.networkComponents(crossmat);
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

      if self.nfish>length(self.tracks)
        self.createNewTracks();
      end
      
      if self.opts.tracks.withTrackDeletion 
        self.deleteLostTracks(); 
      end

    end
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % Results & Plotting helper functions
    %--------------------
    
    function savedTracks = appendSavedTracks(self,savedTracks)
      
      
      if self.opts.tracks.withTrackDeletion
        fish.helper.verbose('WARNING: DO NOT CONVERT SAVEDTRACKS')
        return
      end
      
      s2mat = fish.helper.strucarr2strucmat(savedTracks);
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
    
    
    
    function trackinfo = getCurrentTracks(self)
    % gets the track info 
      
      if self.opts.tracks.keepFullTrackStruc && length(self.tracks)==self.nfish 
        % ignore initial times when not all tracks are created
        self.savedTracksFull = cat(1,self.savedTracksFull,self.tracks);
      end

      
      if isempty(self.saveFieldSegIf)
        self.saveFieldSegIf = zeros(length(self.saveFields),1);
        for i_f = 1:length(self.saveFields)
          f = self.saveFields{i_f};
          idx = find((f=='.'));
          if (length(idx)>1)
            error('Cannot save 2nd level features');
          end

          if ~isempty(idx)
            saveFieldSegIf(i_f) = 1;

            if ~strcmp(f(1:idx(1)-1),'segment') 
              error('expect "segment" if "." is found in saveField ');
            end

            self.saveFieldsIn(1).(f(idx(1)+1:end)) = [];

            fsub = f;
            fsub(idx) = '_';
            self.saveFieldsOut(1).(fsub) = [];
            
          else
            self.saveFieldsIn(1).(f) = [];
            self.saveFieldsOut(1).(f) = [];
          end
        end
      end
      
      % copy with mex
      trackinfo = getCurrentTracks_(self.nfish,self.tracks,self.saveFieldsIn,self.saveFieldsOut,self.saveFieldSegIf);

      % update time
      [trackinfo.t] = deal(self.timeStamp);
    end
    
    
    function checkVideoHandler(self)
    % checks wether the videoHandler is still OK
      
      
      if ~iscell(self.videoFile) && ~exist(self.videoFile)
        fish.helper.verbose('Cannot find videofile: %s',self.videoFile);
        self.videoFile = fish.helper.getVideoFile(self.videoFile);
        self.videoHandler = [];
        self.videoHandler = self.newVideoHandler(self.videoFile,self.timerange);
      end
      

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

      if ~isempty(savepath)
        savemat = [savepath filesep savename '.mat'];
      else
        savemat = [savename '.mat'];
      end
      exists = ~~exist(savemat,'file');
      
    end

  end
  


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  PUBLIC METHODS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Constructor
    %--------------------
    function self = Tracker(vid,varargin) 
    % FISH.TRACKER(VID,...) starts a FishTracker objects on a given video file
    % VID. If VID==[] a uidialog is opened to select the video file. 
    % 
    % Options can the set with FISH.TRACKER(VID,OPTNAME1,OPTVALUE1,...) or
    % FISH.TRACKER(VID,OPTS) where OPTS is a structure with fields identical
    % to the names of the options.
    %
    % A number of options are given in "avgBLC" which is average body
    % length covered and is a measure of the time it takes to cover
    % the distance with the avgVelocity.
    %  
    % To see information about the possible options, run (without arguments):
    %  >> fish.Tracker
    %
    % Example:
    %  >> ft = fish.Tracker([],'nfish',3);
    %  >> ft.track([0,20]); 
    %  >> ft.plot();



      def.vid = [];
      doc.vid = 'Video file name or {CAMIDX,''WRITEFILENAME''}';
      
      %% properties
      
      def.opts.nfish = [];
      doc.nfish = {'Number of fish. Needs to be fixed ' '(*attempt* to be estimated if empty)'};
      
      def.opts.fishlength = [];
      doc.fishlength = {'Approx length of fish in pixel (estimated if empty)'};
      def.opts.fishwidth = [];
      doc.fishwidth = {'Approx width of fish in pixel (estimated if empty)'};

      def.opts.maxVelocity = [];
      doc.maxVelocity = {'Maximal velocity in BL/sec (estimated if empty)'};

      def.opts.avgVelocity = 4;
      doc.maxVelocity = {'Approx. avg velocity. Important parameter ', ...
                          'to set the time scale of the tracking problem'};
      
      def.opts.useScaledFormat = false;
      doc.useScaledFormat = {'Use adaptive scaled gray format (EXPERIMENTAL)'};
      
      def.opts.stmif = false;
      doc.stmif = 'Use online visual stimulation';

      def.opts.useMex = true;
      doc.useMex = {'Use the C++ version FishVideoHandler (FAST)'};

      def.opts.useOpenCV = true;
      doc.useOpenCV = 'uses mexopencv library (FAST)';

      def.opts.useKNN = false;
      doc.useKNN = {'Use the KNN instead of thresholding ' 'for background segmentation (SLOW)'};

      def.opts.writefile = '';
      doc.writefile = 'For saving the tracking progress to video';

      def.opts.verbosity = 2;
      doc.verbosity = {['Sets verbosity level. 0:off, 1:moderate, >1: ' ...
                        'very verbose'],''};

      


      %% cost options

      % 0 turns off cost
      def.opts.cost.Location = 10;
      doc.cost.Location = {'Centroid comparison','Relative cost weighting. 0 turns cost type off.'};
      
      def.opts.cost.Overlap = 5;
      doc.cost.Overlap = {'BBox overlap'};

      def.opts.cost.CenterLine = 3;
      doc.cost.CenterLine = {'CenterLine feature comparison'};

      def.opts.cost.Classifier = 2;
      doc.cost.Classifier = {'Classification prob comparison'};

      def.opts.cost.Size = 1;
      doc.cost.Size = {'Major/Minor axis comparison'};
      
      def.opts.cost.Area = 1;
      doc.cost.Area = {'blob pixel area comparison'};

      def.opts.cost.BoundingBox = 0; 
      doc.cost.BoundingBox = {'[x,y,w,h] bbox comparison',''};

      
      %% detector options
      def.opts.detector(1).history = 500;  %250 [nframes]
      doc.detector(1).history = 'Background update time constant [nFrames]';
      
      def.opts.detector.inverted = false;  
      doc.detector.inverted = {'Set 1 for IR videos (white fish on dark background)'};
      
      def.opts.detector.adjustThresScale = 0.95;   
      doc.detector.adjustThresScale = {'0.8..1 : reduce when detections too noisy (useKNN=0)',''};
      

      %% reader
      def.opts.reader(1).resizeif = false;
      doc.reader(1).resizeif = {'Whether to resize the frame before tracking'};

      def.opts.reader.resizescale = 0.75; 
      doc.reader.resizescale = {'Fraction to resizing the frame for RESIZEIF = true',''};
      
      %% blob anaylser

      def.opts.blob(1).colorfeature = false; 
      doc.blob.colorfeature = {'Use color information for fish feature'};

      def.opts.blob.headprop = 0.6; 
      doc.blob.headprop = {'Proportion of object length to use as feature',''};

      
      %% classification
      def.opts.classifier.npca = 40; 
      doc.classifier.npca = 'Number of PCA components';

      def.opts.classifier.nlfd = -1; 
      doc.classifier.nlfd = {'Number of LFD components. Set to 0 to turn off.',...
                          '-1 for auto setting.'};

      def.opts.classifier.tau = 5000; 
      doc.classifier.tau = {'Slow time constant of classifier [nFrames].'};

      def.opts.classifier.reassignProbThres = 0.2; %0.45
      doc.classifier.reassignProbThres = {'minimal probability for reassignments'};

      def.opts.classifier.forcedReassignProbThres = 0.6; %0.45
      doc.classifier.reassignProbThres = {'minimal probability for reassignments'};

      
      def.opts.classifier.handledProbThres = 0.2; %0.45
      doc.classifier.handledProbThres = {'minimal diff probability for crossing exits'};


      def.opts.classifier.crossCostThres = 3; 
      doc.classifier.crossCostThres = {'candidates for crossings: scales mean assignment cost',''};


      %% tracks
      def.opts(1).tracks.costtau = 100;
      doc.tracks.costtau = {'Time constant for computing the mean cost [avgBLC]'};

      def.opts.tracks.crossBoxLengthScale = 1; 
      doc.tracks.crossBoxLengthScale = {'How many times the bbox is regarded as a crossing'};
      
      def.opts.tracks.adjustCostDuringCrossing = true; 
      doc.tracks.adjustCostDuringCrossing = {'Whether to scale non-assignment cost during crossings'};

      def.opts.tracks.keepFullTrackStruc = false;
      doc.tracks.keepFullTrackStruc = {'Whether to keep full track structure. ONLY FOR DEBUG!'};
      
      def.opts.tracks.costOfNonAssignment = 4;
      doc.tracks.costOfNonAssignment =  {'Scales the threshold for','cost of non assignment'};

      def.opts.tracks.invisibleCostScale = 1;%1 
      doc.tracks.invisibleCostScale = {'Factor for nonAssignmentCost per frame with','invisible count'};


      def.opts.tracks.crossingCostScale =  1;
      doc.tracks.crossingCostScale = {'Scaling of non-assignment cost during crossings'};

      def.opts.tracks.probThresForFish = 0.1; 
      doc.tracks.probThresForFish = {'Classification probability to ' 'assume a fish feature'};

      def.opts.tracks.useDagResults = 1;
      doc.tracks.useDagResults = {'Sets default output results to ' 'DAG (1) or Switch (0) method',''};
      
      def.opts.tracks.kalmanFilterPredcition = false; 
      doc.tracks.kalmanFilterPredcition = {'Whether to use Kalman filter'};
      
      def.opts.classifier.timeToInit = 20;  
      doc.classifier.timeToInit = {'Time to initialize the classifier [avgBLC]'};

      def.opts.classifier.timeAfterCrossing =  2; 
      doc.classifier.timeAfterCrossing = {['When to check for permutations after ' ...
                          'crossings [avgBLC]']};
      
      def.opts.classifier.timeForUniqueUpdate = 12; 
      doc.classifier.nFramesForUniqueUpdate = {['Unique frames needed for update all ' ...
                          'fish simultaneously [avgBLC]']};

      def.opts.classifier.clpMovAvgTau = 1; 
      doc.classifier.clpMovAvgTau = {'Time constant of class prob','moving average [avgBLC].'};

      def.opts.tracks.tauVelocity = 1; 
      doc.tracks.tauVelocity = {'Time constant to compute the ' 'velocity [avgBLC]'};


      
      %% dag
      def.opts.dag.probScale = 0.5;
      doc.dag.probScale = {'DAGraph probScale. 1 means 50/50 weighting of ', ...
                          'classprob with distance if points a fishLength apart',''};
      
      %% stimulus
      def.opts.stimulus.presenter = [];
      doc.stimulus.presenter = {'Stimulus presenter object  (for STMIF=true)'};
      
      def.opts.stimulus.screen = 1;
      doc.stimulus.screen = 'Screen number to use.';

      def.opts.stimulus.screenBoundingBox = [];
      doc.stimulus.screenBoundingBox = {['Bbox of stimulus screen in frame ' ...
                          'pixels'],'use calibrateScreen to get estimate.'};

      
      
      %% display opts
      def.opts.displayif = 3;
      doc.displayif = {'Turn on/off all displaying'};

      def.opts.display.displayEveryNFrame = 20;
      doc.display.displayEveryNFrame = {'How often to update the track display window (SLOW)'};

      def.opts.display.tracks = true;
      doc.display.tracks = {'Whether to plot the tracking process'};

      def.opts.display.level = 3;
      doc.display.level = {'Level of details of the track plot'};
      
      def.opts.display.fishSearchResults = false;
      doc.display.fishSearchResults = {'Info plot nfish auto-search'};

      def.opts.display.stimulusProgress = true;
      doc.display.stimulusProgress = {'ProgressBar in case of stimulation'};
      
      def.opts.display.switchFish = false;
      doc.display.switchFish = {'Switch fish info plot (for DEBUGGING) '};

      def.opts.display.videoHandler = false;
      doc.display.videoHandler = {'Raw frames and bwmsk MEX only (for DEBUGGING) '};
      

      def.opts.display.assignments = false;
      doc.display.assignments = {'Assignment info plot (for DEBUGGING) '};

      def.opts.display.calibration = true;
      doc.display.calibration = {'Stimulus calibration results',''};
      
      
      STRICT = 1; % be not too strict. Some settings of the videoHandler are not
                  % explicitly given here
      VERBOSELEVEL = 0;
      ALLFIELDNAMES = 1;
      
      MFILENAME = 'fish.Tracker';
      fish.helper.parseInputs;
      if HELP;self.videoHandler = [];return;end

      global VERBOSELEVEL;
      VERBOSELEVEL = self.verbosity;

      
      %%options depending on avgVelocity (will be set in init)

      opts.display.leakyAvgFrameTau = 50;
      doc.display.leakyAvgFrameTau = {'Tau for switch fish plot'};


      %%options not really important or debug
      opts.detector.fixedSize = 0;  
      doc.detector.fixedSize = {'Set 1 for saving the fixedSize image'};

      opts.detector.nskip = max(ceil(opts.detector.history/50),1); 
      doc.detector.nskip = 'Skip frames for background (useKNN=0)';

      opts.display.detectedObjects = false;
      doc.display.detectedObjects = {'Detection info plot (for DEBUGGING) '};
      
      opts.display.crossings = false;
      doc.display.crossings = {'Crossings info plot (for DEBUGGING) '};

      opts.display.classifier = false;
      doc.display.classifier = {'Classifier plot (for DEBUGGING) '};

      opts.display.assignmentCost = false;
      doc.display.assignmentCost = {'Assignment cost info plot (for DEBUGGING) '};

      opts.blob.computeMSERthres= 2; % just for init str
      doc.blob(1).computeMSERthres =  {'When to compute MSER (SLOW; only ' 'for useMex=0)'};

      opts.blob.computeMSER= 0;
      doc.blob.computeMSER = {'Whether to compute MSER'};

      opts.blob.difffeature = true; 
      doc.blob.difffeature = {'Use background subtracted gray' 'images for fish feature'};

      opts.blob.interpif = 1;
      doc.blob.interpif = {'Whether to interpolate while','debending'};

      opts.tracks.crossBoxLengthScalePreInit = max(opts.tracks.crossBoxLengthScale*0.75,0.5); 
      doc.tracks.crossBoxLengthScalePreInit = {'Before classifier init ' '(reduced somewhat)'};

      opts.classifier.discardPotentialReversed = opts.blob.headprop~=1;
      doc.classifier.discardPotentialReversed = {'Discard potential reversed during','classification'};
      

      %% other developmental parameters

      %lost tracks
      opts.tracks.invisibleForTooLong = 15; % [nFrames] on track basis. Only for deletion
      opts.tracks.ageThreshold = 10; % [nFrames]
      opts.tracks.withTrackDeletion = false; % BUG !!! TURN OFF. maybe needed later 
      
      opts.classifier.onlyDAGMethod = 0;
      if opts.classifier.onlyDAGMethod
        fish.helper.verbose('Switch method: DAG');
      end
      
      %% parameter checking

      if isempty(vid)
        vid = fish.helper.getVideoFile();
      end

      if iscell(vid)
        if length(vid)==1
          vid{2} = '';
        end
        vname = vid{2};
      else
        vname = vid;
      end
      
      if ~isempty(vname)
        if all(vname(1:2)=='~/') && isunix()
          [~,home] = unix('eval echo ~$USER');
          vname = [home(1:end-1) vname(2:end)];
        elseif vname(1)=='~'
          error('provide full path name. Cannot start with "~"');        
        end

        if ~exist(vname) && ~iscell(vid)
          error(sprintf('Video file "%s" not found',vname));
        end
      end
      
      if ~isempty(opts.nfish) && opts.nfish<0 && ~iscell(vid)
        opts.nfish = fish.helper.chooseNFish(vname,1);
      end
      
      if ~iscell(vid)
        vid = vname;
      else
        vid{2} = vname;
      end
      
      % construct object
      self.opts = opts;
      self.setProps();
      self.setupSystemObjects(vid); % will call setOpts
      
    end
    
    
    function disp(self)
      if isempty(self.videoHandler)
        fprintf('Empty fish.Tracker Object\n\n');
      else
        disp@handle(self);
      end
    end
    
    
    function trackStep(self)
    % Same commands as one step in tracking loop of SELF.TRACK(). However,
    % can be called from outside (DEPRECIATED) to continue tracking
    % in a stepwise manner for plotting reasons. However, NO ERROR
    % CHECKING is provided!
      
      self.currentFrame = self.currentFrame + 1;
      self.segments = self.videoHandler.step();
      
      self.detectObjects();
      self.handleTracks();
    end

    
    
    function track(self,trange,saveif,writefile)
    % TRACK([TRANGE],[SAVEIF],[WRITEFILE]). Detect objects, and track them
    % across video frames. TRANGE=[TSTART,TEND] specifies the time
    % range of tracking (in seconds). Set TRANGE=[] for tracking the
    % while video file.  SAVEIF==1 turns on automatic results
    % saving. If WRITEFILE is set, the video-output (if DISPLAYIF>0)
    % will be saved to the given video file.
    %  
    % CAUTION:
    % Previous tracking data will be overwritten.
    %
    % EXAMPLE:
    %  >> ft = fish.Tracker(videoFile,'nfish',3,'displayif',2);
    %  >> ft.track([0,100]);
    %  >> ft.plot();
    %  >> ft.save();
      
      if exist('writefile','var') && ~isempty(writefile)
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
      
      %% initializing
      self.initTracking();

      if ~self.stmif
        fish.helper.verbose('Start tracking of %d fish  in the time range [%1.0f %1.0f]...',...
                            self.nfish,self.timerange(1),self.timerange(2));
      else
        fish.helper.verbose('Start tracking of %d fish with stimulation', self.nfish);
      end
      
      %% tracking loop
      localTimeReference = tic;
      localTime = toc(localTimeReference);
      verboseTime = toc(localTimeReference);
      s = 0;
      while  hasFrame(self.videoHandler) && ...
          (~self.opts.display.tracks || ~self.displayif || isOpen(self.videoPlayer)) ...
          && (~self.stmif || ~self.stimulusPresenter.isFinished(localTime))

        
        %% main tracking step
        s = s+1;
        %self.trackStep(); % avoid calling this..just paste for performance.
        %unnecessary function call
        self.currentFrame = self.currentFrame + 1;
        self.lastTimeStamp = self.timeStamp;
        if self.displayif && self.opts.display.switchFish
          [self.segments, self.timeStamp, frame] = self.videoHandler.step();
          self.computeLeakyAvgFrame(frame);
        else
          [self.segments,self.timeStamp] = self.videoHandler.step(); % faster..
        end
        
        self.detectObjects();
        self.handleTracks();

        
        %% stimulus presenting
        localTime = toc(localTimeReference);        
        if self.stmif && ~isempty(self.stimulusPresenter) 
          self.tracks = self.stimulusPresenter.step(self.tracks,self.videoHandler.frameSize,localTime);
        end
        
        %% saving 
        if self.opts.tracks.withTrackDeletion
          savedTracks(1:length(self.tracks),s) = self.getCurrentTracks();
        else
          savedTracks(1:self.nfish,s) = self.getCurrentTracks();
        end

        if ~mod(s,self.nFramesAppendSavedTracks) 
          savedTracks = self.appendSavedTracks(savedTracks);
          s = size(savedTracks,2);
          
          if ~self.stmif &&  ~mod(self.currentFrame,self.nStepsSaveProgress*self.nFramesAppendSavedTracks) && saveif
            % occasionally save for long videos 
            fish.helper.verbose('Save tracking progress to disk..');
            self.save();
          end
        end
        
        %% displaying
        if self.displayif && self.opts.display.tracks
          if ~mod(self.currentFrame-1,self.opts.display.displayEveryNFrame)
            self.displayCurrentTracks();
          end
        end
        
        %% verbose
        if ~mod(self.currentFrame,self.nFramesVerbose) && ~isempty(self.tracks)
          tt = toc(localTimeReference);
          verboseDuration = tt - verboseTime;
          verboseTime = tt;        
          trackDuration = self.tabs(self.currentFrame)...
              - self.tabs(max(self.currentFrame-self.nFramesVerbose,1));
          t = datevec(seconds(self.timeStamp));
          if ~self.stmif
            fish.helper.verbose(['%1.0fh %1.0fm %1.1fs  [%1.2f x real ' ...
                                'time]             \r'],t(4),t(5),t(6),trackDuration/verboseDuration);
          else
            fish.helper.verbose(['%1.0fh %1.0fm %1.1fs  [%1.2f FPS] ' ...
                                '                  \r'],t(4),t(5),t(6),self.nFramesVerbose/verboseDuration) ;  
          end
        end
      end

      %% finishing
      self.cleanup();

      % make correct output structure
      if exist('savedTracks','var')
        self.appendSavedTracks(savedTracks);
        generateResults(self);

      end
      self.closeVideoObjects();
      
      if saveif
        % finally save results
        fish.helper.verbose('Save tracking results to disk...');
        self.save();
      end

    end
    
    
    %% 
    function st = getSavedTracksFull(self)
    % returns the savedTracksFull structure
      if ~self.opts.tracks.keepFullTrackStruc
        error(['To access the full structure, track with ' ...
               'opts.tracks.keepFullTrackStruc  turned on']);
      else
        st = self.savedTracksFull;
      end
    end
    
  end
  
end % class




