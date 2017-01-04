classdef Tracker < handle;


  properties 
    timerange= [];    
    maxVelocity = [];        
    avgVelocity = [];
    bodylength = [];
    bodywidth = [];
    nindiv = [];
    opts = struct([]);
  end

  properties (SetAccess = private)


    videoFile = '';
    writefile = '';
    currentFrame = 0;

    stmif = 0;
    useMex = 1;
    useOpenCV = 1;
    useScaledFormat = 0;
    useKNN = 0;
    hasVision = 0;
    
    verbosity = 1;
    displayif = 3;    

    
    videoHandler = [];
    videoPlayer = [];
    videoWriter = [];
    stimulusPresenter = [];
    
    identityClassifier = [];
    daGraph = [];
    tracks = [];
    res = [];

    meanAssignmentCost =1;
    maxAssignmentCost =1;
    maxClassificationProb =0;
    meanClassificationProb =0;
    meanClassificationProbOther =0;
    meanNonCrossingFrames = 0;

    avgTimeScale = [];

    nFramesForInit  = [];
    nFramesAfterCrossing = [];
    nFramesForUniqueUpdate = [];
    clpMovAvgTau = [];
    nFramesForSingleUpdate = [];

  end
  
  properties (SetAccess = private, GetAccess=private);

    saveFieldSegIf = [];
    saveFieldsIn = struct([]);
    saveFieldsOut = struct([]);

    costinfo= {}; 
    scalecost = [];

    saveFields = {};
    
    nextId = 1; 
    isInitClassifier = [];

    centerLine = [];
    thickness = [];
    bboxes = [];
    
    centroids = [];
    locations = [];
    idfeatures = [];
    features = [];
    cost = [];
    segments = [];
    
    classProb = [];
    classProbNoise = [];
    crossings = [];

    identityId2TrackId = [];
    meanCost = ones(4,1);
    pos = [];
    tabs = [];
    uniqueIdentityFrames = 0;
    nonCrossingFrames = 0;
    
    assignments = [];
    unassignedTracks = [];
    unassignedDetections = [];
    featurecosttypes = {};

    leakyAvgFrame = [];
    savedTracks = [];
    savedTracksFull = [];

    currentCrossBoxLengthScale = 1; 
    
    % constants
    nFramesExtendMemory = 100000;
    nFramesAppendSavedTracks = 100;
    nStepsSaveProgress = 250; % times nFramesAppendSavedTracks
    nFramesVerbose = 40;
    
    % autoinit
    maxFramesPerBatch = [];
    minBatchN = [];
    
    
    timeStamp = [];
    lastTimeStamp = [];
    dt = [];
    
    
    %    MAXCROSSBODY = 10; % MAX number of cross for many body....DEVELOPMENTAL
  end

  
  methods(Static)
    
    function success=checkCompiled()
      success = ~~exist('getCurrentTracks_');
    end

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

        opts.stmif = 0; % set to zero;
        opts.tracks.initBackground = 0;

        if iscell(S.videoFile)
          xy.helper.verbose('Realtime-stimulus experiment. Load grabbed video file.')
          S.videoFile = S.videoFile{2};
        end

        if ~exist(S.videoFile,'file') 
          [~,b,c] = fileparts(S.videoFile);
          if exist([b,c],'file') 
            S.videoFile = [b c];
          else
            xy.helper.verbose(['WARNING: could not find video file "%s". Please ' 'select!!'],[b c]);
            S.videoFile = '';
          end
        end

        obj = xy.Tracker(S.videoFile,'STRICT',0,opts);
        for f = fieldnames(S)'
          if isprop(obj,f{1}) 
            if ~any(strcmp(f{1},{'videoHandler','videoWriter','videoPlayer','videoFile'}))
              obj.(f{1}) = S.(f{1});
            end
          else
            if ~any(strcmp(f{1},deleted_fields))
              xy.helper.verbose('Cannot initialize field %s',f{1});
            end
          end
          
        end
      else
        obj = S;
      end
      if isa(obj.stimulusPresenter,'xy.stimulus.Presenter');
        obj.stimulusPresenter.closeProgressBar();
      end
    end
  
    
    %% parameter variations
          
    function [d,ddag, fL, t_elapsed] = parameterfun(opts1,parname,parvalue,tmax,dfname)

      if nargin>4 && ~isempty(dfname)
        global VERBOSEDIARY; %#ok<TLEV>
        VERBOSEDIARY = dfname;
      end
        
      
      maxNumCompThreads(4);
      opts = opts1;
      
      opts = xy.helper.setfield(opts,parname,parvalue);
      opts.nindiv = 5;
      opts.verbosity = 3;
      
      xy.helper.verbose('Start parameter variation: %s = %f',parname,parvalue);
     
      [~,t_elapsed,~,idpos, T] = xy.Tracker.runTest(tmax,opts,[],[],0);
      
      r = T.getDagTrackingResults();
      xyposdag = r.pos;
      
      ddag = squeeze(sqrt(sum((xyposdag - idpos).^2,2)));
      r = T.getSwitchBasedTrackingResults();
      xypos = r.pos;
      d = squeeze(sqrt(xy.helper.nansum((xypos - idpos).^2,2)));
      fL = T.bodylength;
      
      if nargin>4 && ~isempty(dfname)
        VERBOSEDIARY = 0;
      end
    end

    function res = parameterVariations(pv,xyopts,varargin)
    % RES=PARAMETERVARIATIONS(PV,XYTOPTS,...) computes parameter variations for
    % given parameters and checks how the runTest performs. PV is a struc of
    % field names like parmeters in xy.Tracker. The field contains the array
    % to loop over. XYTOPTS is an options structure for other paramters to be
    % fixed.
    % 
    % EXAMPLE: 
    %   pv.avgVelocity  = 1:10;
    %   pv.classifier.reassignProbThres = 0:0.1:0.5;
    %   res = parameterVariations(pv);  % submit tasks
    %   res = parameterVariations(res); % fetch and plot

      
      def.pv = [];
      def.xyopts = [];
      def.opts.diaryif = 0;
      def.opts.diarypath = 'diary';
      def.opts.tmax = [];

      VERBOSELEVEL = 0;
      MFILENAME = 'xy.Tracker.parameterVariations';      
      xy.helper.parseInputs;
      if HELP, return;end;


      if isfield(pv,'MFILENAME')
        SUBMIT = 0;
        res = pv;
       else
        res = INPUTS;
        SUBMIT = 1;
        res.translate = {};
        res.dat = [];
      end
      
      % re-format
      pars= xy.helper.allfieldnames(res.pv)';
      pararr = {};
      for i = 1:length(pars)
        pararr{i} = xy.helper.getfield(res.pv,pars{i});
      end


      if res.opts.diaryif
        basepath = [res.opts.diarypath filesep];
      end

      if SUBMIT
        s = 0;
        for i = 1:length(pars)
          for j = 1:length(pararr{i})
            
            s= s+1;
            
            if res.opts.diaryif
              fn = sprintf('%sdiary.%d.log',basepath,s);
            else
              fn = [];
            end
            
            
            res.f(s) = parfeval(@xy.Tracker.parameterfun,4,res.xyopts,pars{i},pararr{i}(j),res.opts.tmax,fn);
            res.translate{s} = [i,j];
          end
        end

        xy.helper.verbose('Submitted %d tasks.',s);
        xy.helper.verbose('Catch outputs by recalling the function with RES as argument...');
        return
      end

      if isempty(res.dat) || numel(cat(1,res.dat.bodyLength))<numel(res.translate)
        % get output
        xy.helper.verbose('Waiting for results..');
        for k = 1:length(res.f)
          
          ONLINE = 0;%isempty(res.dat);
          if ONLINE
              [completedIdx, d, ddag, bodyLength,t_elapsed] = fetchNext(res.f);
          else
            completedIdx = k;
            try
              [d, ddag , bodyLength,t_elapsed] = fetchOutputs(res.f(k));
            catch
              disp(res.f(k))
              continue;
            end
          end
          
          idx = res.translate{completedIdx};
          i = idx(1);
          j = idx(2);
          res.dat(i,j).d = d;
          res.dat(i,j).ddag = ddag;
          res.dat(i,j).t_elapsed = t_elapsed;
          res.dat(i,j).bodyLength = bodyLength;
          res.dat(i,j).parval = pararr{i}(j);
          res.dat(i,j).parname = pars{i};
          if opts.diaryif
            res.dat(i,j).log = sprintf('%sdiary.%d.log',basepath,completedIdx);
            %read log with: textread(res(i,j).log,'%s','delimiter','\n');
          end
          xy.helper.verbose('Fetched %d/%d results.\r',k,length(res.f));
        end
      
        if isempty(res.dat)
          error('Cannot find any results');
        end
      end
      

      
      %% plotting
      figure;

      init = 100;
      
      [r1,r2] = xy.helper.getsubplotnumber(size(res,1));
      %fun = @(x)max(x,[],2);
      %fun = @(x)xy.helper.nanmean(x,2);
      fun = @(x)xy.helper.nansum(x>0.5,2); % across body;
      
      dat = res.dat;
      nframes = size(dat(1,1).d,1);
      nindiv = size(dat(1,1).d,2);
      n = nframes*nindiv;
      
      for i = 1:size(dat,1)
        fL = cat(1,dat(i,:).bodyLength); 
        parr = cat(2,dat(i,:).parval);
        x  = bsxfun(@rdivide,cat(3,dat(i,:).d),permute(fL,[3,2,1]));
        x = x(init:end,:,:);
        
        xdag  = bsxfun(@rdivide,cat(3,dat(i,:).ddag),permute(fL,[3,2,1]));
        xdag = xdag(init:end,:,:);

        
        md = shiftdim(xy.helper.nanmean(fun(x),1));
        sd = shiftdim(xy.helper.stderr(fun(x),1));
        
        mdag = shiftdim(xy.helper.nanmean(fun(xdag),1));
        sdag = shiftdim(xy.helper.stderr(fun(xdag),1));

        nnan = shiftdim(sum(sum(isnan(cat(3,dat(i,:).d)),1),2))/n;
        nnandag = shiftdim(sum(sum(isnan(cat(3,dat(i,:).ddag)),1),2))/n;
        
        a = subplot(r1,r2,i);
        nc = 1;
        xy.helper.errorbarpatch(parr,conv2([md,mdag],ones(nc,1)/nc,'same'),[sd,sdag]);
        ylabel('% frames > 1 BL');
        legend('SWB','DAG');
        b = xy.helper.submargin(a,'MARGIN',0.01,'SPACE',0.2);
        hold on;
        plot(b,parr,[nnan,nnandag])
        xlabel(unique(cat(1,dat(i,:).parname),'rows'));    
        
        linkaxes([a,b],'x');
      end
    
    

    
    end
    
      

    %% test
    function [self] = runSimpleTest(tmax,plotif)
        if nargin<1
            tmax = 100;
        end
        if nargin<2
            plotif = 5;
        end
      [~,~,~,~,self] = xy.Tracker.runTest(tmax,[],[],[],plotif);
    end
    
    function [success,t_elapsed,varargout] = runTest(tmax,opts,pathToVideo,pathToMat,plotif)
    %[SUCCESS,T_ELAPSED] = RUNTEST(TMAX,OPTS,PATHTOVIDEO,PATHTOMAT,PLOTIF) makes a validation versus
    %the idTracker [1] with the video provided by idTracker. OPTS is a struct
    %with fields given to initialize xy.Tracker. PATHTOMAT is the trajectory
    %file from idTracker tracking the same video. PATHTOVIDEO is the video file
    %from idTracker (5 Zebrafish). TMAX sets the trange of the tracking
    %comparison (in seconds). Set TMAX to [] to test the whole
    %7.4min video. T_ELAPSED is the required tracking time in
    %seconds. if PLOTIF==0 no output is given. PLOTIF==1 plots a result
    %plot and PLOTIF>1 displays the tracks currently tracked.
    %
    % [..,XYTPOS,IDPOS] = RUNTEST returns additionally the position results.
    % [..,XYT] = RUNTEST returns additionally the xy.Tracker object.
    %  
    % Example: 
    % >> xy.Tracker.runTest(50,[],[],[],2); % with plotting
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
          pathToVideo = [xy.helper.getXYTRoot() 'data' filesep '5Zebrafish_nocover_22min.avi'];   
        end
        if ~exist(pathToVideo,'file')
          error(['Please set PATHTOVIDEO to dowloaded video file available at "http://' ...
                 'www.cajal.csic.es/files/gpolavieja/5Zebrafish_nocover_22min.avi" ' ...
                 'or put the downloaded avi-file into the "+xy/../data" directory.']);
        end
      end
      
      if nargin<4 || isempty(pathToMat)
        pathToMat = [xy.helper.getXYTRoot() 'data' filesep 'trajectories.mat']; 
      end

      if nargin<5
        plotif = 1;
      end
      
      MAXDISTANCE = 0;
      
      % handle ID
      id = load(pathToMat);
      idres.pos = permute(id.trajectories,[1,3,2]);

      % run benchmark
      xy.helper.verbose('Starting run test.');
      T = xy.Tracker(pathToVideo,'avgVelocity',5,'blob.colorfeature',false,args{:});
      T.setDisplay(max(plotif-1,0));
      T.addSaveFields('firstFrameOfCrossing', 'lastFrameOfCrossing');

      tic;
      T.track(trange);
      t_elapsed=toc;
      xy.helper.verbose('Ellapsed time for tracking video: %1.1fs',t_elapsed);
      
      nindiv = T.nindiv;
      xyres = T.getTrackingResults();
      xyres.pos(end,:,:) = NaN; % why ? 
      xyresnan = xyres;
      xyresnan.pos = T.deleteInvisible(xyres,'pos');
      
      dist = zeros(nindiv);
      offs = size(idres.pos,1) - size(xyres.pos,1);
      
      
      for i = 1:nindiv
        for j = 1:nindiv
          dist(i,j) = xy.helper.nanmean(sqrt(sum((xyres.pos(:,:,j) - idres.pos(1:end-offs,:,i)).^2,2)),1);
        end
      end
      assignments = xy.helper.assignDetectionsToTracks(dist,1e3);
      
      xypos = xyres.pos;
      xyposnan = xyresnan.pos;
      idpos = idres.pos(1:end-offs,:,assignments(assignments(:,1),2));
 
      t = xyres.t(:,1);
      d = sqrt(sum((xypos(:,:,:) - idpos(:,:,:)).^2,2))/T.bodylength;
      
      success = mean(max(d,[],3)>1)<0.05;

      varargout{1} = xypos;
      varargout{2} = idpos;
      varargout{3} = T;
      
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
          plot(t,max(d,[],3));
        else
          xy.helper.errorbarpatch(t,xy.helper.nanmean(d,3),xy.helper.stderr(d,3));
        end
        
        hold on;
        plot(t([1,end]),[1,1],'r:');
        
        set(a(s),'fontsize',8);
        if MAXDISTANCE
          ylabel(sprintf('Max distance\n[Body length]'),'fontsize',10);
        else
          ylabel(sprintf('Avg. distance\n[Body length]'),'fontsize',10);
        end
        
        xlim(t([1,end]));
        box off;
        
        %detections
        s = s+1;
        a(end+1) = subplot(r1,r2,s,'align');
        nconv = 75;
        Tnan = conv(sum(isnan(xyposnan(:,1,:)),3),ones(nconv,1)/nconv,'same');
        idnan = conv(sum(isnan(idpos(:,1,:)),3),ones(nconv,1)/nconv,'same');
        plot(t,[Tnan,idnan]*100/T.nindiv);
        xlim(t([1,end]))
        set(a(s),'fontsize',8);
        
        ylabel(sprintf('Avg. lost\n tracks [%%]'),'fontsize',10);
        
        legend({'xy.Tracker','idTracker'},'location','NorthWest','fontsize',8)
        box off;
        
        
        %crossings
        s = s+1;
        a(end+1) = subplot(r1,r2,s,'align');
        
        
        cross = diff(xyres.tracks.lastFrameOfCrossing-xyres.tracks.firstFrameOfCrossing)>0;
        ccross = conv2(sum(double(cross),2),ones(nconv,1)/nconv,'same');
        plot(t(1:end-1),ccross,'k')
        set(a(s),'fontsize',8);
        ylabel(sprintf('Avg. # of \ncrossing events'),'fontsize',10);
        xlim(t([1,end]));
        box off;
        
        
        % probability
        s = s+1;
        a(end+1) = subplot(r1,r2,s,'align');
        clp = xyres.tracks.classProb;
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
        
        h1 = xy.helper.errorbarpatch(t,mean(dclp,2),xy.helper.stderr(mclp,2),'color',[0,0.5,0]);
        hold on;
        h2 = xy.helper.errorbarpatch(t,mean(mclp,2),xy.helper.stderr(dclp,2),'color',[1.0000, 0.6445, 0]);
        xlim(t([1,end]));
        box off;
        ylabel(sprintf('Class\nprobability'),'fontsize',10);
        legend([h1(1),h2(1)],{'Avg.','Other'},'location','NorthWest','fontsize',8)
        set(a(end),'fontsize',8)
        
        
        xlabel('Time [sec]','fontsize',10);
        b = xy.helper.labelsubplot(gcf);
        xy.helper.shiftaxes(b,0.02)
        drawnow;
      end

      if nargout>1
        if ~success
          xy.helper.verbose('WARNING: RunTest failed due to high inaccuracies in the tracking.');
        end
      else
        assert(success,'RunTest failed due to high inaccuracies in the tracking.');
      end
    end

    
    function [screenBoundingBox] = calibrateStimulusScreen(camIdx,screenIdx,plotif)
    % [SCREENBOUNDINGBOX] =
    % CALIBRATESTIMULUSSCREEN(CAMIDX,SCREENIDX gets the
    % SCREENBOUNDINGBOX used for stimulus presentation.
    
      opts = [];
      opts.nindiv = 4;
      opts.bodylength = 100;
      opts.bodywidth = 50;
      opts.stmif = 1;
      opts.stimulus.presenter = 'xy.stimulus.PresenterCalibration';
      opts.stimulus.screen =   screenIdx;
      opts.classifier.timeToInit = Inf; 
      opts.detector.inverted = 1;
      opts.tracks.useDagResults = 0;
      opts.tracks.initBackground = 0;
      
      T = xy.Tracker({camIdx,''},opts);

      if nargin<3
        plotif = 1;
      end
      
      if plotif>1
        T.setDisplay(1);
        T.setDisplay('tracks',1,'stimulusProgress',0);      
      else
        T.setDisplay(0);
      end
      
      T.addSaveFields('bbox');

      T.stimulusPresenter.width = 50; % maybe needs to be adjusted
      T.stimulusPresenter.tmax = 30;
      T.stimulusPresenter.freq = 0.2;

      
      xy.helper.verbose(['\n\n****** \tStarting calibration. ' ...
                          'Make sure that the IR filter is NOT ' ...
                          'installed! \n\tHit Enter to proceed!\n']);
      pause;
      
      % start detecting
      T.track(); % track the markings
      T.stimulusPresenter.flip(); % turn stim off
      
      res = T.getTrackingResults();
      
      pos = permute(T.interpolateInvisible(res,'centroid'),[1,3,2]);
      bbox = T.interpolateInvisible(res,'bbox');
      [screenBoundingBox,xyframe] = getCalibrationBox(pos,bbox,T.videoHandler.frameSize);
      
      if plotif
        figure;
        imagesc(T.videoHandler.getCurrentFrame());
        colormap('gray')
        hold on;
        plot(xyframe(:,:,1),xyframe(:,:,2),'linewidth',1);
        rectangle('position',screenBoundingBox,'linewidth',1, ...
                  'edgecolor','r','facecolor','none')
        title('Estimated Screen size');
      end

      xy.helper.verbose('Found Bounding Box [%d,%d,%d,%d]',round(screenBoundingBox));

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
        if ~xy.helper.hasOpenCV() || ~self.useOpenCV
          relaese(self.videoWriter);
        else
          v = self.videoWriter;
          clear v;
        end
        self.videoWriter = [];
        self.writefile = '';
      end
      
    end
    
    function stepVideoWriter(self,frame)
      if ~isempty(self.videoWriter)
        if xy.helper.hasOpenCV() && self.useOpenCV
          self.videoWriter.write(frame);
        else
          self.videoWriter.step(frame);
        end
      end
    end
    
    
    function [writer] = newVideoWriter(self,vidname)
      xy.helper.verbose('Init VideoWriter "%s"',vidname);
      if (xy.helper.hasOpenCV() && self.useOpenCV) || ~self.hasVision 
        writer = cv.VideoWriter(vidname,self.videoHandler.frameSize([2,1]),'FPS',self.videoHandler.frameRate,'fourcc','X264');
      else
        writer = vision.VideoFileWriter(vidname,'FrameRate',self.videoHandler.frameRate);
      end
    end
    
    
    function [player] = newVideoPlayer(self,vidname)
      
      if self.hasVision
        player = vision.VideoPlayer();
      else
        player = xy.core.VideoPlayer();
      end
      
    end
    
    
    function [handler, timerange] = newVideoHandler(self,vidname,timerange,opts)
    % note: does not set the current time. 
      if nargin<3
        timerange = [];
      end
      if nargin<4
        opts = self.opts;
      end
      if self.useMex && xy.core.VideoHandlerMex.installed()
        self.useMex = 1;
        self.useOpenCV = 1;
        handler = xy.core.VideoHandlerMex(vidname,timerange,self.useKNN,opts);
      elseif xy.helper.hasOpenCV() && self.useOpenCV
        self.useMex = 0;
        self.useOpenCV = 1;
        handler = xy.core.VideoHandler(vidname,timerange,self.useKNN,opts);
      else
        self.useMex = 0;
        self.useOpenCV = 0;
        handler = xy.core.VideoHandlerMatlab(vidname,timerange,self.useKNN,opts);
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
      self.identityClassifier.plotif = ~~dopts.classifier && self.displayif;
      self.videoHandler.plotting(~~dopts.videoHandler && self.displayif);
                
      if self.stmif
        self.stimulusPresenter.progressBar = self.displayif && ...
            self.opts.display.stimulusProgress;
      end
          
    end
    
    
    function predictor = newPredictor(self,centroid)
    % Create a Kalman filter object.
      if self.hasVision
          vel = 1/(self.avgTimeScale)*self.bodylength;
          predictor = configureKalmanFilter('ConstantVelocity', centroid, [self.bodylength, self.bodylength], [vel,vel], 1);
      else
          predictor = [];
      end
      
    end
    

    
    
    function self=setupSystemObjects(self,vid)
    % Initialize Video I/O
    % Create objects for reading a video from a file, drawing the tracked
    % objects in each frame, and playing the video.
      
      if ~xy.Tracker.checkCompiled()
        error(['Please compile the code. Consult the Readme. ("make ' ...
               'clean;make" on linux)']);
      end
      %self.hasVision = ~isempty(which('vision.VideoPlayer'));
    
      
      %% Create a video file reader.
      self.videoHandler = [];
      [self.videoHandler,self.timerange] = self.newVideoHandler(vid,self.timerange);
      self.videoFile = vid;
      self.removeSaveFields();
      
      if self.stmif
        if isempty(self.opts.stimulus.presenter) 
          self.stimulusPresenter = xy.stimulus.Presenter(self.opts.stimulus);
        else
          if ischar(self.opts.stimulus.presenter)
            self.stimulusPresenter = eval(sprintf('%s(self.opts.stimulus)',self.opts.stimulus.presenter));
          elseif isa(self.opts.stimulus.presenter,'xy.stimulus.Presenter')
            self.stimulusPresenter = self.opts.stimulus.presenter;
          else
            error('Expect xy.stimulus.Presenter object');
          end
        end
        self.stimulusPresenter.setOpts(self.opts.stimulus);
        self.stimulusPresenter.init();
        self.addSaveFields('stmInfo','identityId','predIdentityId');
      end
      
      %% set all options
      self.initBackgroundAndObjects();
      self.checkOpts();
      self.setOpts();
      
    end
 
    
    function initBackgroundAndObjects(self)

      %% look for objects 
      if isempty(self.nindiv) || isempty(self.bodylength) || isempty(self.bodywidth) || self.useScaledFormat  

        [nindiv,bodySize] = self.findObjectSizes();
        
        
        if (nindiv>100 || ~nindiv) && isempty(self.nindiv)% too many... something wrong
          xy.helper.verbose('WARNING: The bodies size and number cannot be determined')
          if self.displayif && self.opts.display.bodySearchResults
            self.nindiv = xy.helper.chooseNIndiv(vid,1); % only if interactively
            bodySize = [100,20]; % wild guess;
            close(gcf);
          else
            error('Please manual provide bodylength, bodywidth and nindiv');
          end
        end
        
        
        if isempty(self.opts.bodylength) % otherwise already set by hand
          self.bodylength = bodySize(1);
        end
        if isempty(self.opts.bodywidth) 
          self.bodywidth = bodySize(2); 
        end
        if isempty(self.opts.nindiv) 
          self.nindiv = nindiv;
        end
        
      else
        % init background  
        if self.opts.tracks.initBackground
          n = min(self.videoHandler.history,floor(self.videoHandler.timeRange(2)*self.videoHandler.frameRate));
          n = min(n,500);
          self.videoHandler.reset(); % resets reader to timerange(1)
          self.videoHandler.resetBkg();
          self.videoHandler.initialize(0);
          self.videoHandler.computeSegments = false;
          
          for i = 1:n
            self.videoHandler.step();
            xy.helper.verbose('%1.1f%%\r',i/n*100); % some output
          end
          self.videoHandler.computeSegments = true;
          self.videoHandler.reset();
        end
        
      end
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
      self.videoHandler.bodylength = self.bodylength;
      self.videoHandler.bodywidth = self.bodywidth;
      self.videoHandler.headprop = self.opts.blob.headprop;

      self.timeStamp = self.videoHandler.getCurrentTime();
      
      if isempty(self.maxVelocity)
        self.maxVelocity = 15*self.avgVelocity;
      end
      
      self.avgTimeScale = 1/(self.avgVelocity/self.videoHandler.frameRate);
      xy.helper.verbose('Set time scale to %1.2f [frame/BL]',self.avgTimeScale);
      
      % initialize graph
      self.daGraph = [];
      self.daGraph = xy.core.DAGraph(self.nindiv,self.nindiv);

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
          xy.helper.verbose('WARNING: display tracks for writing a video file!');
        end
        self.videoWriter = [];
      end

      %% get new identity classifier 
      self.identityClassifier = newClassifier(self,self.opts.classifier);
      self.isInitClassifier = isInit(self.identityClassifier);

      
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

      self.identityId2TrackId = nan(250,self.nindiv);
      self.pos = [];
      self.res = [];
      self.tracks = [];
      self.savedTracks = [];
      self.savedTracksFull = [];
      self.leakyAvgFrame = [];
      

      if self.opts.tracks.keepFullTrackStruc
        xy.helper.verbose('WARNING: Will keep full tracks-structure. This will cost huge amount of memory!');
      end

      self.tracks = self.initializeTracks(); % Create an empty array of tracks.
      self.nextId = 1;
      self.uniqueIdentityFrames = 0;
      self.nonCrossingFrames = 0;

      self.meanCost(:) = 1;
      self.meanAssignmentCost = 1;

      self.segments = [];
      self.bboxes = [];
      self.centroids = [];
      self.locations = [];
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
      assert(~isempty(self.scalecost));

      FEATCOSTTYPES = {'Area','Size','BoundingBox'};
      self.featurecosttypes = {};
      for i = 1:length(FEATCOSTTYPES)
        if any(strcmp(FEATCOSTTYPES{i},self.costinfo));
          self.featurecosttypes(end+1) = FEATCOSTTYPES(i);
        end
      end
      
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
      self.nFramesForSingleUpdate = floor(2 * self.nFramesForUniqueUpdate);
      self.maxFramesPerBatch = 2*self.nFramesForSingleUpdate;

      self.meanNonCrossingFrames = self.nFramesForUniqueUpdate;
      
      if ~isinf(self.opts.classifier.timeToInit)
        self.nFramesForInit = min(max(ceil(self.opts.classifier.timeToInit*self.avgTimeScale),1),200);
      else
        self.nFramesForInit = Inf; % disables the classifier
      end

      xy.helper.verbose('nFramesForInit: %d',self.nFramesForInit);
      xy.helper.verbose('nFramesAfterCrossing: %d',self.nFramesAfterCrossing);
      xy.helper.verbose('nFramesForUniqueUpdate: %d',self.nFramesForUniqueUpdate);
      xy.helper.verbose('nFramesForSingleUpdate: %d',self.nFramesForSingleUpdate);

      xy.helper.verbose('clpMovAvgTau: %d',self.clpMovAvgTau);
      xy.helper.verbose('minBatchN: %d',self.minBatchN );

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
        
        while ~isempty(i) && length(self.segments)>self.nindiv

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
          overlap = getBBoxOverlap(self.bboxes,self.bboxes,self.bodywidth);
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
        
        % update location information
        self.locations = self.centroids;
        if ~self.useMex
          for i = 1:length(segm)
            % better take MSER regions if existent
            if ~isempty(segm(i).MSERregions)
              self.locations(i,:) = double(segm(i).MSERregionsOffset + segm(i).MSERregions(1).Location);
            end
          end
        else
        
          % use the mean of center line as centroid instead
          mcl = permute(sum(self.centerLine,1)/size(self.centerLine,1),[3,2,1]);
          ncl = ~isnan(mcl(:,1));
          self.locations(ncl,:) = mcl(ncl,:);
        end
        
        
        % features
        idfeat = permute(cat(4,segm.IdentityFeature),[4,1,2,3]);
        self.idfeatures =  idfeat;
        
        for ct = self.featurecosttypes
          self.features.(ct{1}) = cat(1,segm.(ct{1}));
        end
        
        
        self.classProb = predict(self.identityClassifier,self.idfeatures(:,:));
        self.classProbNoise = double(cat(1,segm.bendingStdValue));
      
      else
        self.centroids = [];
        self.locations = [];
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
    
    function cla = newClassifier(self,opts)
    % Create Classifier objects
      featdim = self.videoHandler.getFeatureSize();
      cla = xy.core.BatchClassifier(self.nindiv,featdim);
      
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

    
    function [nc,ncidx] = connectedComponents(self,trackIndices,assignedIdentityIds)
      oldIdentityIds = [self.tracks(trackIndices).identityId];
      A = zeros(self.nindiv);
      A(sub2ind(size(A),oldIdentityIds,assignedIdentityIds)) = 1;
      nc = xy.helper.networkComponents(A);
      % delete self-components
      nc(cellfun('length',nc)==1) = []; 
      
      % recode into the trackIndex
      ncidx = cell(size(nc));
      for i = 1:length(nc)
        ncidx{i} = find(ismember(oldIdentityIds,nc{i}));
      end
      
      % add those single components which are actually present
      idx = find(oldIdentityIds == assignedIdentityIds);
      for i = 1:length(idx)
        ncidx{end+1} = idx(i);
        nc{end+1} = oldIdentityIds(idx(i));
      end
      
    end
    
        
    
    function resetBatchIdx(self,trackidx)
      [self.tracks(trackidx).nextBatchIdx] = deal(1);
    end

    
    function changed = identityClassifierUpdate(self,trackIndices,updateif)
    % updates and tests the collected features according to the current identityIDs. It does NOT
    % force and update, update is only done if the batch is likly to come from the true
    % identity. In case a valid permutation is detected and reassignProb is higher than a
    % threshold, identityIDs are exchanged.
    %
    % resets the batchIdxs after update
      changed = 0;
      identityIds = cat(2,self.tracks(trackIndices).identityId);

      % test on the whole set of possible identityIds
      [assignedIdentityIds, prob, steps, probdiag] = self.predictIdentity(trackIndices,1:self.nindiv,...
                                                        self.nFramesForSingleUpdate);
      same = assignedIdentityIds==identityIds; 

      if updateif && (all(same) || self.enoughEvidenceForForcedUpdate(prob,steps,probdiag))
        % do not update if classes might be mixed up (wait for more data)
        % change classifier
        feat = self.getFeatureDataFromTracks(trackIndices);
        self.identityClassifier.batchUpdate(identityIds,feat,1); % force
        self.resetBatchIdx(trackIndices);
        self.uniqueIdentityFrames = 0;
        changed  = 1;
        return
      end
      
      %we here also check whether some unaccounted switching occurred (maybe outside of
      %a crossing)
      if length(same)>1 && ~all(same) && self.enoughEvidenceForAllIdentitySwitch(prob,steps,probdiag)
        % we have some mixed classes 

        % might be some NaNs...
        msk = ~isnan(assignedIdentityIds);
        trackIndices = trackIndices(msk);
        assignedIdentityIds = assignedIdentityIds(msk);
        prob = prob(msk);
        probdiag = probdiag(msk);
        steps = steps(msk);
        identityIds = cat(2,self.tracks(trackIndices).identityId);
        
        [nc,ncidx] = connectedComponents(self,trackIndices,assignedIdentityIds);
        
        for i = 1:length(ncidx)
          if length(ncidx{i})==1
            continue;
          end
          idx = ncidx{i};
          if all(ismember(nc{i},identityIds)) ...
              && (self.enoughEvidenceForReassignment(prob(idx),steps(idx),probdiag(idx)))

            % only permute if true permutation and enough evidence
            xy.helper.verbose('Switching identity with identityClassifierUpdate [%1.2f,%d]',min(prob(idx)),min(steps(idx)))

            switchIdentity(self,trackIndices(idx),assignedIdentityIds(idx),0);
            changed = 1;
            self.uniqueIdentityFrames = 0;
          end
        end
      end
      
    end

    
    function updatePos(self)
    % save current track->identityID
      t = self.currentFrame;
      identityIds = cat(1,self.tracks.identityId);
      tracks = self.tracks(~isnan(identityIds));

      if size(self.identityId2TrackId,1)< t
        self.identityId2TrackId = cat(1,self.identityId2TrackId,nan(self.nFramesExtendMemory,self.nindiv));
      end
      % save track positions
      if isempty(self.pos)
        self.pos = nan(2,self.nindiv,self.nFramesExtendMemory);
        self.tabs = nan(self.nFramesExtendMemory,1);
      end

      if size(self.pos,3)<t
        self.pos = cat(3,self.pos,nan(2,self.nindiv,self.nFramesExtendMemory));
      end

      if size(self.tabs,1)<t
        self.tabs = cat(1,self.tabs,nan(self.nFramesExtendMemory,1));
      end

      self.tabs(t,1) = self.timeStamp;
      
      if isempty(tracks) % no tracks
        return;
      end


      identityIds = cat(2,tracks.identityId);      
      trackIds = cat(2,tracks.id);      
      self.identityId2TrackId(t,identityIds) = trackIds;
      self.pos(1:2,identityIds,t) = cat(1,tracks.location)';
      
    end
    
    
    function success = initClassifier(self)
    % checks condition and inits classifier (if conds are met) and resets relevant counters

      success = false;
      thres = 3;
      if length(self.tracks)<self.nindiv 
        if  self.currentFrame > thres*self.nFramesForInit
          xy.helper.verbose('nindiv setting might be wrong or some identity are lost\r')
        end
        return
      end

      if self.uniqueIdentityFrames>self.nFramesForInit  || ...
          self.currentFrame> thres*self.nFramesForInit % cannot wait for ages..
        
        identityIds = cat(2,self.tracks.identityId);      

        feat = self.getFeatureDataFromTracks(1:length(self.tracks));
        [~,sidx] = sort(identityIds);
        feat = feat(sidx);
        batchsample = cell(1,self.nindiv);
        [batchsample{:}] = deal([]);
        batchsample(1:length(feat)) = feat;
        
        self.identityClassifier.init(batchsample);
        self.currentCrossBoxLengthScale = self.opts.tracks.crossBoxLengthScale; % update
        self.resetBatchIdx(1:length(self.tracks));
        
        % reset all potentials previous crossings
        [self.tracks.lastFrameOfCrossing] = deal(self.currentFrame);
        [self.tracks.firstFrameOfCrossing] = deal(self.currentFrame);
        [self.tracks.crossedTrackIds] = deal([]);
        
        % initialize DAG
        self.daGraph.reset(self.pos(:,:,self.currentFrame),self.currentFrame);
        
        
        self.uniqueIdentityFrames  = 0;
        success = true;
        
      elseif self.currentFrame < thres/2*self.nFramesForInit % cannot wait for ages..
        if ~isempty(self.crossings)
          self.resetBatchIdx(1:length(self.tracks));% we do not want a mixture at the beginning
        end
      end
      
      self.isInitClassifier = success;
      
    end
    
    function  [assignedIdentityIds, prob, steps, probdiag] = predictIdentity(self,trackIndices,identityIdSet,maxSteps)
    % predicts the identity according to the running probhist mean 
      
      C = zeros(length(trackIndices),length(identityIdSet));

      assignedIdentityIds = nan(size(trackIndices));
      prob = nan(size(assignedIdentityIds));
      probdiag = nan(size(assignedIdentityIds));
      steps = nan(size(assignedIdentityIds));

      oldIdentityIds = cat(1,self.tracks(trackIndices).identityId);
      for i = 1:length(trackIndices)
        track = self.tracks(trackIndices(i));
        if ~isempty(track.classProbHistory)
          nsteps = min(max(self.currentFrame-track.lastFrameOfCrossing,1),maxSteps);
          [p, nsteps] = track.classProbHistory.mean(nsteps);
          steps(i) = nsteps;        
        else
          p = track.clpMovAvg;
          steps(i) = maxSteps; % fake;
        end

        C(i,:) = 1-p(identityIdSet);
        probdiag(i) = p(oldIdentityIds(i));
      end
      C(isnan(C)) = Inf;
      
      assignments = xy.helper.assignDetectionsToTracks(C',2);
      
      assignedIdentityIds(assignments(:,2)) = identityIdSet(assignments(:,1)); % identityIDs have 1:nindiv order
      for i = 1:size(assignments,1)
        prob(assignments(i,2)) = 1-C(assignments(i,2),assignments(i,1));
      end
      prob(isinf(prob)) = nan;
    
    end
    
    function bool = enoughEvidenceForAllIdentitySwitch(self,prob,steps,probdiag)
      pother = self.meanClassificationProb-self.meanClassificationProbOther;
      if any(isnan(probdiag)) || pother<self.opts.classifier.minOtherProbDiff
        bool = false;
        if pother<self.opts.classifier.minOtherProbDiff
          xy.helper.verbose('Would like to switch but pother too low\r');
        end
        
      else
        p = max(prob-probdiag,0);
        n = sum(p>0);
        bool = n > 1 ...
               && sum(p)/n> self.opts.classifier.allSwitchProbThres*self.maxClassificationProb ...
               && min(steps)>=self.minBatchN...
               &&  min(steps) >= self.nFramesAfterCrossing;
        if bool
          xy.helper.verbose('Enough evidence: %1.2f',sum(prob-probdiag));
        end

      end
    end

    
    function bool = enoughEvidenceForReassignment(self,prob,steps,probdiag)
      if any(isnan(probdiag)) || sum(prob)/length(prob)<self.opts.tracks.probThresForIdentity*self.maxClassificationProb
        bool = false;
      else
        bool = max(prob-probdiag)>self.opts.classifier.reassignProbThres*self.maxClassificationProb ...
               && min(steps)>=self.minBatchN;

        if bool
          xy.helper.verbose('Enough evidence: %1.2f',sum(prob-probdiag));
        end

      end
    end
    
    function bool = enoughEvidenceForForcedUpdate(self,prob,steps,probdiag)
      if any(isnan(probdiag)) || (length(prob)<self.nindiv && length(prob)>1) ...
          || sum(prob)/length(prob)<self.opts.tracks.probThresForIdentity*self.maxClassificationProb
        bool = false;
      elseif length(prob)==1 %&& self.nindiv>self.MAXCROSSIDENTITY
        bool = true; % single update;
      else
        bool1 = min(steps) >= self.nFramesAfterCrossing;
        bool = bool1 && max(prob)<=self.opts.classifier.forcedUpdateProbThres*self.maxClassificationProb;
        bool = bool || (bool1 && self.meanClassificationProb-self.meanClassificationProbOther<self.opts.classifier.minOtherProbDiff);
      end
    end
    
    function bool = enoughEvidenceForBeingHandled(self,prob,steps,probdiag)
      if min(steps)<self.nFramesAfterCrossing || any(isnan(probdiag)) 
        bool = false;
      else
        m =sum(probdiag)/length(probdiag); 
        if m>self.opts.tracks.probThresForIdentity*self.maxClassificationProb
          
          bool = m>self.meanClassificationProb*self.opts.classifier.handledProbThresScale ...
                 || max(steps)>=self.nFramesForUniqueUpdate;
          % bool = bool || (self.nindiv>9 && length(prob)>self.nindiv/2); % crossing too big to be handled
        else
          bool = true; % not enough evidence, but to inaccurate to be sure,  release
        end
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
      identityIds = cat(2,self.tracks.identityId);      
      
      
      crossedTrackIdStrs = arrayfun(@(x)num2str(sort(x.crossedTrackIds)), self.tracks(trackIndices),'uni',0);
      [u,idxct,idxu] = unique(crossedTrackIdStrs);

      
      % we should distinguish the following cases
      % 1)  only one trackindex. (numel(trackIndices)==1)
      %    in this case, we should check :
      %      a) whether the track is very likely a certain
      %      identityID. If yes, we can erase it from the other tracks
      %      crossed track fields,
      %      otherwise, we should wait for the other tracks to
      %      finish. (for instance by just increasing the
      %      lastCrossed by one). [If this track does come into a
      %      new crossing, then all crossed tracks will also enter this new
      %      crossing and the tracks will be handled later]
      % 2)  if there are more than one trackindex. We test them in
      %     a permuted way, allowing for all identity currently in the
      %     crossing. If it is not clear which identity left the crossing,
      %     we just leave it and wait for the other to finish. If
      %      instead we find a true permutation of the identity that are
      %     leaving, we delete them out of the crossed fields of all
      %     others. 
      % 3)  if there is only one identity in the crossed id field, then,
      %     if it is the same trackid, it's OK. just delete it. If it
      %     is another trackid, we have a problem. One should force a
      %     all identity update (maybe we could just put ALL identity into the
      %     crossID field and let it handle the next time.)
      
      for i = 1:length(u)
        crossedTrackIds =  self.tracks(trackIndices(idxct(i))).crossedTrackIds;
        thisTrackIndices = trackIndices(idxu==i);

        
        tstart = min([self.tracks(thisTrackIndices).firstFrameOfCrossing]);
        tend = max([self.tracks(thisTrackIndices).lastFrameOfCrossing]);
        
        [~,crossedIdentityIds] = find(ismember(self.identityId2TrackId(tstart:tend,:),crossedTrackIds));
        crossedIdentityIds = unique(crossedIdentityIds)';

        assumedIdentityIds = identityIds(thisTrackIndices);
        
        %this should no happen because identityIDs swapping should be blocked during crossing
        assert(all(ismember(assumedIdentityIds,crossedIdentityIds)));
        

        % we force the permutation to be in the valid identity only (otherwise too many errors for many identity)
        [assignedIdentityIds, prob, steps, probdiag] = self.predictIdentity(thisTrackIndices,crossedIdentityIds,...
                                                          self.nFramesForUniqueUpdate); 

        MAXCROSSIDENTITY = max(10,self.nindiv/2);
        if length(assumedIdentityIds)>MAXCROSSIDENTITY && ~self.enoughEvidenceForReassignment(prob,steps,probdiag)
          % cannot deal with large crossings
          subDeleteCrossedTrackIds(thisTrackIndices); 
          continue;
        end

        OLD = 0;
        if OLD

          if all(ismember(assignedIdentityIds,assumedIdentityIds)) 
            % we have a valid assignment (note that it is always a permutation inside
            % those that cross because that is guaranteed in predictIdentity)

            
            [nc ncidx] = self.connectedComponents(thisTrackIndices,assignedIdentityIds);
            
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
                  xy.helper.verbose('valid permutation [%1.2f,%d]: switch...',min(prob(idx)),min(steps(idx)))
                  switchIdentity(self,thisTrackIndices(idx),assignedIdentityIds(idx),true)
                  subDeleteCrossedTrackIds(thisTrackIndices(idx)); % always handle
                  
                elseif self.enoughEvidenceForBeingHandled(prob(idx),steps(idx),probdiag(idx)) 
                  subDeleteCrossedTrackIds(thisTrackIndices(idx));
                end

              end
            end
            
          else
            if self.enoughEvidenceForBeingHandled(prob,steps,probdiag) 
              subDeleteCrossedTrackIds(thisTrackIndices);
            else
            
              % Still inside the identityID previously involved in the crossing. 
              % We better put
              % the track back into the crossing
              if ~any(isnan(prob))
                xy.helper.verbose('Put tracks back into crossings!')
                identityIds = [self.tracks.identityId];
                ids = [self.tracks.id];
                self.crossings{end +1 } = ids(find(ismember(identityIds,crossedIdentityIds)));
              end
            end
          end

            
        else % NEW (no exit of individual identity from crossings)

          if length(assignedIdentityIds)==1
            subDeleteCrossedTrackIds(thisTrackIndices); % always handle;
          elseif all(ismember(assignedIdentityIds,assumedIdentityIds)) 
            if any(assignedIdentityIds ~= assumedIdentityIds)
              
              if self.enoughEvidenceForReassignment(prob,steps,probdiag)
                % enough evidence switch tracks and delete from others. 
                xy.helper.verbose('valid permutation [%1.2f,%d]: switch...',min(prob),min(steps))

                switchIdentity(self,thisTrackIndices,assignedIdentityIds,true)
                subDeleteCrossedTrackIds(thisTrackIndices); % always handle;
                
              elseif self.enoughEvidenceForBeingHandled(prob,steps,probdiag) % might need to be handled because too long
                subDeleteCrossedTrackIds(thisTrackIndices);
              end
              
            else
              if  self.enoughEvidenceForBeingHandled(prob,steps,probdiag)
                subDeleteCrossedTrackIds(thisTrackIndices);
              end
              
            end
          else
            % should be always the case  ? 
            if self.enoughEvidenceForBeingHandled(prob,steps,probdiag) % might need to be handled because too long
              subDeleteCrossedTrackIds(thisTrackIndices);
            end
          end
        end
      end
      
      
      if self.displayif && self.opts.display.crossings
        plotCrossings_(self,2,trackIndices);
      end

      
      function subDeleteCrossedTrackIds(localTrackIndices)
        ids = [self.tracks.id];
        for ii = 1:length(localTrackIndices)
          
          thisTrackId = self.tracks(localTrackIndices(ii)).id;
          crossedIds = self.tracks(localTrackIndices(ii)).crossedTrackIds;
          crossedTrackIndices = find(ismember(ids,crossedIds));
          
          crossedIds = setdiff(crossedIds,thisTrackId);
          if isempty(crossedIds)
            crossedIds = [];
          end
          
          [self.tracks(crossedTrackIndices).crossedTrackIds] = deal(crossedIds); %#ok<FNDSB>
          self.tracks(localTrackIndices(ii)).crossedTrackIds = [];
        end

      end % nested function
        
    end
      
    function computeLeakyAvgFrame(self,frame)
      
      if isempty(self.leakyAvgFrame)
        self.leakyAvgFrame = im2double(frame);
      else
        p = self.opts.display.leakyAvgFrameTau;
        self.leakyAvgFrame = (1-1/p)*self.leakyAvgFrame + im2double(frame);
      end
    end
    
    
    function handled = testHandled(self)
      handled = false(1,length(self.tracks));
      for i = 1:length(self.tracks)
        handled(i) = isempty(self.tracks(i).crossedTrackIds);
      end
      
      %handled = arrayfun(@(x)isempty(x.crossedTrackIds),self.tracks);
    end
    
    function  hitBound = testBeyondBound(self)
    % hit bound will be true even if already handled!!
      hitBound = self.currentFrame-[self.tracks.lastFrameOfCrossing] >= self.nFramesAfterCrossing;
    end
    
    function newCrossing = testNewCrossing(self)
      newCrossing  = ismember(1:length(self.tracks), cat(2,self.crossings{:}));
    end
    
    function [bool,pdiff] = testTrackMisalignment(self,trackIndices)
    %      if length(trackIndices)<self.MAXCROSSIDENTITY
        identityIds = cat(1,self.tracks(trackIndices).identityId);
        clp = cat(1,self.tracks(trackIndices).clpMovAvg);
        [sz1,~] = size(clp);
        idx =  (identityIds-1)*sz1 + (1:sz1)';
        prob_correct = clp(idx); %  probability that current settings correct
        clp(idx) = 0;
        prob_notcorrect = max(clp,[],2); % maximal offset (might be nan)
        d = prob_notcorrect - prob_correct;
        n = sum(d>0);
        if n>1 % at least two
          pdiff = sum(d(d>0))/n;
          bool = pdiff > self.opts.classifier.allSwitchProbThres*self.maxClassificationProb;
        else
          pdiff = 0;
          bool = false;
        end
        
        
    % $$$       else
% $$$         bool = false;
% $$$         pdiff = 0;
% $$$       end
      
    end
    
    
    
    function updateClassifier(self)
    % updates the identity classifier and tracks the identity to
    % trackID assisgments
      
      self.crossings = self.detectTrackCrossings(); % only detects crossings


      mt = 20*self.nFramesForUniqueUpdate;
      self.meanNonCrossingFrames = self.meanNonCrossingFrames*(mt-1)/mt + 1/mt*self.nonCrossingFrames;

      % adjust update policy if necessary
      if self.currentFrame>mt
        if self.meanNonCrossingFrames<self.nFramesForUniqueUpdate/2 && self.nFramesForUniqueUpdate>2*self.minBatchN
          self.nFramesForUniqueUpdate = self.nFramesForUniqueUpdate-1;
          self.nFramesForSingleUpdate = self.nFramesForSingleUpdate-2;
        elseif self.meanNonCrossingFrames>self.nFramesForUniqueUpdate && self.nFramesForSingleUpdate<self.maxFramesPerBatch
          self.nFramesForUniqueUpdate = self.nFramesForUniqueUpdate+1;
          self.nFramesForSingleUpdate = self.nFramesForSingleUpdate+2;
        end
      end

      
      % update unique frames
      crossif = ~isempty(self.crossings);
      if ~crossif && length(self.tracks)==self.nindiv
        self.nonCrossingFrames = self.nonCrossingFrames + 1;
        self.uniqueIdentityFrames = self.uniqueIdentityFrames + 1;
      else
        self.nonCrossingFrames = 0;
        self.uniqueIdentityFrames = 0;
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
      if sum(handled)>1 && mod(self.currentFrame,2)
        handledIndices = find(handled);
        [misif,pdiff] = self.testTrackMisalignment(handledIndices);
        if misif 
          xy.helper.verbose('Probably misaligned (%1.2f)..\r',pdiff);
          self.identityClassifierUpdate(handledIndices,0); % use function for testing/switching 
        end
      end
      
          
      
      %% unique identity update
      if (self.uniqueIdentityFrames-1 > self.nFramesForUniqueUpdate)  && allhandled
        if self.identityClassifierUpdate(1:length(self.tracks),1);
          xy.helper.verbose('Performed unique frames update.\r')
        end
        self.uniqueIdentityFrames = 0; % anyway reset;
      end

      
      %% single identity update
      enoughData =  [self.tracks.nextBatchIdx] > self.nFramesForSingleUpdate;
      singleUpdateIdx = find(handled & enoughData);
      if ~isempty(singleUpdateIdx)
        if ~self.identityClassifierUpdate(singleUpdateIdx,1) 
          if  allhandled
            self.identityClassifierUpdate(1:length(self.tracks),0); % something might be wrong. Try all
          end
          self.resetBatchIdx(singleUpdateIdx); % anyway reset          
        else
          xy.helper.verbose('Performed single identity update.\r');          
        end
      end

    end

    
    function cleanup(self) 
    %switch if necessary
      
      trackIndices = 1:length(self.tracks);
      if ~isempty(trackIndices)
        self.identityClassifierUpdate(trackIndices,0);
      end

    end
    

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % Track related functions
    %---------------------------------------------------
    
    
    function tracks = initializeTracks(self) %#ok<MANU>
    % create an empty array of tracks
      tracks = struct(...
        'id', {}, ...
        'identityId', {}, ...
        'predIdentityId', {}, ...
        'bbox', {}, ...
        'centroid',{},...
        'location',{},...
        'velocity',{},... % in px/sec
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
        'consecutiveVisibleCount', {},...
        'switchedFrames',{},...
        'assignmentCost', {},...
        'nIdFeaturesLeftOut', {},...
        'stmInfo',{});
    end
    
    
    function predictNewLocationsOfTracks(self)

      if self.opts.tracks.kalmanFilterPredcition
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
          self.tracks(i).location = predictedCentroid;
        end
      elseif self.useMex && self.opts.tracks.centerLinePrediction
        % use centerLine head position
        vel = cat(1,self.tracks.velocity)*self.dt;
        if size(vel,1)~=self.nindiv
          return;
        end
        cl = cat(3,self.tracks.centerLine);
        d = squeeze(sqrt((cl(1,1,:)-cl(end,1,:)).^2 + (cl(1,2,:)-cl(end,2,:)).^2)); % better integrate along line in MEX?
        advancedClIdx = round(sqrt(vel(:,1).^2 + vel(:,2).^2)./d*size(cl,1));
        nextIdx =  min(max(floor(size(cl,1)/2)-advancedClIdx,1),size(cl,1));
        for i = find(advancedClIdx)'
          if ~isnan(cl(1,1,i))
            self.tracks(i).location = cl(nextIdx(i),:,i);
          end
        end
        
      end
    end
    
    
    function fcost = computeCostMatrix(self)
      
      nTracks = length(self.tracks);
      nDetections = size(self.locations, 1);
      fcost = zeros(nTracks, nDetections,length(self.costinfo));
      highCost = 5;
      somewhatCostly = 2;
      
      if nDetections>0
        maxDist = self.maxVelocity*self.dt*self.bodylength; % dist per frame
        thresDist = maxDist;

        
        % Compute the cost of assigning each detection to each track.
        for k = 1:length(self.costinfo)
          switch self.costinfo{k}

            case {'Location','Centroid'}

              %% center position 
              centers = cat(1,self.tracks.location);
              dst = xy.helper.pdist2Euclidean(centers,self.locations);

              dst1 = (dst/thresDist);

              fcost(:,:,k) = dst1;
            
            case {'Velocity'}
              
              centers = cat(3,self.tracks.location);
              velocity = cat(3,self.tracks.velocity);
              velmat = bsxfun(@minus,self.locations,centers)/self.dt;
              dvelmat = permute(sqrt(sum(bsxfun(@minus,velmat,velocity).^2,2)),[3,1,2]);
              fcost(:,:,k) = dvelmat;

            case {'CenterLine'}
              if self.useMex
                
                %% use the mex-file
                dst = xy.helper.pdist2CenterLine(cat(3,self.tracks.centerLine),self.centerLine);
                nidx = any(isnan(dst),2);
                dst = (dst/thresDist);
                dst(nidx,:) = 1;
                fcost(:,:,k) = dst;
                
              else
                fcost(:,:,k) = 1;
              end

              
            case {'Overlap'}

              bbsegs = double(cat(1,self.segments.BoundingBox));
              bbtracks = double(cat(1,self.tracks.bbox));
              overlap = getBBoxOverlap(bbtracks,bbsegs,0);
              fcost(:,:,k) = 1-overlap ;%+ (overlap==0)*somewhatCostly;
              
            case {'Classifier'}

              if self.isInitClassifier
                clProb = cat(1,self.tracks(:).classProb);
                dst = xy.helper.pdist2Euclidean(clProb,self.classProb);%correlation??!??
                dst(isnan(dst)) = 1; 
                
                % prob for identity
                msk = max(self.classProb,[],2)< self.opts.tracks.probThresForIdentity;
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
              fcost(:,:,k) = xy.helper.pdist2Euclidean(trackFeatures,detectedFeatures);
          end % switch

        end % costtype k loop
        
      end
    end
    
    
    
    
    function detectionToTrackAssignment(self)

      nDetections = size(self.locations, 1);

      self.cost = [];
      self.assignments= [];
      self.unassignedTracks = [];
      self.unassignedDetections = 1:nDetections;


      if nDetections && ~isempty(self.tracks)
        
        scalecost = self.scalecost/sum(self.scalecost);
        scale = scalecost./self.meanCost;

        % is there is currently a crossing then prob of non-assignement should be scaled down 
        costOfNonAssignment =  self.maxAssignmentCost + self.meanAssignmentCost*self.opts.tracks.costOfNonAssignment;

        % determine cost of each assigment to tracks
        fcost = self.computeCostMatrix();

        % scale the cost
        sfcost = sum(bsxfun(@times,fcost,permute(scale(:),[3,2,1])),3);

        % Solve the assignment problem. First for the visible
        nvis = cat(1,self.tracks.consecutiveInvisibleCount);
        validIdx = find(~nvis);%<self.avgTimeScale); !!!!!!!!!!!!!!!

        [self.assignments, self.unassignedTracks, self.unassignedDetections] = ...
            xy.helper.assignDetectionsToTracks(sfcost(validIdx,:),  costOfNonAssignment);

        self.unassignedTracks = validIdx(self.unassignedTracks);
        self.assignments(:,1) = validIdx(self.assignments(:,1)); 

        invisibleIdx = find(nvis);       

        if ~isempty(self.unassignedDetections) && ~isempty(invisibleIdx)
          % Solve the assignment problem for the invisible tracks
          sfcost1 = sfcost(invisibleIdx,self.unassignedDetections);

          costOfNonAssignment1 = costOfNonAssignment + nvis(invisibleIdx)*self.meanAssignmentCost;

          if self.opts.tracks.adjustCostDuringCrossing
            crossMsk = ~self.testHandled();
            crossMsk = crossMsk(invisibleIdx);
            costOfNonAssignment1(crossMsk) = costOfNonAssignment1(crossMsk)*self.opts.tracks.crossingCostScale;
          end

          [assignments, unassignedTracks, unassignedDetections] = ...
              xy.helper.assignDetectionsToTracks(sfcost1,  costOfNonAssignment1',max(costOfNonAssignment1)*ones(1,size(sfcost1,2)));

          self.assignments = [self.assignments;...
                              [invisibleIdx(assignments(:,1)),...
                              self.unassignedDetections(assignments(:,2))]];
          self.unassignedDetections = self.unassignedDetections(unassignedDetections);
          self.unassignedTracks = [self.unassignedTracks;invisibleIdx(unassignedTracks)];

        end
        

        % check max velocity of the assignments
        if ~isempty(self.assignments)
          detectionIdx = self.assignments(:, 2);
          locations = self.locations(detectionIdx, :);
          trackIdx = self.assignments(:, 1);
          latestLocations = cat(1,self.tracks(trackIdx).location);
          velocity = sqrt(sum((locations - latestLocations).^2,2))./(nvis(trackIdx) +1)/self.dt;
          
          idx = find(velocity>self.maxVelocity*self.bodylength);

          if ~isempty(idx) && self.dt
            self.unassignedTracks = [self.unassignedTracks; trackIdx(idx)];
            self.unassignedDetections = [self.unassignedDetections; detectionIdx(idx)];
            self.assignments(idx,:) = [];
            xy.helper.verbose('%d assignment(s) exceeded velocity.\r',length(idx));
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
        tmp = reshape(fcost,[],size(fcost,3));
        self.meanCost = self.meanCost*(mt-1)/mt  + 1/mt/size(tmp,1)*sum(tmp,1);

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
          xy.helper.verbose('there were %d unassigned detections\r',length(self.unassignedDetections))
        end

        
      end
      
    end
    
    
    
    function updateTracks(self)

    %% update assigned
      assignments = self.assignments;
      numAssignedTracks = size(assignments, 1);
      trackopts = self.opts.tracks;
      maxClass = 0;
      meanClass = 0;
      meanClassOther = 0;
      imeanClass = 0;
      
      for i = 1:numAssignedTracks
        trackIdx = assignments(i, 1);
        
        detectionIdx = assignments(i, 2);
        location = self.locations(detectionIdx, :);
        centroid = self.centroids(detectionIdx, :);

        %% Correct the estimate of the object's location  using the new detection.
        if trackopts.kalmanFilterPredcition
          correct(self.tracks(trackIdx).predictor, location); 
        end
        

        tau = trackopts.tauVelocity*self.avgTimeScale;
        self.tracks(trackIdx).velocity = self.tracks(trackIdx).velocity*(1-1/tau) ...
            + 1/tau*(self.tracks(trackIdx).location-location)/self.dt;
        
        %% detect whether the identity feature might be reversed
        reversed = 0;
        if self.useMex
          thickness = self.thickness(detectionIdx,:);
          centerLine = self.centerLine(:,:,detectionIdx);
          n2 = ceil(size(centerLine,1)/2);
          centeredCenterLine = bsxfun(@minus,centerLine,centerLine(n2,:));
          vel = self.tracks(trackIdx).velocity*self.dt; % already updated velocity 
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
          classprob = nan(1,self.nindiv);
          probNoise = NaN;
        end


        if reversed  && self.opts.classifier.discardPotentialReversed
          % not update.. identityfeature reversed (cannot be changed outside of mex)
          probNoise = NaN;
          classprob = nan(1,self.nindiv);
        else
          self.tracks(trackIdx).classProb =  classprob; % only update classprob in tracks if not reversed
        end

        % always update history, though (weight will be zero for NaN)
        [reasonable, w]  = self.tracks(trackIdx).classProbHistory.update(classprob, probNoise);

        
        if w>0
          %update moving average
          tmp = (1-1/self.clpMovAvgTau)*self.tracks(trackIdx).clpMovAvg;
          self.tracks(trackIdx).clpMovAvg =  tmp  + (w/self.clpMovAvgTau) * classprob;

          maxClass = max(maxClass,max(classprob));
          
          meanClass = meanClass + classprob(self.tracks(trackIdx).identityId);
          cl = classprob;
          cl(self.tracks(trackIdx).identityId) = 0;
          meanClassOther = meanClassOther + max(cl);
          imeanClass = imeanClass+1;
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
        
        self.tracks(trackIdx).assignmentCost =  self.cost(detectionIdx);
        
        %update locations
        self.tracks(trackIdx).centroid = centroid;
        self.tracks(trackIdx).location = location;

        
        %% Update visibility.
        self.tracks(trackIdx).totalVisibleCount =  self.tracks(trackIdx).totalVisibleCount + 1;

        % length of last visbility
        if self.tracks(trackIdx).consecutiveInvisibleCount>0
          self.tracks(trackIdx).consecutiveVisibleCount = 1;
        else
          self.tracks(trackIdx).consecutiveVisibleCount = self.tracks(trackIdx).consecutiveVisibleCount + 1;
        end
        
        self.tracks(trackIdx).consecutiveInvisibleCount = 0;
      end

      %% update max classprob for reassignment scaling
      if maxClass
        mt = min(self.opts.tracks.costtau*self.avgTimeScale,self.currentFrame-self.nFramesForInit);
        self.maxClassificationProb = self.maxClassificationProb*(mt-1)/mt + 1/mt*maxClass;
        self.meanClassificationProb = self.meanClassificationProb*(mt-1)/mt + 1/mt*meanClass/imeanClass;
        self.meanClassificationProbOther = self.meanClassificationProbOther*(mt-1)/mt + 1/mt*meanClassOther/imeanClass;
      end
      
      
      %% update unassigned
      for i = 1:length(self.unassignedTracks)
        ind = self.unassignedTracks(i);
        self.tracks(ind).age = self.tracks(ind).age + 1;
        self.tracks(ind).consecutiveInvisibleCount =  self.tracks(ind).consecutiveInvisibleCount + 1;

        % do not update the tracks with nan, but use clp average instead
        self.tracks(ind).classProb = self.tracks(ind).clpMovAvg;
        

        % only nan into history:
        self.tracks(ind).classProbHistory.update(nan(1,self.nindiv),NaN);
        
      end
      
      %% update DAG
      if self.isInitClassifier % classifier should be init
        self.daGraph.updateFromTracks(self.tracks,self.bodylength);

        
        %% predict identityIDs according to DAG
        currentPredIdentityIds = [self.tracks.predIdentityId];

        predIdentityIds = self.daGraph.predictIdentityIds();
        trackIndices = find(currentPredIdentityIds ~= predIdentityIds);
        if ~isempty(trackIndices)
          % switching!
          self.resetBatchIdx(trackIndices); 
          if self.verbosity>1 && self.nindiv<10
            xy.helper.verbose('Dag switched an identity...\r');
          end
        end
        for i = 1:length(trackIndices)
          self.tracks(trackIndices(i)).predIdentityId = predIdentityIds(trackIndices(i));
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
        identityIds = [self.tracks.identityId];
        if any(~isnan(identityIds(lostInds)))
          error('want to delete a true track');
        end
        
        self.tracks = self.tracks(~lostInds);

        % WE SHOULD PROB AVOID DELTEING TRACKS CURRENTLY CROSSING !
        %delete the track ids from the crossings
        for i = 1:length(self.tracks)
          self.tracks(i).crossedTrackIds = setdiff(self.tracks(i).crossedTrackIds,lostTrackIds);
        end
        
        xy.helper.verbose('Deleted %d tracks',sum(lostInds));
      end
      
    end
    
    
    function createNewTracks(self)
    % create new tracks for unassigndetections (if less the number of identity)

      if length(self.tracks)>=self.nindiv
        return
      end
      
      availableidentityids = setdiff(1:self.nindiv,cat(2,self.tracks.identityId));
      if  ~self.opts.tracks.withTrackDeletion && isempty(availableidentityids)
        return;
      end
      
      
      if self.isInitClassifier
        % already have identity.  take only one detection with least
        % assignment cost to the other classes
        unassignedDetections = self.unassignedDetections;
      else
        %take all at the beginning
        if isempty(availableidentityids)
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
        if s<=length(availableidentityids)
          newidentityid = availableidentityids(s);
        else
          newidentityid = NaN; 
          % just take the first detections. Maybe sort which to take?
          if ~self.opts.tracks.withTrackDeletion
            break
          end
        end
        
        if self.opts.tracks.kalmanFilterPredcition
          pred = self.newPredictor(self.locations(i,:));
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
          clprob = nan(1,self.nindiv);
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
        
        history = xy.core.ClassProbHistory(self.nindiv);
        % update the classprobhistory and set the parameters
        history.lambda = self.bodylength/self.bodywidth;
        [~,w] = history.update(clprob,clprobnoise);

        
        % Create a new track.
        newTrack = struct(...
          'id',        self.nextId,       ...
          'identityId',   newidentityid,...
          'predIdentityId',   newidentityid,...
          'bbox',      self.bboxes(i,:),  ...
          'centroid',  self.centroids(i,:),  ...
          'location',  self.locations(i,:),  ...
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
          'consecutiveVisibleCount', 1, ...
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

      locations = cat(1,self.tracks.location);
      dmsk = xy.helper.pdist2Euclidean(locations,locations)<self.currentCrossBoxLengthScale*self.bodylength;

      bbox = cat(1,self.tracks.bbox);
      bBoxXOverlap = bsxfun(@ge,bbox(:,1) + bbox(:,3), bbox(:,1)') & bsxfun(@lt,bbox(:,1), bbox(:,1)');
      bBoxYOverlap = bsxfun(@ge,bbox(:,2) + bbox(:,4), bbox(:,2)') & bsxfun(@lt,bbox(:,2), bbox(:,2)');
      bBoxOverlap = bBoxXOverlap & bBoxYOverlap;
      bBoxOverlap = bBoxOverlap | bBoxOverlap';
      
      costThres = self.opts.classifier.crossCostThresScale*self.meanAssignmentCost;
      %costThres = self.meanAssignmentCost;
      
      assignmentCost = [self.tracks.assignmentCost];
      nvis = [self.tracks.consecutiveInvisibleCount];
      vmsk = (assignmentCost>costThres) | (nvis>0);
      
      if self.nindiv<=5
        crossmat = bsxfun(@and,(bBoxOverlap | dmsk), vmsk);
      else
        crossmat = (bBoxOverlap | dmsk) & bsxfun(@and,vmsk, vmsk'); % only when all above thres
      end
      
      if any(crossmat(:))
        %crossmat = crossmat | crossmat'; % will be done with in the networkComponent function
        crossings = xy.helper.networkComponents(crossmat);
        crossings(cellfun('length',crossings)==1) = []; % delete self-components
      end
      

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
      

      self.predictNewLocationsOfTracks();
      
      self.detectionToTrackAssignment();
      
      self.updateTracks();
      self.updatePos();
      self.updateClassifier();

      if self.nindiv>length(self.tracks)
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
        xy.helper.verbose('WARNING: DO NOT CONVERT SAVEDTRACKS')
        return
      end
      
      s2mat = xy.helper.strucarr2strucmat(savedTracks);
     
      % no CAT yet !! too slow for large vectors (memory is copied). Use
      % cell temporarily, then memory block does not have to be continous
      if isempty(self.savedTracks)
        for f = fieldnames(s2mat)'
          self.savedTracks(1).(f{1}) = {};
          self.savedTracks.(f{1}){1} = s2mat.(f{1});
        end
      else
        for f = fieldnames(s2mat)'
          
          self.savedTracks.(f{1}){end+1} = s2mat.(f{1});
        end
      end
      savedTracks(:) =  [];
    end
    
    
    
    function trackinfo = getCurrentTracks(self)
    % gets the track info 
      
      if self.opts.tracks.keepFullTrackStruc && length(self.tracks)==self.nindiv 
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
      trackinfo = getCurrentTracks_(self.nindiv,self.tracks,self.saveFieldsIn,self.saveFieldsOut,self.saveFieldSegIf);

    end
    
    
    function checkVideoHandler(self)
    % checks wether the videoHandler is still OK
      
      
      if ~iscell(self.videoFile) && ~exist(self.videoFile)
        xy.helper.verbose('Cannot find videofile: %s',self.videoFile);
        self.videoFile = xy.helper.getVideoFile(self.videoFile);
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
        savename = 'TrackerSave';
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
    % XY.TRACKER(VID,...) starts a Tracker objects on a given video file
    % VID. If VID==[] a uidialog is opened to select the video file. 
    % 
    % Options can the set with XY.TRACKER(VID,OPTNAME1,OPTVALUE1,...) or
    % XY.TRACKER(VID,OPTS) where OPTS is a structure with fields identical
    % to the names of the options.
    %
    % A number of options are given in "avgBLC" which is average body
    % length covered and is a measure of the time it takes to cover
    % the distance with the avgVelocity.
    %  
    % To see information about the possible options, run (without arguments):
    %  >> xy.Tracker
    %
    % Example:
    %  >> T = xy.Tracker([],'nindiv',3);
    %  >> T.track([0,20]); 
    %  >> T.plot();



      def.vid = [];
      doc.vid = 'Video file name or {CAMIDX,''WRITEFILENAME''}';
      
      %% properties
      
      def.opts.nindiv = [];
      doc.nindiv = {'Number of individuals. Needs to be fixed ' '(*attempt* to be estimated if empty)'};
      
      def.opts.bodylength = [];
      doc.bodylength = {'Approx length of individuals in pixel (estimated if empty)'};
      def.opts.bodywidth = [];
      doc.bodywidth = {'Approx width of individuals in pixel (estimated if empty)'};

      def.opts.maxVelocity = [];
      doc.maxVelocity = {'Maximal velocity in BL/sec (estimated if empty)'};

      def.opts.avgVelocity = 4;
      doc.maxVelocity = {'Approx. avg velocity [BL/sec]. Important parameter ', ...
                          'to set the time scale of the tracking problem'};
      
      def.opts.useScaledFormat = false;
      doc.useScaledFormat = {'Use adaptive scaled gray format (EXPERIMENTAL)'};

      def.opts.timerange  = [];
      doc.timerange = {'Restrict tracking to time range [in seconds].'};
      
      def.opts.stmif = false;
      doc.stmif = 'Use online visual stimulation';

      def.opts.useMex = true;
      doc.useMex = {'Use the C++ version VideoHandler (FAST)'};

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
      doc.cost.Location = {'Location comparison','Relative cost weighting. 0 turns cost type off.'};

      def.opts.cost.Velocity = 3;
      doc.cost.Location = {'Velocity direction comparison','Relative cost weighting. 0 turns cost type off.'};
      
      def.opts.cost.Overlap = 2;
      doc.cost.Overlap = {'BBox overlap'};

      def.opts.cost.CenterLine = 5;
      doc.cost.CenterLine = {'CenterLine feature comparison'};

      def.opts.cost.Classifier = 0.1;
      doc.cost.Classifier = {'Classification prob comparison'};

      %      def.opts.cost.Size = 0; % !!! seems to have bad effect!
      %doc.cost.Size = {'Major/Minor axis comparison'};
      
      % def.opts.cost.Area = 0; % !!! seems to have bad effect!

      %doc.cost.Area = {'blob pixel area comparison'};

      def.opts.cost.BoundingBox = 0;
      doc.cost.BoundingBox = {'[x,y,w,h] bbox comparison',''};

      
      %% detector options
      def.opts.detector(1).history = 500;  %250 [nframes]
      doc.detector(1).history = 'Background update time constant [nFrames]';
      
      def.opts.detector.inverted = false;  
      doc.detector.inverted = {'Set 1 for IR videos (white bodies on dark background)'};
      
      def.opts.detector.adjustThresScale = 1;   
      doc.detector.adjustThresScale = {'0.8..1.2 : change thres when detections too noisy (useKNN=0)',''};


      %% reader
      def.opts.reader(1).resizeif = false;
      doc.reader(1).resizeif = {'Whether to resize the frame before tracking'};

      def.opts.reader.resizescale = 0.75; 
      doc.reader.resizescale = {'Fraction to resizing the frame for RESIZEIF = true',''};
      
      %% blob anaylser

      def.opts.blob(1).colorfeature = true; 
      doc.blob.colorfeature = {'Use color information for identityfeature'};

      def.opts.blob.headprop = 0.6; 
      doc.blob.headprop = {'Proportion of object length to use as feature',''};

      
      %% classification
      def.opts.classifier.npca = 40; 
      doc.classifier.npca = 'Number of PCA components';

      def.opts.classifier.nlfd = 10; 
      doc.classifier.nlfd = {'Number of LFD components. Set to 0 to turn off.',...
                          '-1 for auto setting.'};

      def.opts.classifier.tau = 5000; 
      doc.classifier.tau = {'Slow time constant of classifier [nFrames].'};

      def.opts.classifier.reassignProbThres = 0.2; %0.2%0.45
      doc.classifier.reassignProbThres = {'minimal diff probability for reassignments'};

      def.opts.classifier.allSwitchProbThres = 0.4; %0.6
      doc.classifier.allSwitchProbThres = {['minimal diff probability for ' ...
                          'all identity reassignments']};

      def.opts.classifier.forcedUpdateProbThres = 0.4; %0.45
      doc.classifier.reassignProbThres = {'minimal diff probability for reassignments'};

      def.opts.classifier.minOtherProbDiff = 0.2;
      doc.classifier.minOtherProbDiff = {'minimal runingavg class diff probability',...
                          'to use all identity switches and updates'};

      def.opts.classifier.handledProbThresScale = 0.5;
      doc.classifier.handledProbThresScale = {'mean class probability SCALE for crossing exits'};


      def.opts.classifier.crossCostThresScale = 2; 
      doc.classifier.crossCostThresScale = {'candidates for crossings: scales mean assignment cost',''};


      %% tracks
      def.opts(1).tracks.costtau = 500;
      doc.tracks.costtau = {'Time constant for computing the mean cost [avgBLC]'};

      def.opts.tracks.crossBoxLengthScale = 1; 
      doc.tracks.crossBoxLengthScale = {'How many times the bbox is regarded as a crossing'};
      
      def.opts.tracks.adjustCostDuringCrossing = true; 
      doc.tracks.adjustCostDuringCrossing = {'Whether to scale non-assignment cost during crossings'};

      def.opts.tracks.keepFullTrackStruc = false;
      doc.tracks.keepFullTrackStruc = {'Whether to keep full track structure. ONLY FOR DEBUG!'};
      
      def.opts.tracks.costOfNonAssignment = 2;
      doc.tracks.costOfNonAssignment =  {'Scales the threshold for','cost of non assignment'};

      def.opts.tracks.crossingCostScale =  1;
      doc.tracks.crossingCostScale = {'Scaling of non-assignment cost during crossings'};

      def.opts.tracks.probThresForIdentity = 0.05; 
      doc.tracks.probThresForIdentity = {'Classification probability to ' 'assume a identity feature'};

      def.opts.tracks.useDagResults = 0;
      doc.tracks.useDagResults = {'Sets default output results to ' 'DAG (1) or Switch (0) method',''};
      
      def.opts.tracks.kalmanFilterPredcition = false; 
      doc.tracks.kalmanFilterPredcition = {'Whether to use Kalman filter (overwrites CL prediction)'};

      def.opts.tracks.centerLinePrediction = true; 
      doc.tracks.centerLinePrediction = {'Whether to use prediction along center line'};

      
      def.opts.classifier.timeToInit = 20;  
      doc.classifier.timeToInit = {'Time to initialize the classifier [avgBLC]'};

      def.opts.classifier.timeAfterCrossing =  1; 
      doc.classifier.timeAfterCrossing = {['When to check for permutations after ' ...
                          'crossings [avgBLC]']};
      
      def.opts.classifier.timeForUniqueUpdate = 12; 
      doc.classifier.timeForUniqueUpdate = {['Unique frames needed for update all ' ...
                          'identity simultaneously [avgBLC]']};

      def.opts.classifier.clpMovAvgTau = 0.5; 
      doc.classifier.clpMovAvgTau = {'Time constant of class prob','moving average [avgBLC].'};

      def.opts.tracks.tauVelocity = 0.5; 
      doc.tracks.tauVelocity = {'Time constant to compute the ' 'velocity [avgBLC]'};

      def.opts.tracks.initBackground = 1;
      
      %% dag
      def.opts.dag.probScale = 0.5;
      doc.dag.probScale = {'DAGraph probScale. 1 means 50/50 weighting of ', ...
                          'classprob with distance if points a bodylength apart',''};
      
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

      def.opts.display.BWImg = false;
      doc.display.BWImg = {'Whether to plot the BWImg during the tracking process'};
      
      def.opts.display.level = 3;
      doc.display.level = {'Level of details of the track plot'};
      
      def.opts.display.bodySearchResults = false;
      doc.display.bodySearchResults = {'Info plot nindiv auto-search'};

      def.opts.display.stimulusProgress = true;
      doc.display.stimulusProgress = {'ProgressBar in case of stimulation'};
      
      def.opts.display.switchIdentity = false;
      doc.display.switchIdentity = {'Switch identity info plot (for DEBUGGING) '};

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
      
      MFILENAME = 'xy.Tracker';
      xy.helper.parseInputs;
      if HELP;self.videoHandler = [];return;end

      global VERBOSELEVEL;
      VERBOSELEVEL = self.verbosity;

      
      %%options depending on avgVelocity (will be set in init)

      opts.display.leakyAvgFrameTau = 50;
      doc.display.leakyAvgFrameTau = {'Tau for switch identity plot'};


      %%options not really important or debug
      opts.detector.fixedSize = 0;  
      doc.detector.fixedSize = {'Set 1 for saving the fixedSize image'};

      opts.detector.nskip = 5; 
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
      doc.blob.difffeature = {'Use background subtracted gray' 'images for identity feature'};

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
      
      
      %% parameter checking

      if isempty(vid)
        vid = xy.helper.getVideoFile();
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
        if length(vname)>1 && all(vname(1:2)=='~/') && isunix()
          [~,home] = unix('eval echo ~$USER');
          vname = [home(1:end-1) vname(2:end)];
        elseif vname(1)=='~'
          error('provide full path name. Cannot start with "~"');        
        end

        if ~exist(vname) && ~iscell(vid) && length(vid)>1
          error(sprintf('Video file "%s" not found',vname));
        end
      end
      
      if ~isempty(opts.nindiv) && opts.nindiv<0 && ~iscell(vid)
        opts.nindiv = xy.helper.chooseNIndiv(vname,1);
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
        fprintf('Empty xy.Tracker Object\n\n');
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
    %  >> T = xy.Tracker(videoFile,'nindiv',3,'displayif',2);
    %  >> T.track([0,100]);
    %  >> T.plot();
    %  >> T.save();
      
      if exist('writefile','var') && ~isempty(writefile)
        self.writefile = writefile;
      end

      if ~exist('saveif','var') || isempty(saveif)
        saveif = false;
      end
      
      if exist('trange','var') 
        self.timerange = trange;
      else
        self.timerange = self.videoHandler.timeRange; % take all
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
        xy.helper.verbose('Start tracking of %d identities in the time range [%1.0f %1.0f]...',...
                            self.nindiv,self.timerange(1),self.timerange(2));
      else
        xy.helper.verbose('Start tracking of %d identities with stimulation', self.nindiv);
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
        if self.displayif && self.opts.display.switchIdentity
          [self.segments, self.timeStamp, frame] = self.videoHandler.step();
          self.computeLeakyAvgFrame(frame);
        else
          [self.segments,self.timeStamp] = self.videoHandler.step(); % faster..
        end 
        self.dt = self.timeStamp -self.lastTimeStamp;
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
          savedTracks(1:self.nindiv,s) = self.getCurrentTracks();
        end

        if ~mod(s,self.nFramesAppendSavedTracks) 
          savedTracks = self.appendSavedTracks(savedTracks);
          s = size(savedTracks,2);
          
          if ~self.stmif &&  ~mod(self.currentFrame,self.nStepsSaveProgress*self.nFramesAppendSavedTracks) && saveif
            % occasionally save for long videos 
            xy.helper.verbose('Save tracking progress to disk..');
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
            xy.helper.verbose(['%1.0fh %1.0fm %1.1fs  [%1.2f x real ' ...
                                'time]             \r'],t(4),t(5),t(6),trackDuration/verboseDuration);
          else
            xy.helper.verbose(['%1.0fh %1.0fm %1.1fs  [%1.2f FPS] ' ...
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
        xy.helper.verbose('Save tracking results to disk...');
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




