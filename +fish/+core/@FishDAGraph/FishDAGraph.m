classdef FishDAGraph < handle;
% FISHDAGRAPH is a class for implementing the Directed Acyclic
% Graph based MHT tracking approach. It looks for the shorted path
% from the beginning of a crossing to the end, where identities are
% again known. It basically accumulates the classprobs (and
% distances) on the way and finds the shortest and likeliest
% (according to the classprobs of individual fish). 
%
% update is fast and linear in N = ndet*nhyp*nfish. per time step. 
% Backtracing to find the tracks needs a single through the past
% time steps for each fish. 
%
% CAUTION: it might be that individual detections are shared among
% fish.
  
  
 
    
  properties
    nhyp = [];
    dim = 2;
    probScale = 0.5; % 1 means 50/50 weighting of classprob with distance
                   % if points a fishLength apart
    useClassProbNoise = false; % weights class prob with estimated noise
    saveDagIif = 0;
  end

  properties (GetAccess=private)
    startPenalty = 1e5; % for start
    seqMissedPenalty = 100;
    memoryBlock = 250; % in frames
  
    dagIsave = [];
  end
  
  
  properties (SetAccess=private,GetAccess=private)
    nfish = [];
    frameCounter = 0;
    dagPos = [];
    dagIdx = [];
    dagI = [];
    
    offsetFrame = 0;

  end

  methods(Access=private)    


    function setOpts(self,opts)

      if isfield(opts,'tracks')
        opts = opts.tracks;
      end
      
      for f = fieldnames(opts)
        if isprop(self,f{1});
          self.(f{1}) = opts.(f{1});
        end
      end
    end
    
  end
  

  
  methods    
    
    function test(self);
      
      T = 0.3;
      dt = 0.004;
      nt = floor(T/dt)+1;
      velocity = 10;
      sig = 0.05;
      pfalsealarm = 0.0;
      pmissed = 0.00;
      background = 0.1;
      foreground = 0.12;

      nfish = self.nfish;
      nobj = self.nhyp;
      assert(self.dim==2,'Test only available for dim==2');
      
      cutoff = velocity*dt;
      [Hb,Ha] = butter(4,cutoff,'low');
      trace = filtfilt(Hb,Ha,randn(nt,2,nfish));
      trace = trace/std(trace(:));

      % noisy observations of the tracks
      otrace = trace + randn(size(trace))*sig;
      
      % add some noise detections
      if nfish<nobj
        nadd = nobj-nfish; 
        pos =randn(nt,self.dim,nadd);
        
        otrace = cat(3,otrace, pos);
        for i = nfish+1:nobj
          ridx = find(rand(nt,1)>pfalsealarm);
          otrace(ridx,:,i) = NaN;
        end
      end

      
      % some missed detections
      for j = 1:nobj
        ridx = find(rand(nt-1,1)<pmissed);
        otrace(ridx,:,j) = NaN;
      end
      

      br = background*rand(nt,nfish,nobj);
      classprob = br;
      for i = 1:nfish
        classprob(:,i,i) = classprob(:,i,i) + foreground*rand(nt,1,1);
      end
      for j = 1:nobj
        classprob(isnan(otrace(:,1,j)),:,j) = br(isnan(otrace(:,1,j)),:,j);
      end
      

      % distance for tracking
      dsig = sig;
      col = parula(nfish+1);
      sym = 'xos<>^vp';
      
      % permute the positions (obj index is random)
      [~,ridx] = sort(rand(nt,1,nobj),3);
      dim = 2;
      idx = sub2ind(size(otrace),repmat((1:nt)',[1,dim,nobj]),...
                    [ones(nt,1,nobj),2*ones(nt,1,nobj)],...
                    repmat(ridx,[1,dim,1]));
      rotrace = otrace(idx);
      
      idx = sub2ind(size(classprob),repmat((1:nt)',[1,nfish,nobj]),...
                    repmat(1:nfish,[nt,1,nobj]),...
                    repmat(ridx,[1,nfish,1]));
      rclassprob = classprob(idx);
      

      clf;
      for i = 1:nfish
        plot(squeeze(trace(:,1,i)),squeeze(trace(:,2,i)),'color',col(i,:));
        hold on;
      end

      plot(squeeze(rotrace(:,1,:)),squeeze(rotrace(:,2,:)),'k.');


      %tracking loop
      otrace = permute(otrace,[2,3,1]); % dim x nobj x  nt
      rotrace = permute(rotrace,[2,3,1]); % dim x nobj x  nt
      rotrace(:,:,nt) = otrace(:,:,nt); % set last in correct order !!

      rclassprob = permute(rclassprob,[2,3,1]); % dim x nobj x  nt
      classprob = permute(classprob,[2,3,1]); % dim x nobj x  nt
      lastcl = classprob(:,:,1);
      
      self.reset(otrace(:,:,1));

      h = [];
      for i_t = 2:nt
  
        detections = rotrace(:,:,i_t);
        cl = rclassprob(:,:,i_t);
  
        msk = all(isnan(detections),1);
        detections(:,msk) = [];
        cl(:,msk) = [];
        
        % compute the DAG shortest path
        self.update(cl,detections);
        
        %plot whole traces 
        ifish = 1;
        delete(h);
        h = self.plotTraces(ifish*ones(1,self.nhyp),1:self.nhyp);
        drawnow;

      end
      
      self.checkOverlap();
      self.plotTraces(1:self.nfish,1:self.nfish);

    end

    
    function varargout = plotTraceAssignments(self,fishIds,assignments,stepsback);
    % plots the traces for a current assignemts. Called from fish.Tracker
      
      detectionIdx = assignments(:, 2);
      fishIdx = fishIds(assignments(:, 1));
      
      h = self.plotTraces(fishIdx,detectionIdx,[],stepsback);
      if nargout
        varargout{1} = h;
      end
      
    end
    
      
    function varargout=plotTraces(self,ifish,mobj,mt,stepsback)
    % PLOTTRACES(SELF,IFISH,MOBJ,[MT]) plot traces for all
    % combinations of IFISH(i) and MOBJ(i), backtraced starting from time MT. 
    % 
    % NOTE: if ifish and ihyp have the same length ifish(i),
    % ihyp(i) is plotted instead

      col = parula(self.nfish+1);
      sym = 'xos<>^vp';
      
      if nargin<4
        mt = self.frameCounter;
      end
      if nargin<5
        stepsback = mt;
      end

      [rtrace,objidx] = self.backtrace(ifish,mobj,mt,stepsback);

      hold on;
      h = [];
      for i = 1:length(mobj)
        h(i) = plot(rtrace(1,:,i),rtrace(2,:,i),sym(mod(mobj(i)-1,length(sym))+1),'color',col(ifish(i),:));
      end
      
      if nargout
        varargout{1} = h;
      end
    
    end
    
    
    function varargout =checkOverlap(self,objidx,verbosity);

      if nargin<2 || isempty(objidx)
        [rtrace,objidx] = self.backtrace();
      end
      objidx = reshape(objidx,[],self.nfish,self.nhyp);

      if nargin<3
        verbosity = 2;
      end
      
      % check for potential overlaps
      n = self.nfish;
      eqmsk = zeros([size(objidx,1),n*(n-1)/2]);
      s = 0;allL = [];
      for i = 1:self.nfish
        for j = i+1:self.nfish
          s = s+1;
          eqmsk(:,s) = (objidx(:,i,i) == objidx(:,j,j));
          if sum(eqmsk(:,s))>1
            d = diff([0;eqmsk(:,s);0]);
            stop = find(d==-1);
            start = find(d==1);
            L = stop - start; 
            sL = sort(L,'descend');
            allL = [allL;L];
            if verbosity>1
              fish.helper.verbose('[%d,%d] Found %1.2f%% [%d] overlap. ', i,j,mean(eqmsk(:,s))*100,sum(eqmsk(:,s)))
              fish.helper.verbose('Median length %d. Longest: [ %s] \n', median(L),sprintf('%d ',sL(1:min(end,3))));
            end
          end
          
        end
      end

      if verbosity==1
        sallL = sort(allL,'descend');
        fish.helper.verbose('Found %1.2f%% [%d] overlap. ', mean(eqmsk(:))*100,sum(eqmsk(:)))
        fish.helper.verbose('Median length %d. Longest: [ %s]', median(allL),sprintf('%d ',sallL(1:min(end,8))));

      end
      
      if nargout
        varargout{1} = eqmsk;
      end
    end
    
    
    function res = getVariables(self)
      res = struct(self);
    end
    
    function [rtrace,objidx] = backtrace(self,ifish,mobj,mt,stepsback)
    % [RTRACE,OBJIDX] = BACKTRACE(SELF,IFISH,MOBJ,MT) backtraces for the
    % fishids IFISH by assuming the corresponding MOBJ(j) at time MT
    % belongs to IFISH(i).
    %
    % If IFISH is empty all possible combinations are returned.
      
      if nargin<4 || isempty(mt)
        mt = self.frameCounter;
      end
      mt = min(mt,self.frameCounter);
      
      if nargin<5 || isempty(stepsback)
        stepsback = mt;
      end
      if stepsback>mt
        stepsback = mt;
      end
      
      if nargin<2 ||  isempty(ifish) 
        [ifish,mobj] = ndgrid(1:self.nfish,1:self.nhyp);
        ifish = ifish(:)';
        mobj = mobj(:)';
      end
      
      assert(length(mobj)==length(ifish),['IFISH and MOBJ has to be of same length']);
      
      [rtrace,objidx] = backtrace_(self.dagIdx,self.dagPos,ifish,mobj,mt,stepsback);
      
      DEBUG = 0;
      if DEBUG
        % test
        n = length(ifish);
        objidxt = zeros(mt,n);
        rtracet = zeros(self.dim,mt,n);
        
        for i = 1:n
          
          rtracet(:,mt,i) = self.dagPos(:,mobj(i),mt);
          objidxt(mt,i) = mobj(i);
          
          for t = mt-1:-1:1
            objidxt(t,i) = self.dagIdx(ifish(i),objidxt(t+1,i),t+1);
            rtracet(:,t,i) = self.dagPos(:,objidxt(t,i),t);
          end
          
        end
        
        assert(all(objidx(:)==objidxt(:)));
        assert(all(rtrace(:)==rtracet(:)));
        
      end
    end
    
    
    function predFishIds = predictFishIds(self);

      assignments = fish.helper.assignDetectionsToTracks(self.dagI,2*max(self.dagI(:)));
      predFishIds = assignments(assignments(:,1),2)';
    end
    
    
    
    function reset(self,initPos,tframe);
      % init pos has to be in fishID order...
      
      self.dagI = self.startPenalty*ones(self.nfish,self.nhyp);
      self.dagI(find(eye(self.nfish))) = 0;
      
      self.dagIdx = nan(self.nfish,self.nhyp,self.memoryBlock);
      self.dagPos = nan(self.dim,self.nhyp,self.memoryBlock);

      if self.saveDagIif
        self.dagIsave = nan(self.nfish,self.nhyp,self.memoryBlock);
        self.dagIsave(:,:,1) = self.dagI;
      end
      
      self.frameCounter = 1; % next step if 2

      for i = 1:self.nfish
        self.dagIdx(i,:,1) = i;
      end

      if nargin>1
        assert(size(initPos,1)==self.dim,'Dimension of position mismatch!')
        self.dagPos(:,1:size(initPos,2),1) = initPos;
      end
      
      if nargin>2
        self.offsetFrame = tframe;
      end
    end
  
    
    function updateFromTracks(self,tracks,fishLength)
    % UPDATEFROMTRACKS(SELF,TRACKS) as UPDATE but takes position
    % information directly from tracks strcuture. Length of tracks
    % should be equal to the number of hypothesis. Expects position
    % information for all tracks every frame (could be a predicted value)

      pos = cat(1,tracks.centroid)'; 
      cl = cat(1,tracks.classProb)';

      if size(cl,2)~=size(pos,2)
        % some classprob were empty. Do one by one
        cl = zeros(size(pos,2));
        for i = 1:length(tracks)
          if ~isempty(tracks(i).classProb)
            cl(:,i) = tracks(i).classProb;
          end
        end
      end
      
      cl(isnan(cl)) = 0;
      if self.useClassProbNoise
        w = cat(1,tracks.classProbW)';
        cl = bsxfun(@times,cl,w);
      end
      
      dist = fish.helper.pdist2Euclidean(self.dagPos(:,:,self.frameCounter)',pos')/fishLength;

      self.updateWCost(dist,cl,pos);
    
    end
    
    
    
    function update(self,classProb,detections,costOfNonAssignment)
    % UPDATE(SELF,CLASSPROB,DETECTIONS,[COSTOFNONASSIGNMENT]) computes the
    % Directed Acyclic Graph minimal path problem in an online manner
    % (Nodes are in topological order) based on distance and class
    % probs. DETECTIONS are the positions/features of the detections
    % in format NDIM x NDETECT. CLASSPROB is NFISH x NDETECT

        
      dist = fish.helper.pdist2Euclidean(self.dagPos(:,:,self.frameCounter)',detections');

      if nargin<4 || isempty(costOfNonAssignment)
        costOfNonAssignment = nanmedian(dist(:));
      end
      
      
      dist(all(isnan(dist),2),:) = costOfNonAssignment*self.seqMissedPenalty; % missed second expensive
      dist(isnan(dist)) = costOfNonAssignment; % missed one cheap

      self.updateWCost(dist,classProb,detections,costOfNonAssignment);
    
    end
    
    
    
    function updateWCost(self,sfcost,classProb,detections,costOfNonAssignment);
    % computes the Directed Acyclic Graph minimal path problem in
    % an online manner (Nodes are in topological order). SFCOST is
    % cost matrix of size NHYP x NDETECT
    % DETECTIONS shoue be NDIM x NDETECT. CLASSPROB is NFISH x
    % NDETECT 
      
      ndetect = size(sfcost,2);

      % construct correct sized matrizes
      if ndetect==self.nhyp
        dist = sfcost;
        cl = classProb;
        pos = detections;
      elseif ndetect<self.nhyp
        dist = costOfNonAssignment*ones(self.nhyp,self.nhyp);
        dist(:,1:ndetect) = sfcost;
        cl = zeros(self.nfish,self.nhyp);
        if ~isempty(classProb)
          cl(:,1:ndetect) = classProb;
        end
        pos = nan(2,self.nhyp);
        pos(:,1:ndetect) = detections;
      else % ndetect>nhyp 
        [~,rank] =  sort(min(sfcost,[],1),'ascend');
        msk = rank<=nhyp;
        dist = sfcost(:,msk);
        if ~isempty(classProb)
          cl = classProb(:,msk);
        end
        pos = detections(:,msk);
      end
      
      % compute min path
      dist = bsxfun(@minus,dist,min(dist,[],1));
      [tmpdagI,idx] = min(bsxfun(@plus,self.dagI,permute(dist,[3,1,2])),[],2);
      
      self.dagI = permute(tmpdagI,[1,3,2]);
      % add cl to distance if available
      if ~isempty(cl)
        self.dagI = self.dagI + (1-cl)*self.probScale;
      end

      % save
      self.frameCounter = self.frameCounter + 1;

      if size(self.dagIdx,3)<self.frameCounter
        self.dagIdx = cat(3,self.dagIdx,nan(self.nfish,self.nhyp,self.memoryBlock));
        
        self.dagPos = cat(3,self.dagPos,nan(2,self.nhyp,self.memoryBlock));        

        % to avoid overflow
        self.dagI = self.dagI - min(self.dagI(:));

        if self.saveDagIif
          self.dagIsave = cat(3,self.dagIsave,nan(self.nfish,self.nhyp,self.memoryBlock));          
        end
      
      end
      
      self.dagIdx(:,:,self.frameCounter) = permute(idx,[1,3,2]);
      self.dagPos(:,:,self.frameCounter) = pos;

      if self.saveDagIif
        self.dagIsave(:,:,self.frameCounter) = self.dagI;
      end

    end

    function dagI = getSavedDagI(self);
      if ~self.saveDagIif
        error('Need to turn on option "saveDagIif" to save the dagI');
      else
        dagI = self.dagIsave;
      end
    end
    
    function self = FishDAGraph(nfish,nhyp,opts)
    % self = FishDAGraph(nfish,nhyp)
      self.nhyp = nhyp;
      self.nfish = nfish;
      
      if nargin>2
        self.setOpts(opts);
      end
      self.reset();
      
    end
    
  end
end
