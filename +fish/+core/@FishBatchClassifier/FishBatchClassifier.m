classdef FishBatchClassifier < handle;
% Classifies the Fish 

  properties 
    npca = 15;
    nlfd = 0;
    nfish = 2;
    minBatchN = 50; % minimal number of features for updateing the batch
    noveltyThres = 0.10;
    regularizer = 1e-10;
    featdim = [];
    outliersif = 0;
    plotif = 0;
    tau = Inf; % in fact 1/tau in frames
  end

  properties (SetAccess = private)

    Sigma = [];
    invSigma = [];
    mu = [];
    norm = [];
    n = [];
    d = [];
    pcapc = [];
    lfdpc = [];
    pcamu = [];
  end
  
  
  methods(Static)
    function obj = loadobj(S);
      if isstruct(S)

        obj = fish.core.FishBatchClassifier(S.nfish,S.featdim);
        for f = fieldnames(S)'
          if isprop(obj,f{1}) 
            obj.(f{1}) = S.(f{1});
          end
        end
      else
        obj = S;
      end
    end
  end
  
  
  
  methods (Access=private);
    
    
    
    function initialize(self,batchsample)
    % gets the data from the first sample. 
    % CAUTION: also computes the pca from it....

      
      if ~iscell(batchsample) || length(batchsample)~=self.nfish || ...
            all(cellfun('isempty',batchsample))
        fish.helper.verbose('WARNING: provide a batch for each fish!');
        keyboard
        return
      end
      
      
      if prod(self.featdim)==size(batchsample{1},2)
      elseif 3*prod(self.featdim)==size(batchsample{1},2)
        self.featdim = [self.featdim ,3];        
      else
        fish.helper.verbose('WARNING: Set feature dimensions correctly, otherwise plotting will be impiared');
        self.featdim = [size(batchsample{1},2),1];        
      end
      
      
      fish.helper.verbose('Initialize FishBatchClassifier...');
      
      if self.npca>0 || self.nlfd>0
        Z = cat(1,batchsample{:});
        % do some weak oulyier detection
        [pc,a] = fish.helper.emclustering(Z,self.nfish+2);
        [~,cl] = max(pc,[],2);
        idx = find(a<self.noveltyThres/self.nfish);
        goodidx = ~ismember(cl,idx);
        
        cl = cellfun(@(x,k)k*ones(size(x,1),1),batchsample,num2cell(1:length(batchsample)),'uni',0);
        cl = cat(1,cl{:});
        [self.pcapc,~,proj,self.pcamu] = fish.helper.pca1(Z(goodidx,:),max(self.nlfd,self.npca));
        self.npca = max(self.nlfd,self.npca);

        if self.nlfd>0 && self.npca>self.nlfd
          [self.lfdpc,~,proj] = fish.helper.lfd(proj,cl(goodidx),self.nlfd,0);
        end
          
        self.pcamu = self.pcamu';
        d = size(proj,2);
      else
        d = size(batchsample{1},2);
      end

      % initialize 
      self.d = d;
      self.n = zeros(self.nfish,1);
      self.mu = zeros(self.nfish,d);
      self.Sigma = repmat(eye(d),[1,1,self.nfish]);
      self.invSigma = repmat(eye(d),[1,1,self.nfish]); % dummy
      self.norm = ones(self.nfish,1);
      
      % convert the batchsample
      l = cellfun('isempty',batchsample);
      X = self.convertBatch(batchsample(~l));

      % initial setting
      s = 0;
      for i = 1:self.nfish
        if ~l(i)
          s = s+1;
          assert(d==size(X{s},2));
          self.n(i) = size(X{s},1);
          self.Sigma(:,:,i) = cov(X{s},1) + eye(d)*self.regularizer;
          self.mu(i,:) = mean(X{s},1);
        end
        self.updatePars(i);
      end
    end
    
    
    
        
    function [classprob] = predictWithoutPreProcessing(self,X)

      b = size(X,1);
      classprob = zeros(b,self.nfish);

      for j = 1:self.nfish
        if self.n(j)
          xm = bsxfun(@minus,X,self.mu(j,:));
          xms = xm*self.invSigma(:,:,j);
          classprob(:,j) = exp(-sum(xms.*xm,2)/2);%/self.norm(j);
        end
      end
      
      classprob = bsxfun(@rdivide,classprob,sum(classprob,2)+eps);
    end


    
    function updatePars(self,j)
    % update 
      self.norm(j) = sqrt(det(2*pi*self.Sigma(:,:,j)));
      self.invSigma(:,:,j) = inv(self.Sigma(:,:,j)+self.regularizer*eye(self.d));
    end
    
      
    function X = discardOutliers(self,X)
      % to avoid outliers train a few EM-cluster
      if self.outliersif
        [pc,a] = fish.helper.emclustering(X,3);
        [~,cl] = max(pc,[],2);
        idx = find(a<self.noveltyThres);
        if ~isempty(idx)
          delidx = ismember(cl,idx);
          X(delidx,:);
        end
      end
    end
    
  
  end
  
  methods 
    
    function reset(self,batchsample)
      self.Sigma = [];
      
    end

    
    function init(self,batchsample)
      self.initialize(batchsample)
    end
    
    function self = FishBatchClassifier(nfish,dim,varargin) % constructor
      
      self = self@handle();

      nargs = length(varargin);
      if nargs>0 && mod(nargs,2)
        error('expected arguments of the type ("PropName",pvalue)');
      end
      for i = 1:2:nargs
        if ~ischar(varargin{i})
          error('expected arguments of the type ("PropName",pvalue)');
        else
          self.(varargin{i}) = varargin{i+1};
        end
      end
      
      self.nfish = nfish;
      self.featdim = dim;
    end
    
    
    function X = preprocess(self,sample)
    % just simple pca for now
      if self.npca>0
        X = bsxfun(@minus,sample, self.pcamu)*self.pcapc; 
      else
        X = sample;
      end;
        
      if self.nlfd>0
        X = X*self.lfdpc; 
      end
      
    end
    
    
    function X = convertBatch(self,batchsample)
    
      if ~self.isInit()
        error('First initialize classifier!');
      end
      assert(iscell(batchsample));

      for s = 1:length(batchsample)

        X{s} = self.preprocess(batchsample{s});
        d = size(X{s},2);        

        if size(batchsample{s},1)>=self.minBatchN
          % to avoid outliers train a few EM-cluster
          X{s} = discardOutliers(self,X{s});
        end
      end

    end
    
    
        
    function [assignedFishIdx prob]= batchUpdate(self,fishidx,batchsample,varargin)
    % updates the fatures of the fishidx with the cell-array;
    % returns "nan" if not updated.
      
      if nargin>3
        forceif = varargin{1};
      else
        forceif = 0;
      end
      minBatchN = self.minBatchN;
      fishidx = fishidx(:)'; %row vector;
      
      assignedFishIdx = nan(size(fishidx));
      prob = nan(size(fishidx));
      l = cellfun(@(x)(size(x,1)),batchsample);
      if ~self.isInit() % not yet initialized
        if  all(l>=minBatchN) && length(unique(fishidx))==self.nfish
          self.initialize(batchsample);
          assignedFishIdx = fishidx;
        end
        return
      end
      

      validclasses = l(:)>=minBatchN;
      if ~any(validclasses) || isempty(fishidx)
        fish.helper.verbose('WARNING: no valid classes')
        return
      end
      
        
      X = self.convertBatch(batchsample(validclasses));
      if forceif 
        assignedFishIdx(validclasses) = fishidx(validclasses);
      else
        [assignedFishIdx(validclasses) prob(validclasses)] = self.predictPermutedAssignments(X,fishidx(validclasses),0);
        
        if any(assignedFishIdx(validclasses)~=fishidx(validclasses))
          fish.helper.verbose('WARNING: mixed up classes.')
          %return !?!?!?
        end
      end
      
      %update 
      s = 0;
      for i = 1:length(assignedFishIdx)
        
        if ~validclasses(i)
          continue;
        end

        s = s+1;
        if isnan(assignedFishIdx(s))
          continue;
        end
        idx = assignedFishIdx(s);
        
        mu_old = self.mu(idx,:);
        Sigma_old = self.Sigma(:,:,idx);
        b = size(X{s},1);
        n = min(self.n(idx),self.tau);
        
        % UPDATE WITH WEIGHT ACCORDING TO PROB WOULD BE BETTER (LIKE KALMAN GAIN)
        mu_new = n/(n+b)*mu_old + sum(X{s},1)/(b+n);
        self.Sigma(:,:,idx)  = n/(n+b)*(Sigma_old + mu_old'*mu_old) + X{s}'*X{s}/(n+b) - mu_new'*mu_new;
        self.mu(idx,:) = mu_new;
        self.n(idx) = n + b;
        
        
        self.updatePars(idx);
      end
    end


    function [assignedClassIdx prob] = predictPermutedAssignments(self,batchsample,assumedClassIdx,varargin)
      
      if nargin>3
        preprocessif = varargin{1};
      else
        preprocessif = 1;
      end
      assumedClassIdx = assumedClassIdx(:)'; % row vector
      
      if preprocessif
        X = self.convertBatch(batchsample);
      else
        X = batchsample;
      end
      
      
      cost = zeros(self.nfish,self.nfish);
      s = 0;
      validclasses = true(size(assumedClassIdx));
      for i = assumedClassIdx
        % check the distance of the new points to the
        % existing/newly created  clusters
        s = s+1;
        if ~isempty(X{s})
          cost(i,:) = mean(1-self.predictWithoutPreProcessing(X{s}),1); % HOW TO BEST COMBINE IT?.. mean should be OK
        else
          validclasses(s) = false;
        end
      end
      
      cost = cost(assumedClassIdx(validclasses),assumedClassIdx(validclasses))'; % only within
                                                    
      assignments =  fish.helper.assignDetectionsToTracks(cost,1e10);
      
      assignedClassIdx = nan(size(assumedClassIdx));
      assignedClassIdx(assignments(:,2)) = assumedClassIdx(assignments(:,1));
      prob = nan(size(assumedClassIdx));
      for i = 1:size(assignments,1)
        prob(assignments(i,2)) = 1-cost(assignments(i,1),assignments(i,2));
      end
      
    end
    
    function [classprob] = predictBatch(self,batchsamples)
      
      if ~iscell(batchsamples)
        classprob = mean(self.predict(batchsamples),1);
      else
      
        for i = 1:length(batchsamples)
          classprob{i} = mean(self.predict(batchsamples{i}),1);
        end
      end
    end
    
      
      
    function [classprob] = predict(self,sample)
    % sample needs to be n x dim. Calculates the class probabilities for
    % each sample according to a Gaussian mixture model (naive bayesian)
      if isempty(self.Sigma) % once updated.. classes  established 
        classprob = [];
      else
        classprob = self.predictWithoutPreProcessing(self.preprocess(sample));
      end
    end

    function mu = getMeans(self);
    % returns the means of the classes in "readable" format
      mu = zeros([self.nfish,self.featdim]);

      for i = 1:self.nfish
        if self.nlfd
          x = self.lfdpc*self.mu(i,:)';
        else
          x = self.mu(i,:)';
        end
        if self.npca
          x = bsxfun(@plus,self.pcapc*x,self.pcamu');
        end
        mu(i,:)  = x;
      end
      
    end
    
    function bool = isInit(self)
      bool = ~isempty(self.Sigma);
    end
    
      
    
    function plotMeans(self)
    % plots the means

      if ~self.isInit()
        return
      end
      
      clf;
      [r1,r2] = fish.helper.getsubplotnumber(self.nfish+(self.nfish>1));
      mu = self.getMeans();
      
      for i = 1:self.nfish
        a = subplot(r1,r2,i,'align');
      
        imagesc(squeeze(mu(i,:,:)));
        daspect([1,1,1])        
        axis off
        colorbar;

        title(i)
      end

      % difference of the last 2 fish
      if self.nfish>1
        a = subplot(r1,r2,i+1,'align'); 
        z = squeeze(mu(end-1,:,:) - mu(end,:,:));
        imagesc(z);
        daspect([1,1,1])        
        axis off
        colorbar;
        title('Difference')
      end      
      drawnow
    end
    
    
  end
end
