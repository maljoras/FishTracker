classdef BatchClassifier < handle;
% Classifies the  

  properties 
    npca = 15;
    nlfd = 0;
    nbody = 2;
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

        obj = xy.core.BatchClassifier(S.nbody,S.featdim);
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

      
      if ~iscell(batchsample) || length(batchsample)~=self.nbody || ...
            all(cellfun('isempty',batchsample))
        xy.helper.verbose('WARNING: provide a batch for each fish!');
        keyboard
        return
      end
      
      
      if prod(self.featdim)==size(batchsample{1},2)
      elseif 3*prod(self.featdim)==size(batchsample{1},2)
        self.featdim = [self.featdim ,3];        
      else
        xy.helper.verbose('WARNING: Set feature dimensions correctly, otherwise plotting will be impiared');
        self.featdim = [size(batchsample{1},2),1];        
      end
      
      
      xy.helper.verbose('Initialize BatchClassifier...');
      
      if self.npca>0 || self.nlfd>0
        Z = cat(1,batchsample{:});
        % do some weak oulyier detection
        [pc,a] = xy.helper.emclustering(Z,self.nbody+2);
        [~,cl] = max(pc,[],2);
        idx = find(a<self.noveltyThres/self.nbody);
        goodidx = ~ismember(cl,idx);
        
        cl = cellfun(@(x,k)k*ones(size(x,1),1),batchsample,num2cell(1:length(batchsample)),'uni',0);
        cl = cat(1,cl{:});
        [self.pcapc,~,proj,self.pcamu] = xy.helper.pca1(Z(goodidx,:),max(self.nlfd,self.npca));
        self.npca = max(self.nlfd,self.npca);
        if size(Z,2) < self.nlfd % too many very small fish 
          self.nlfd = 0; % do not do any lfd
        end
        
        if self.nlfd>0 && self.npca>self.nlfd 
          [self.lfdpc,~,proj] = xy.helper.lfd(proj,cl(goodidx),self.nlfd,0);
        end
          
        self.pcamu = self.pcamu';
        d = size(proj,2);
      else
        d = size(batchsample{1},2);
      end

      % initialize 
      self.d = d;
      self.n = zeros(self.nbody,1);
      self.mu = zeros(self.nbody,d);
      self.Sigma = repmat(eye(d),[1,1,self.nbody]);
      self.invSigma = repmat(eye(d),[1,1,self.nbody]); % dummy
      self.norm = ones(self.nbody,1);
      
      % convert the batchsample
      l = cellfun('isempty',batchsample);
      X = self.convertBatch(batchsample(~l));

      % initial setting
      s = 0;
      for i = 1:self.nbody
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
      classprob = zeros(b,self.nbody);

      for j = 1:self.nbody
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
        [pc,a] = xy.helper.emclustering(X,3);
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
    
    function self = BatchClassifier(nbody,dim,varargin) % constructor
      
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
      
      self.nbody = nbody;
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
    
    
        
    function [assignedIdentityIdx prob]= batchUpdate(self,fishidx,batchsample,varargin)
    % updates the fatures of the fishidx with the cell-array;
    % returns "nan" if not updated.
      
      if nargin>3
        forceif = varargin{1};
      else
        forceif = 0;
      end
      minBatchN = self.minBatchN;
      fishidx = fishidx(:)'; %row vector;
      
      assignedIdentityIdx = nan(size(fishidx));
      prob = nan(size(fishidx));
      l = cellfun(@(x)(size(x,1)),batchsample);
      if ~self.isInit() % not yet initialized
        if  all(l>=minBatchN) && length(unique(fishidx))==self.nbody
          self.initialize(batchsample);
          assignedIdentityIdx = fishidx;
        end
        return
      end
      

      validclasses = l(:)>=minBatchN;
      if ~any(validclasses) || isempty(fishidx)
        xy.helper.verbose('WARNING: no valid classes')
        return
      end
      
        
      X = self.convertBatch(batchsample(validclasses));
      if forceif 
        assignedIdentityIdx(validclasses) = fishidx(validclasses);
      else
        [assignedIdentityIdx(validclasses) prob(validclasses)] = self.predictPermutedAssignments(X,fishidx(validclasses),0);
        
        if any(assignedIdentityIdx(validclasses)~=fishidx(validclasses))
          xy.helper.verbose('WARNING: mixed up classes.')
          %return !?!?!?
        end
      end
      
      %update 
      s = 0;
      for i = 1:length(assignedIdentityIdx)
        
        if ~validclasses(i)
          continue;
        end

        s = s+1;
        if isnan(assignedIdentityIdx(s))
          continue;
        end
        idx = assignedIdentityIdx(s);
        
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
      
      
      cost = zeros(self.nbody,self.nbody);
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
                                                    
      assignments =  xy.helper.assignDetectionsToTracks(cost,1e10);
      
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
      mu = zeros([self.nbody,self.featdim]);

      for i = 1:self.nbody
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
      [r1,r2] = xy.helper.getsubplotnumber(self.nbody+(self.nbody>1));
      mu = self.getMeans();
      
      for i = 1:self.nbody
        a = subplot(r1,r2,i,'align');
      
        imagesc(squeeze(mu(i,:,:)));
        daspect([1,1,1])        
        axis off
        colorbar;

        title(i)
      end

      % difference of the last 2 fish
      if self.nbody>1
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
