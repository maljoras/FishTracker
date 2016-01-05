classdef FishClassProbHistory < handle;
  
  
  
  properties 
    nHistory = 500;
    lambda = 1.5;
    tau = 10; % for moving average (step==frame)
    taulambda = 100; % for lambda adjustements. Set to Inf if not wanted
    reasonableThres =  0.1353/2; % 
  end
  
  
  properties(SetAccess=private);

    buffer = [];
    weightbuffer = [];
    currentIdx = 0;
    age =0 ;
    movage = 0;
    movavg = [];
    movavgLast = [];
  end
  
  
  
  
  methods
    
    function self = FishClassProbHistory(nfish, varargin) 
    % constructor
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

      self.buffer = nan(self.nHistory,nfish);
      self.weightbuffer = nan(self.nHistory,1);
      self.movavg = zeros(1,nfish);
      self.movavgLast = zeros(1,nfish);
      self.currentIdx = 0;
      self.age = 0;
      self.movage = 0;
    end

    
    
    function update(self,vector,somenoiseval)
      self.currentIdx = mod(self.currentIdx,self.nHistory) + 1;
      self.age = self.age +1;
      self.movage = self.movage +1;
      self.buffer(self.currentIdx,:) = vector(:)';

      if isnan(somenoiseval) || any(isnan(vector))
        self.weightbuffer(self.currentIdx,1) = 0;
      else
        self.weightbuffer(self.currentIdx,1) = exp(-somenoiseval/self.lambda);
      end
      
      if self.movage >= self.tau
        p =  self.weightbuffer(self.currentIdx)/self.tau;
      else
        p =  1/min(self.tau,self.movage);
      end
      
      self.movavgLast = self.movavg;
      self.movavg = nansum([self.movavg*(1-p); p*self.buffer(self.currentIdx,:)]);
    
      % also update lambda 
      if ~isinf(self.taulambda)
        p =  1/self.taulambda;
        self.lambda = nansum([self.lambda*(1-p),p*somenoiseval]);
      end
      
    end

    function bool = currentNoiseReasonable(self)
    % checks whether the noise of the current sample exceeds a thresold
      bool = self.weightbuffer(self.currentIdx) >  self.reasonableThres; % 2 x lambda
    end
    
    
    function reset(self)
      self.movage = 0;
      self.movavgLast(:) = NaN;
      self.movavg(:) = NaN;
    
      % also reset the other stuff ? 
    end
    
    
    function [classProbs, weights] = getData(self,N)
    % returns available data. The order is from past to preset

      t = self.currentIdx;
      N = min(mod(N-1,self.nHistory)+1,self.age); 
      idx = mod(t-N:t-1,self.nHistory)+1;

      weights = self.weightbuffer(idx,1);
      classProbs = self.buffer(idx,:);
      
    end

    function mclassProb = leakyMaxProb(self,omitLatestIf)
      if nargin>1 && omitLatestIf
        mclassProb = max(self.movavgLast);
      else
        mclassProb = max(self.movavg);
      end
    end
    
    function mclassProb = leakyMean(self,omitLatestIf)
      if nargin>1 && omitLatestIf
        mclassProb = self.movavgLast;
      else
        mclassProb = self.movavg;
      end
    end
    
    function id = leakyID(self)
      mClassProb = self.leakyMean();
      [~,id] = max(mClassProb);
    end
      
    function [mclassProb usedN] = mean(self,N,omitLatestIf)
    % MCLASSPROB = MEAN(SELF,N,OMITLATESTIF) returns the weighted mean
    % class prob for the last N steps (omitting the latest prob if
    % OMITLATESTIF is set, default: OMITLATESTIF==0).
      
      if nargin<3
        omitLatestIf = 0;
      end
      
      [cl,w] = self.getData(N);
      if omitLatestIf
        cl = cl(1:end-1,:);
        w = w(1:end-1,:);
      end
      msk = w>self.reasonableThres;
      w = w(msk);
      
      if isempty(w)
        mclassProb = NaN*self.movavg;
        usedN = 0;
      else
        w = w/sum(w);
        mclassProb = sum(bsxfun(@times,cl(msk,:),w),1);
        usedN = length(w);
      end
      
    end
    
      
    
  end
  
    
  
end
