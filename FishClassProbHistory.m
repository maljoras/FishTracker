classdef FishClassProbHistory < handle;
  
  
  
  properties 
    nHistory = 2000;
    lambda = 1.5;
    taulambda = 100; % for lambda adjustements. Set to Inf if not wanted
    reasonableThres =  0.1353/2; % 
  end
  
  
  properties(SetAccess=private);

    buffer = [];
    weightbuffer = [];
    currentIdx = 0;
    age =0 ;
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
      self.currentIdx = 0;
      self.age = 0;
    end

    
    
    function [reasonable w] = update(self,vector,somenoiseval)
      
      self.age = self.age +1;
      self.currentIdx = self.currentIdx + 1;

      if self.currentIdx>self.nHistory
        self.currentIdx = 1;
      end
      % store vector
      self.buffer(self.currentIdx,:) = vector;

      w = exp(-somenoiseval/self.lambda);

      if isnan(w*sum(vector))
        w = 0;
      end
      % store w
      self.weightbuffer(self.currentIdx,1) = w;
    
      % also update lambda 
      if w && ~mod(self.age,5)
        p =  1/self.taulambda;
        self.lambda = self.lambda*(1-p) + p*somenoiseval;
      end

      reasonable = w>self.reasonableThres;      
    end

    function bool = currentNoiseReasonable(self)
    % checks whether the noise of the current sample exceeds a thresold
      bool = self.weightbuffer(self.currentIdx) >  self.reasonableThres; % 2 x lambda
    end
    
    
    function reset(self)
      self.currentIdx = 0;
      self.age = 0;
    end
    
    
    function [classProbs, weights] = getData(self,N)
    % returns available data. The order is from past to preset

      t = self.currentIdx;
      if (N>=self.nHistory) || N>self.age
        error(['Too long history requested. Consider incresing ' ...
               'nHistory']);
      end
      N = min(mod(N-1,self.nHistory)+1,self.age); 
      idx = mod(t-N:t-1,self.nHistory)+1;

      weights = self.weightbuffer(idx,1);
      classProbs = self.buffer(idx,:);
      
    end

      
    function [mclassProb usedN] = mean(self,N)
    % MCLASSPROB = MEAN(SELF,N) returns the weighted mean
    % class prob for the last N steps 
      
      [cl,w] = self.getData(N);

      msk = w>self.reasonableThres;
      w = w(msk);
      
      if isempty(w)
        mclassProb = NaN([1,size(self.buffer,2)]);
        usedN = 0;
      else
        w = w/sum(w);
        mclassProb = sum(bsxfun(@times,cl(msk,:),w),1);
        usedN = length(w);
      end
      
    end
    
      
    
  end
  
    
  
end
