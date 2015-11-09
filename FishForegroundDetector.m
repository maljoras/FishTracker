classdef FishForegroundDetector < handle;
  

  
  properties
    expectedFrameFormat = 'U';

    useScaledFormat = 0;
    inverse = 0;  
    history = 500;
  
  end

  
  
  properties(Access=private)
    detector = [];
    detectShadows = 0;
  end
  
    
  
  methods(Access=protected)
  
    function a_init(self);
    %self.detector = cv.BackgroundSubtractorMOG2();
      self.detector = cv.BackgroundSubtractorKNN();
      self.detector.DetectShadows = self.detectShadows;
      self.detector.History = self.history;
    end
    
    function a_setHistory(self,value)
      self.detector.History = self.history;
    end
    
    function  bwmsk = a_step(self,frame);
      bwmsk = self.detector.apply(frame);
% $$$       if self.inverse
% $$$         bwmsk = ~cv.adaptiveThreshold(frame, 1, 'AdaptiveMethod', 'Mean','ThresholdType','BinaryInv');
% $$$       else
% $$$         bwmsk = ~cv.adaptiveThreshold(frame, 1, 'AdaptiveMethod', 'Mean','ThresholdType','Binary');
% $$$       end
    end
  end
  
  
  methods
    
    
    function set.history(self,value)
      self.a_setHistory(value);
      self.history = value;
    end
    
    function self = FishForegroundDetector(varargin)
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
    
      self.a_init();
    end
    
    function bwmsk = step(self,frame);
      bwmsk = self.a_step(frame);
    end
  
  end
  
  
  
end
