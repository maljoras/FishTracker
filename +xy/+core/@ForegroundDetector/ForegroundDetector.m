classdef ForegroundDetector < handle;
  

  
  properties
    expectedFrameFormat = 'U';

    useScaledFormat = 0;
    inverted = 0;  
    history = 500;
    
    adjustThresScale = 1;
  end


  
  properties(Access=private)
    detector = [];
    detectShadows = 0;
    bwmsk;
  end
  
    
  
  methods(Access=protected)
  
    function a_init(self);
    %self.detector = cv.BackgroundSubtractorMOG2();
      self.detector = [];
      self.detector = cv.BackgroundSubtractorKNN();
      self.detector.DetectShadows = self.detectShadows;
      self.detector.History = self.history;
    end
    
    function a_setHistory(self,value)
      self.detector.History = self.history;
    end

    function a_reset(self)
      self.a_init();
    end
    
    
    function  bwmsk = a_step(self,frame);
      if self.inverted
        if isa(frame,'uint8')
          bwmsk = self.detector.apply(uint8(255)-frame);
        else
          bwmsk = self.detector.apply(1-frame);
        end
      else
        bwmsk = self.detector.apply(frame);
      end
    end
  end
  
  
  methods
   

    function set.history(self,value)
      self.a_setHistory(value);
      self.history = value;
    end
    
    function self = ForegroundDetector(varargin)
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
    
    function reset(self)
      self.a_reset();
    end

    
    function bwmsk = step(self,frame);
      bwmsk = self.a_step(frame);
      self.bwmsk = bwmsk;
    end

    function bwmsk = getCurrentBWImg(self);
      bwmsk =  self.bwmsk;
    end

    
  end
  
  
  
end
