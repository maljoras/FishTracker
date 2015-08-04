classdef FishForegroundDetector < handle;
  

  properties(Abstract)
    expectedFrameFormat;
    history;
  end
  
  properties;
    useScaledFormat = 0;
    inverse = 0;
  end
  
    
  methods(Abstract);
    
    bwmsk = a_step(self,frame);
    a_init(self);
    
  end
  

  
  methods
    
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
