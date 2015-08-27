classdef FishForegroundDetectorCV <FishForegroundDetector;
  
  properties
    detector = [];
    detectShadows = 0;
    history = 500;
    expectedFrameFormat = 'U';
  end
  
    
  
  
  methods
  
    function a_init(self);
    %self.detector = cv.BackgroundSubtractorMOG2();
      self.detector = cv.BackgroundSubtractorKNN();
      self.detector.DetectShadows = self.detectShadows;
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
    
    
    function self = FishForegroundDetectorCV(varargin)

      self = self@FishForegroundDetector(varargin{:});
    
    end
  
  end
  
end
