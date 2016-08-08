classdef FishForegroundDetectorMatlabCV < xy.core.FishForegroundDetector;
  
  
  
  properties 
    
    plotif = 0;


    nskip = 5;  
    nAutoThres = 30;

  end

  properties (SetAccess = private);

    thres = 0.5; % wild guess. Will be adjusted 
    mframe = 0;
    frame = [];

    framecounter = 0;
    thresholder = [];
  end
  
  
  methods (Access=private)
    
    
    function bwimg = applyThres(self,frame)
      
    % bw image
      if self.inverted
        bwimg = (frame >= self.thres*(2-self.adjustThresScale));
      else
        bwimg = (frame <= self.thres*(self.adjustThresScale));
      end
    end

    function [bwimg, thres] = applyAutoThres(self,frame)
    % use the Autothresholder
      if self.inverted
        [bwimg, thres] = cv.threshold(frame, 'Otsu', 'Method', 'Binary');
      else
        [bwimg, thres] = cv.threshold(frame, 'Otsu', 'Method', 'BinaryInv');
      end
    end
    
  end
  
  
  methods (Access=protected)

    function  bwmsk = a_step(self,frame);

      self.framecounter =  self.framecounter  + 1;

      if self.history>0
        sframe = single(frame);
        frame1 = uint8(sframe - self.mframe + 127);
        self.updateMean(frame);
      else
        frame1 = frame;
      end


      if  self.framecounter>min(self.history,self.nAutoThres)
        bwmsk = self.applyThres(frame1);
      else
        [bwmsk,self.thres] = self.applyAutoThres(frame1);
      end
      
      
      if self.plotif
        figure(1);
        subplot(2,1,1)
        imagesc(frame1);
        title(self.framecounter)
        subplot(2,1,2)
        imagesc(bwmsk);
        drawnow;

      end
      
    end

    
    function a_init(self);
      self.framecounter = 0;
      self.thresholder  = vision.Autothresholder('ThresholdOutputPort',1);
      self.mframe = 0;
      self.expectedFrameFormat = 'U';
    end
    
    
    
  end
  
  
  
  methods    
    
    function self = FishForegroundDetectorMatlabCV(varargin)
      self = self@xy.core.FishForegroundDetector(varargin{:});
      % will call a_init
    end

    
    function  updateMean(self,frame)
      

      if  self.framecounter<self.history
        mtau = self.framecounter;
        self.mframe = (mtau-1)/mtau * self.mframe + 1/mtau*single(frame);
        
      elseif ~mod(self.framecounter,self.nskip)

        mtau = self.history/self.nskip;
        self.mframe = (mtau-1)/mtau * self.mframe + 1/mtau*single(frame); %
        
      else
        % pass
      end
    end
    
  end
end
