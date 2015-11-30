classdef FishForegroundDetectorMatlab < FishForegroundDetector;
  
 
    
  properties
    
    plotif = 0;
    meanSkip = 10;  
    nAutoThres = 30;

    subtractTotalMean = 1; % ALREADY DONE IN C
  end

  properties (SetAccess = private);

    thres = 0.5; % wild guess. Will be adjusted 
    mframe = 0;
    frame = [];
    
    frameMean = 0;
    framecounter = 0;
    thresholder = [];
  end
  
  
  methods (Access=private)
    

    function bwimg = applyThres(self,frame)
      
      % bw image
      if self.inverted
        bwimg = frame >= (self.thres + self.frameMean);
      else
        bwimg = frame <= (self.thres+ self.frameMean);
      end
    end

    function [bwimg, thres] = applyAutoThres(self,frame)
    % use the Autothresholder
      [bwimg thres] =  step(self.thresholder, frame);
    end
    

    

    function  updateMean(self,frame)
      
      if  self.framecounter<self.history
        mtau = self.framecounter;
        self.mframe = (mtau-1)/mtau * self.mframe + 1/mtau*frame;
      elseif ~mod(self.framecounter,self.meanSkip)
        mtau = self.history/self.meanSkip;
        self.mframe = (mtau-1)/mtau * self.mframe + 1/mtau*frame; 
      end
      
    end
  end     

  

  methods(Access=protected)
      
    function a_setHistory(self,value);
    end
    
    function  bwmsk = a_step(self,frame);

      self.framecounter =  self.framecounter  + 1;
      
      if self.history>0
        frame1 = frame - self.mframe;
        self.updateMean(frame);
      end

      if self.subtractTotalMean
        % to exclude general light effects
        self.frameMean =  mean2(frame1);
      else
        self.frameMean = 0;
      end
      

      if  self.framecounter>min(self.history,self.nAutoThres)
        bwmsk = self.applyThres(frame1);
      else
        [bwmsk,self.thres] = self.applyAutoThres(-frame1);
        self.thres = -self.thres;
      end

      self.frame = frame1;
      
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

    function a_reset(self)
      self.framecounter =0;
      self.mframe = 0;
    end

    
    function a_init(self);
      self.expectedFrameFormat = 'S';
      self.framecounter = 0;
      self.thresholder  = vision.Autothresholder('ThresholdOutputPort',1);
      self.mframe = 0;
    end
  end
  
  methods    
    
    function self = FishForegroundDetectorMatlab(varargin)
      self = self@FishForegroundDetector(varargin{:});
      % will call a_init
    end
    
  end
end
