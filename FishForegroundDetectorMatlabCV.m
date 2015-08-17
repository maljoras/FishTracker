classdef FishForegroundDetectorCV < FishForegroundDetector;
  
 
    
  properties 
    
    history = 100; % former mtau
    inverse = 0;

    plotif = 0;
    expectedFrameFormat = 'S';

    meanSkip = 5;  
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
      if self.inverse
        bwimg = (frame >= self.thres);
      else
        bwimg = (frame <= self.thres);
      end
    end

    function [bwimg, thres] = applyAutoThres(self,frame)
    % use the Autothresholder
      if self.inverse
        [bwimg, thresh] = cv.threshold(-frame, 'auto', 'Method', 'Binary');
      else
        [bwimg, thresh] = cv.threshold(-frame, 'auto', 'Method', 'BinaryInv');
      end
    end
    
  end
  

  methods    
    
    function self = FishForegroundDetectorCV(varargin)
      self = self@FishForegroundDetector(varargin{:});
      % will call a_init
    end

    
    function  bwmsk = a_step(self,frame);

      self.framecounter =  self.framecounter  + 1;

      if self.history>0
        frame1 = frame - self.mframe;
        self.updateMean(frame);
      else
        frame1 = frame;
      end

      % to exclude general light effects
      frame1 = frame1 - mean2(frame1);

      if  self.framecounter>min(self.history,self.nAutoThres)
        bwmsk = self.applyThres(frame1);
      else
        [bwmsk,self.thres] = self.applyAutoThres(-frame1);
        self.thres = -self.thres;
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
    end
    
    
     function  updateMean(self,frame)
       

       if  self.framecounter<self.history
         mtau = self.framecounter;
         self.mframe = (mtau-1)/mtau * self.mframe + 1/mtau*frame;
       
       elseif ~mod(self.framecounter,self.meanSkip)

         mtau = self.history/self.meanSkip;
         self.mframe = (mtau-1)/mtau * self.mframe + 1/mtau*frame; %
       
       else
         % pass
       end
     end
     
  end
end
