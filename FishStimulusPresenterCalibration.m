classdef FishStimulusPresenterCalibration < FishStimulusPresenter;
  
  properties 
    freq = 0.1;
    col = [1,1,1]; % stimulus color (RGB [1,1,1] for white)
    width = 100;
  end
  
  methods 
    
    function self = FishStimulusPresenterCalibration(varargin)
    
      self = self@FishStimulusPresenter(varargin{:});
      self.setOpts(varargin{:});
    end
    
      
    function tracks = step(self,tracks,framesize,t)
    % this function will be called from FishTracker after each round

      w = self.width/framesize(1);
      h = self.width/framesize(2);
      
      pos = (sin(2*pi*t*self.freq)+1)/2;
      self.patch((1-pos)*(1-h),0,h,w,self.col);
      self.patch(pos*(1-h),1-w,h,w,self.col);
      self.patch(1-h,(1-pos)*(1-w),h,w,self.col);
      self.patch(0,pos*(1-w),h,w,self.col);


      self.flip();

    end
    
  end
end
