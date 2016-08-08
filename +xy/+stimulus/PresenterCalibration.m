classdef PresenterCalibration < xy.stimulus.Presenter;
  
  properties 
    freq = 0.1;
    col = [1,1,1]; % stimulus color (RGB [1,1,1] for white)
    width = 50;
  end
  
  methods 
    
    function self = PresenterCalibration(varargin)
    
      self = self@xy.stimulus.Presenter(varargin{:});
      self.setOpts(varargin{:});

    end
    
      
    function tracks = step(self,tracks,framesize,t)
    % this function will be called from FishTracker after each round

      w = self.width/self.windowSize(1);
      h = self.width/self.windowSize(2);
      pos = (sin(2*pi*t*self.freq)+1)/2;
      self.patch((1-pos)*(1-h),0,h,w,self.col);
      self.patch(pos*(1-h),1-w,h,w,self.col);
      self.patch(1-h,(1-pos)*(1-w),h,w,self.col);
      self.patch(0,pos*(1-w),h,w,self.col);

      self.flip();

      [tracks(:).stmInfo] = deal(NaN);
      
    end
    
  end
end
