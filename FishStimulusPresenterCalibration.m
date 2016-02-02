classdef FishStimulusPresenterCalibration < FishStimulusPresenter;
  
  properties 
    freq = 0.1;
    col = [1,1,1]; % stimulus color (RGB [1,1,1] for white)
    width = 100;
    onset = 10;
  end
  
  methods 

    
    function tracks = step(self,tracks,framesize,t)
    % this function will be called from FishTracker after each round

      if t<self.onset
        return;
      end
      
      pos = sin(2*pi*t*self.freq)/2+1;
      self.patch(1-pos,0,self.width,self.width,self.col);
      self.patch(pos,1,self.width,self.width,self.col);
      self.patch(1,1-pos,self.width,self.width,self.col);
      self.patch(0,pos,self.width,self.width,self.col);
    
      self.flip();
      
    end
    
  end
end
