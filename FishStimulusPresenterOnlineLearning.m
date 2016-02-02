classdef FishStimulusPresenterPlane < FishStimulusPresenter;
  
  properties 
    switchInterval = 30; % switching the stimulus time (in seconds)
    stmInterval = 150; % switch between blank periods and stimulus presentations
    adaptationTime = 0; % time at the beginning (in seconds)
    col1 = [0,0,0]; % background color (RGB [0,0,0] for black)
    col2 = [1,1,1]; % stimulus color (RGB [1,1,1] for white)
    midline = 0.5;  % position of the line form 0..1
    stmidx = -1;
  end
  
  methods 

    function borderFlash(self,tracks,framesize,t)
      
      
      
      
    end
    
    
    
    
    function tracks = step(self,tracks,framesize,t)
    % this function will be called from FishTracker after each round
    
      oldstmidx = self.stmidx;

      if t> self.adaptationTime
        deltat  = t - self.adaptationTime;
        if ~mod(floor(deltat/self.stmInterval),2)
          
          deltat2 = mod(deltat,self.stmInterval);
          if mod(floor(deltat2/self.switchInterval),2)
            self.stmidx = 1;
          else
            self.stmidx = 2;
          end
        else
          self.stmidx = 3;
        end
      else
        self.stmidx = 0;
      end
      
      if oldstmidx~=self.stmidx
        verbose('Switch Stimulus Plane %d -> %d',oldstmidx,self.stmidx);
        switch self.stmidx
          case 1
            self.plotVPlane(self.midline,self.col1,self.col2);
          case 2 
            self.plotVPlane(self.midline,self.col2,self.col1);
          case {3,0}
            self.plotVPlane(self.midline,self.col1,self.col1);
        end
        self.flip();      
      end
      
      if ~isempty(tracks)
        [tracks.stmInfo] = deal([t, self.stmidx]);
      end
    end
    
  end
end


