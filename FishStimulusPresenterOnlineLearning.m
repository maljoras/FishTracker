classdef FishStimulusPresenterOnlineLearning < FishStimulusPresenter;
  
  properties 
    switchInterval = 30; % switching the stimulus time (in seconds)
    stmInterval = 150; % switch between blank periods and stimulus presentations
    adaptationTime = 0; % time at the beginning (in seconds)
    col1 = [0,0,0]; % background color (RGB [0,0,0] for black)
    col2 = [1,1,1]; % stimulus color (RGB [1,1,1] for white)
    midline = 0.5;  % position of the line form 0..1
    stmidx = -1;
    
    flashCol = [1,1,1];
    flashSize = 100; % in pix
    flashBorder = 0.05; % in percent
  end
  
  methods 

    function borderFlash(self,x,y)
      
      for i = 1:length(x)
        if x(i)>1-self.flashBorder
          self.plotDot(1,y(i),self.flashSize,self.flashCol);
        elseif x(i)<self.flashBorder
          self.plotDot(0,y(i),self.flashSize,self.flashCol);
        elseif y(i)<self.flashBorder
          self.plotDot(x(i),0,self.flashSize,self.flashCol);
        elseif y(i)>1-self.flashBorder
          self.plotDot(x(i),1,self.flashSize,self.flashCol);
        end
      end
      
    end
    
    function tracks = step(self,tracks,framesize,t)
    % this function will be called from FishTracker after each round
    
      if isempty(tracks)
        return;
      end
      
      fishId = [tracks.fishId];

      
      if isempty(self.screenBoundingBox)
        sbbox = [1,1,framesize([2,1])-1];
      else
        sbbox = self.screenBoundingBox;
      end
      
      for i =1:length(tracks)

        x(i) = (tracks(i).centroid(1)-sbbox(1))/sbbox(3);
        y(i) = (tracks(i).centroid(2)-sbbox(2))/sbbox(4);
        col(i,:) = [1-(fishId(i)==1),(fishId(i)==2),1 - (fishId(i)-1)/length(tracks)];
      end
      self.borderFlash(x,y);

      idx = x<0.5;
      if any(idx)
        self.plotDot(x(idx),y(idx),50,col(idx,:));
      end

      % also save the flashs ??
      for i =1:length(tracks)
        if idx(i)
          tracks(i).stmInfo = [x(i),y(i)];
        else
          tracks(i).stmInfo = [NaN,NaN];          
        end
      end
      self.flip();
    end
    
    
    
    
  end
end


