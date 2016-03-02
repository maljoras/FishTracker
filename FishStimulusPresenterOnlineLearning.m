classdef FishStimulusPresenterOnlineLearning < FishStimulusPresenter;
  
  properties 
    stmSize = 50;
    flashCol = [1,1,1];
    flashSize = 100; % in pix
    flashBorder = 0.05; % in percent
    midLine = 0.6
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
    
    function bool = isFinished(self,t);
      bool = t>self.tmax;
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

      idx = x<self.midLine;
      if any(idx)
        self.plotDot(x(idx),y(idx),self.stmSize,col(idx,:));
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


