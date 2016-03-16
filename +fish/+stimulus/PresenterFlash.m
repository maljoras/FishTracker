classdef PresenterFlash < fish.stimulus.Presenter;
  
  properties 

    flashCol = [1,1,1];
    flashSize = 100; % in pix
    flashBorder = 0.05; % in percent
    flashLambda = 0.01;
  
  end
  
  properties(SetAccess=private, GetAccess=private)
    flashState = 1;
  end
  
    
  methods 

    function borderFlash(self,x,y)

      %% change stimulus state with poisson characeteristics
      if self.flashLambda>0
        msk = self.flashLambda>rand(size(self.flashState))
        self.flashState(msk) = ~self.flashState(msk);
      end

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
    
    
    function stmInfo = stepStimulus(self,x,y,t,fishIds)
    % this function will be called from fish.Tracker after each round
    
      self.borderFlash(x,y);

    end
    
       
    
  end
end


