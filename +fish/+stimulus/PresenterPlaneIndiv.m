classdef PresenterPlaneIndiv < fish.stimulus.PresenterOnlineLearningCue

   
  properties
    nRoundsPerGroup = 4; % switch betweem left & right stimulation
    nPauseStmGroupsPerExp = 2;  %   
    nFishStimPerExps = [1,5,3,4,5];
    
  end
  
  
  
  methods  
    
    function init(self,varargin)

      
      init@fish.stimulus.PresenterOnlineLearningCue(self,varargin{:});

      % fix some values
      self.gapTime =0; % no gap      
      self.signalTime = 0;
      self.funRoundStmPause  = @(iround) ~mod(floor((iround-1)/self.nRoundsPerGroup),2); 
      %self.stmBkgType = 'plane';
      self.lrif = 1;
      self.lrSwitchProb =  1;
      self.stmLambda = 0;
      
      self.updateRoundTime();

    end
    
    function bool = isFinished(self,t)
      iExp = self.getExpIdx(t);
      bool = iExp>length(self.nFishStimPerExps);
    end
    
    function iExp = getExpIdx(self,t);
      tt = t - self.adaptationTime;
      if tt<0
        iExp =0;
      else
        iround = floor(tt/self.roundTime)+1;
        nRoundsPerExp = self.nPauseStmGroupsPerExp*self.nRoundsPerGroup*2;
        iExp = floor((iround-1)/nRoundsPerExp) + 1;
      end
    end
    
    
    function stmInfo = stepStimulus(self,x,y,t,fishIds)

      iExp = self.getExpIdx(t);
      if iExp>0
        nFishStm = self.nFishStimPerExps(iExp);
      else
        nFishStm = 0;
      end
      
      self.funLRStm = @(fishIds) fishIds<=nFishStm;
      
      % let the super do the work
      stmInfo1= stepStimulus@fish.stimulus.PresenterOnlineLearningCue(self,x,y,t,fishIds);
      stmInfo  = cat(2,stmInfo1, zeros(length(fishIds),3));
      stmInfo(:,size(stmInfo1,2)+1) = self.iround;
      stmInfo(:,size(stmInfo1,2)+2) = iExp;
      stmInfo(:,size(stmInfo1,2)+3) = nFishStm;
    end
    
    
  end
  


end