classdef PresenterDot < fish.stimulus.Presenter;
  
  properties
    stmTime = 10; %time of stimulus (in seconds)
    gapTime = 50; % time of gap after stim
    adaptationTime =300;%600; % time at the beginning (in seconds
    
    colBackground = [0,0,0]; % background color (RGB [0,0,0] for black)
    
    funStmFishIdProb = @(fishIds) (fishIds==1);

    stmLambda = 1; % this eq to dt*PosissonRate  !
    stmSize =100;
    stmCol = [1,1,1]; % stimulus color (RGB [1,1,1] for white)    

    nRound = 100;
    
  end

  properties(SetAccess=private)
    ID_ADAPTATION = -1;
    ID_GAP = 0;
    ID_STIMULATION = 1;
    
    roundTime = 0;
    iround = 0;
  end
  
  properties(SetAccess=private,GetAccess=private)
    testingif = 0;

    stmIdx = -1;
    stmState = [];
    textures = [];
    stimulatedMsk = [];
  end
    

  
  methods
    
    function init(self,varargin)
      
      init@fish.stimulus.Presenter(self,varargin{:});
      
      updateRoundTime(self);
    end

    
   function [x,y] = plotStimulus(self,x,y,t,fishIds,lrswitch)

     if length(x) ~= length(self.stmState)
       % should only happen once in the beginning. (nfish is constant)
       self.stmState = ones(size(x))';
     end

     % border
     %b = self.stmSize/self.windowSize(1)/2;
     %x = max(x,b);
     %x = min(x,1-b);
     
     
     %% change stimulus state with poisson characeteristics
     if self.stmLambda>0
       msk = self.stmLambda>rand(size(self.stmState));
       self.stmState(msk) = ~self.stmState(msk);
     end

     
     msk = self.stmState(fishIds) & self.stimulatedMsk;
     self.plotDot(x(msk),y(msk),self.stmSize,self.stmCol);


     x(~msk) = NaN;
     y(~msk) = NaN;
   end
   
   function updateRoundTime(self);
     trainingTime = self.stmTime;
     self.roundTime = trainingTime + self.gapTime;
   end
   
   function stmInfo = stepStimulus(self,x,y,t,fishIds)
   % stmInfo = stepStimulus(self,x,y,t,fishIds) this function will be
   % called from fish.stimulus.Presenter/step once for each frame. 
   %  
   % STMINFO: [STMIDX X Y T ]  
   %   
   % X,Y is position of stimulus or NaN if not given. 
   %
   % STMIDX is 
   %       -1 :  Adaptation time
   %       0 :  Gap period
   %       1 :  Stimulation period (+ 8 LR switch)

            
     %% set the stmidx
     oldstmIdx = self.stmIdx;
     lastRound = self.iround;

    
     if isempty(self.stimulatedMsk)
       self.stimulatedMsk = zeros(length(fishIds),1);
     end
     
     if t < self.adaptationTime
       self.stmIdx = self.ID_ADAPTATION;
     else

       % within round
       tt = t - self.adaptationTime;
       ttmod = mod(tt,self.roundTime);
       
       self.iround = floor(tt/self.roundTime)+1;
       
       if ttmod<self.stmTime
         % stimulus time
         if self.iround~=lastRound
           % new stim round
           self.stimulatedMsk = rand(length(fishIds),1) < self.funStmFishIdProb(fishIds(:)));
         end
         self.stmIdx = self.ID_STIMULATION;
       
       else
         self.stmIdx = self.ID_GAP;
       end
       
     
     end

     % background
     self.plotVPlane(0.5,self.colBackground,self.colBackground);
           

     %% plotting
     if self.stmIdx==self.ID_STIMULATION
       [x,y] = self.plotStimulus(x,y,t,fishIds,self.stimulatedMsk); 
     end

     % save stm info
     stmInfo = zeros(length(fishIds),5);
     stmInfo(:,1) = self.stmIdx;
     stmInfo(:,2:3) = [x,y];
     stmInfo(:,4) = t;
     stmInfo(:,5) = self.stimulatedMsk;
   end
   
        
   function bool = isFinished(self,t)
     bool = t>self.adaptationTime+(self.stmTime+self.gapTime)*self.nRound;
   end
   
  end
  
end



    
    