classdef FishStimulusPresenterOnlineLearningCue < FishStimulusPresenter;
  
  properties
    signalTime = 2;  %time of begining and ending signal (in seconds)
    stmTime = 40; %time of stimulus (in seconds)
    gapTime =40; % time of testing (in seconds)
    adaptationTime =600; % time at the beginning (in seconds
    
    switchInterval = 30; % switching the stimulus time (in seconds)
    stmInterval = 150; % switch between blank periods and stimulus presentations
    
    colBackgound = [0,0,0]; % background color (RGB [0,0,0] for black)
    colStimulus = [1,1,1]; % stimulus color (RGB [1,1,1] for white)
    colBeginCue = [255,0,0];% begining signal (RGB[255,0,0]for red)
    colEndCue = [0,255,0]; % ending signal(RGB[0,255,0]for green)

    colSideMarkerLeft = [122,122,0];
    colSideMarkerRight = [0,122,122];
    sideMarkerPos = 0.0;
    sideMarkerWidth = 5; % in pix
    
    midline = 0.5;  % position of the line form 0..1
    stmIdx = -1;
    
    testingProb = 0.1;
    stmSize =150;
    nRound= 100;
  
  
  end

  properties(SetAccess=private)
    ID_ADAPTATION = 0;
    ID_BEGINCUE = 1;
    ID_ENDCUE = 2;
    ID_STIMULUS = 3;
    ID_GAP = 4;
    ID_TEST = 5;
  end
  
  properties(SetAccess=private,GetAccess=private)
    testingif = 0;
  end
    

  
  methods
    
   function [x,y] = plotStimulus(self,tracks,framesize,t)


     if isempty(tracks)
       return;
     end
     
     fishIds = cat(1,tracks.fishId);
     
     if isempty(self.screenBoundingBox)
       sbbox = [1,1,framesize([2,1])-1];
     else
       sbbox = self.screenBoundingBox;
     end

     centroids = cat(1,tracks.centroid);
     x = (centroids(:,1)-sbbox(1))/sbbox(3);
     y = (centroids(:,2)-sbbox(2))/sbbox(4);
     
     msk = mod(fishIds,2);
     xmsk = x<self.midline;

     left = xmsk & msk;
     self.plotDot(x(left),y(left),self.stmSize,self.colStimulus);

     right = ~xmsk & ~msk;
     self.plotDot(x(right),y(right),self.stmSize,self.colStimulus);

     nostmmsk = ~(right | left);
     x(nostmmsk) = NaN;
     y(nostmmsk) = NaN;
   end
   
   function tracks = step(self,tracks,framesize,t)
   % this function will be called from FishTracker after each
   % round. It plots an stimulus below the fish if fishID odd and on left
   % side and vice versa. 
   %  
   % STMINFO: [T SIMIDX X Y]  
   %   
   % X,Y is position of stimulus or NaN if not given. 
   %
   % STMIDX is 
   %       0 :  Adaptation time
   %       1 :  Begin Cue 
   %       2 :  End Cue 
   %       3 :  Stimulation period
   %       4 :  Gap period
   %       5 :  Testing period
     
     oldstmIdx = self.stmIdx;
     trainingTime = 2*self.signalTime + self.stmTime;
     roundTime = trainingTime + self.gapTime;
            
     x = nan(length(tracks),1);
     y = nan(length(tracks),1);

     % set the stmidx
     if t < self.adaptationTime
       self.stmIdx = ID_ADAPTATION;
     else
       % within round
       tt = t - self.adaptationTime;
       ttmod = mod(tt,roundTime);
                
       if ttmod<self.signalTime
         self.stmIdx = ID_BEGINCUE;    
         self.testingif = rand<self.testingProb;
       elseif  ttmod<self.signalTime + self.stmTime
         if ~self.testingif
           self.stmIdx = ID_STIMULUS;
           
         else
           self.stmIdx = ID_TEST;
         end
       elseif ttmod<2*self.signalTime + self.stmTime
         self.stmIdx = ID_ENDCUE;
       else
         self.stmIdx = ID_GAP;
       end
     end
            
     % plot the background only if stmidx changed       
     if oldstmIdx~=self.stmIdx
       verbose('Switch Stimulus Plane %d -> %d\r',oldstmIdx,self.stmIdx);
       switch self.stmIdx
         case ID_BEGINCUE
           self.plotVPlane(0,self.colBeginCue,self.colBeginCue);
         case ID_ENDCUE
           self.plotVPlane(0,self.colEndCue,self.colEndCue);
         otherwise
           self.plotVPlane(0,self.colBackgound,self.colBackgound);
       end
     end
          
     if self.stmIdx==ID_STIMULUS
       [x,y] = self.plotStimulus(tracks,framesize,t);
     end       
     if self.stmIdx==ID_STIMULUS || self.stmIdx==ID_TEST
       % markers and end
       plotVLine(self,self.sideMarkerPos,self.sideMarkerWidth,self.colSideMarkerLeft);
       plotVLine(self,1-self.sideMarkerPos,self.sideMarkerWidth,self.colSideMarkerRight);
     end
       
     % flip each step !
     timestamp = self.flip();
            
     % save stm info
     for i = 1:length(tracks)
       tracks(i).stmInfo = [t, timestamp, self.stmIdx, x(i), y(i)];
     end
   end
   
        
   function bool = isFinished(self,t)
     bool = t>self.adaptationTime+(self.stmTime+2*self.signalTime+self.gapTime)*self.nRound;
   end
   
  end
  
end



    
    