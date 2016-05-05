classdef PresenterOnlineLearningCue < fish.stimulus.Presenter;
  
  properties
    signalTime = 0;  %time of begining and ending signal (in seconds)
    stmTime = 5; %time of stimulus (in seconds)
    gapTime =10; % time of gap after stim
    adaptationTime =1;%600; % time at the beginning (in seconds
    
    
    colBackground = [0,0,0]; % background color (RGB [0,0,0] for black)
    colBeginCue = [1,0,0];% begining signal (RGB[255,0,0]for red)
    colEndCue = [0,1,0]; % ending signal(RGB[0,255,0]for green)

    sideMarkerCol = [1,1,0;1,0,1];
    sideMarkerPos = 0.05;
    sideMarkerWidth = 10; % in pix
    
    midline = 0.5;  % position of the line form 0..1

    testingProb = 0;
    lrSwitchProb =  1; % 0: no change, 1:switching each stmround
    
    % this function is applied to stimulate fishIds either on the left (<0)
    % or on the right (>0) side or not stm (0)
    funLRStm = @(fishIds) -(fishIds>0);
    funCueBlending  = @(t,signalTime) t/signalTime;
    funRoundStmPause  = @(iround) 1; % always stimulus

    % the two background texture
    textureCutoff = [0.01,0.0001]; % right left
    textureOri = [pi,pi/2];
    textureCol = [1,1,1;1,1,1]; % RGB
    
    stmLambda = 0.5; % this eq to dt*PosissonRate  !
    stmBkgType = 'texture'; % one of 'none','plane','texture','border'
    stmBkgBrightness = 0.8; % from 0..1  
    stmSize =100;
    stmCol = [1,1,1]; % stimulus color (RGB [1,1,1] for white)    

    lrif = 1; % start first stim LR (1) or RL(0)
    
  end

  properties(SetAccess=private)
    ID_ADAPTATION = 0;
    ID_BEGINCUE = 1;
    ID_ENDCUE = 2;
    ID_STIMULUS_LR = 3;
    ID_STIMULUS_RL = 3 + 8;
    ID_GAP = 4;
    ID_TEST_LR = 5;
    ID_TEST_RL = 5 + 8;
    ID_PAUSE = 6;
    
    roundTime = 0;
    iround = 0;
  end
  
  properties(SetAccess=private,GetAccess=private)
    testingif = 0;

    stmIdx = -1;
    stmState = [];
    textures = [];
  end
    

  
  methods
    
    function init(self,varargin)
      
      init@fish.stimulus.Presenter(self,varargin{:});
      
      self.initStmBkg();
      updateRoundTime(self);
    end

    
    function initStmBkg(self)
    % inits all possible backgrounds
      
      %  texture background
      l =max(length(self.textureOri),length(self.textureCutoff));
      assert(l>1); % at least 2 textures
      for i = 1:l
        ori = self.textureOri(min(i,end));
        cut = self.textureCutoff(min(i,end));
        col = shiftdim(self.textureCol(min(i,end),:),-1);
        assert(numel(col)==3);
        
        texture = fish.helper.makeTexture('f2ori',max(self.windowSize),'ori',ori,...
                              'cutoff',cut,'contrast',1);
        texture = texture(1:self.windowSize(1)/2,1:self.windowSize(2))';
        coltexture = bsxfun(@times,texture,col*self.stmBkgBrightness);
        self.textures(i) = self.initTexture(coltexture);
      end
    end
    
      
    function plotStmBkg(self,type,lrswitch);

      switch lower(type)
        
        case 'none'
          plotVPlane(self,self.midline,self.colBackground,self.colBackground);
        case 'texture'
          
          textures = self.textures;
          if lrswitch
            textures = textures([2,1]);
          end
          self.drawTexture(textures(1),[self.midline,0,1-self.midline,1]);
          self.drawTexture(textures(2),[0,0,self.midline,1]);
          
        case {'border','plane'}
          col = {self.sideMarkerCol(1,:)*self.stmBkgBrightness,...
                 self.sideMarkerCol(2,:)*self.stmBkgBrightness};
          if lrswitch
            col = col([2,1]);
          end
          if  strcmp(lower(type),'border')
            plotVLine(self,self.sideMarkerPos,self.sideMarkerWidth,col{1});
            plotVLine(self,1-self.sideMarkerPos,self.sideMarkerWidth,col{2});
          else
            plotVPlane(self,self.midline,col{1},col{2});
          end
          
          
        otherwise;
          error('do not know background stimulus type');
      end
    end
    
    
    
    function plotCue(self,id,tsignal)
      bri  = self.funCueBlending(tsignal,self.signalTime);
      if id==self.ID_BEGINCUE 
        self.plotVPlane(0,self.colBeginCue*bri,self.colBeginCue*bri);
      else
        self.plotVPlane(0,self.colEndCue*bri,self.colEndCue*bri);
      end
    end
    
    function plotVarDot(self,x,y,inSize,inColor)
    % plots a dot in normalized coordinates.  
      
      xx = self.toScreenX(x);
      yy = self.toScreenY(y);
      w = self.toScreenWidth(self.midline);

      for i = 1:length(xx)
        s2 = inSize(min(i,end))/2;
        s2 = min(abs(xx(i)-w),s2);
        rect = [xx(i)-s2,yy(i)-s2,xx(i)+s2,yy(i)+s2];
        Screen('FillOval', self.window, inColor(min(i,end),:)*255, rect);
      end
    end

   function [x,y] = plotStimulus(self,x,y,t,fishIds,lrswitch)

     if length(x) ~= length(self.stmState)
       % should only happen once in the beginning. (nfish is constant)
       self.stmState = ones(size(x));
     end

     % border
     %b = self.stmSize/self.windowSize(1)/2;
     %x = max(x,b);
     %x = min(x,1-b);

     xleft = x<self.midline; % definition of left
     
     %% fetch stimulate whether left or right stm
     if lrswitch
       stmleft = self.funLRStm(fishIds(:))>0;
       stmright = self.funLRStm(fishIds(:))<0;
     else
       stmleft = self.funLRStm(fishIds(:))<0;
       stmright = self.funLRStm(fishIds(:))>0;
     end

     
     %% change stimulus state with poisson characeteristics
     if self.stmLambda>0
       msk = self.stmLambda>rand(size(self.stmState));
       self.stmState(msk) = ~self.stmState(msk);
     end

     left = xleft & stmleft & self.stmState(fishIds);
     self.plotVarDot(x(left),y(left),self.stmSize,self.stmCol);

     right = ~xleft & stmright & self.stmState(fishIds);
     self.plotVarDot(x(right),y(right),self.stmSize,self.stmCol);

     nostmmsk = ~(right | left);
     x(nostmmsk) = NaN;
     y(nostmmsk) = NaN;
   end
   
   function updateRoundTime(self);
     trainingTime = 2*self.signalTime + self.stmTime;
     self.roundTime = trainingTime + self.gapTime;
   end
   
   function stmInfo = stepStimulus(self,x,y,t,fishIds)
   % stmInfo = stepStimulus(self,x,y,t,fishIds) this function will be
   % called from fish.stimulus.Presenter/step once for each frame. It
   % plots an stimulus below the fish if fishID odd and on left side
   % and vice versa.
   %  
   % STMINFO: [STMIDX X Y T ]  
   %   
   % X,Y is position of stimulus or NaN if not given. 
   %
   % STMIDX is 
   %       0 :  Adaptation time
   %       1 :  Begin Cue 
   %       2 :  End Cue 
   %       3 :  Stimulation period (+ 8 LR switch)
   %       4 :  Gap period
   %       5 :  Testing period (+ 8 LR switch)
   %       6 :  Pause (if funRoundStmPause(iround)==0)  
            
     %% set the stmidx
     oldstmIdx = self.stmIdx;
     lastRound = self.iround;
     
     if t < self.adaptationTime
       self.stmIdx = self.ID_ADAPTATION;
     else

       % within round
       tt = t - self.adaptationTime;
       ttmod = mod(tt,self.roundTime);
       
       self.iround = floor(tt/self.roundTime)+1;
       
       if self.funRoundStmPause(self.iround)
       
         if ttmod<self.signalTime
           self.stmIdx = self.ID_BEGINCUE;    
         elseif  ttmod<self.signalTime + self.stmTime
           % stimulus time
           if self.iround~=lastRound
             % new stim round
             self.testingif = rand<self.testingProb;
             if rand<self.lrSwitchProb
               self.lrif = ~self.lrif;
             end
           end
           
           if self.lrif
             if ~self.testingif
               self.stmIdx = self.ID_STIMULUS_LR;
             else
               self.stmIdx = self.ID_TEST_LR;
             end
           else
             if ~self.testingif
               self.stmIdx = self.ID_STIMULUS_RL;
             else
               self.stmIdx = self.ID_TEST_RL;
             end
           end
       
         elseif ttmod<2*self.signalTime + self.stmTime
           self.stmIdx = self.ID_ENDCUE;
         else
           self.stmIdx = self.ID_GAP;
         end
       else
         self.stmIdx = self.ID_PAUSE;    
       end
       
     
     end
       
     %% plotting
     switch self.stmIdx
       case {self.ID_BEGINCUE,self.ID_ENDCUE}
         self.plotCue(self.stmIdx,min(abs(ttmod-self.stmTime-self.signalTime),ttmod));
       case {self.ID_GAP,self.ID_ADAPTATION,self.ID_PAUSE}
         self.plotVPlane(0.5,self.colBackground,self.colBackground);
       otherwise
         % during stimulus
         lrswitch = self.stmIdx==self.ID_TEST_RL || self.stmIdx==self.ID_STIMULUS_RL;
         self.plotStmBkg(self.stmBkgType,lrswitch);

         if self.stmIdx==self.ID_STIMULUS_LR || self.stmIdx==self.ID_STIMULUS_RL
           [x,y] = self.plotStimulus(x,y,t,fishIds,lrswitch); 
         end
                   
     end

     % save stm info
     stmInfo = zeros(length(fishIds),4);
     stmInfo(:,1) = self.stmIdx;
     stmInfo(:,2:3) = [x,y];
     stmInfo(:,4) = t;
   end
   
        
   function bool = isFinished(self,t)
     bool = t>self.adaptationTime+(self.stmTime+2*self.signalTime+self.gapTime)*self.nRound;
   end
   
  end
  
end



    
    