classdef PresenterTrackTextures < xy.stimulus.Presenter;
  
  properties

    adaptationTime = 0;

    
    stmOnInt = 0.5; % this means average  switch interval in sec;
    stmOnIntCV = 0.1; % needs to be >0. near 0 means almost not randomness
                      % and regular intervals.
    stmOffInt = 2; % this means average  switch interval in sec;
    stmOffIntCV = 0.1; % needs to be >0. near 0 means almost not randomness
    
    stmSize =100;
    stmCol = [1,1,1]; % stimulus color (RGB [1,1,1] for white)    

    stmBorderThres = 0; % minimal distance to border for stimulation
    stmVelThres = 0;
    
    
    %% mean parmeters. Can be vector with elements for each identityId
    stmSizeFactor = 1; 
    stmShift = 50; % in px of xy.tracker frame
    stmShiftOri = 0; % in RAD=0 means in front

    %% randomness. 
    % coefficient of varaiton of gamma-dist (SCALAR)
    % 0 turns off randomness. 1 means var equals mean. 

    stmSizeFactorCV = 1; % >1 higher variance. (var = mu*scale)
    stmShiftCV = 1;  
    stmShiftOriSTD = -1; % negative for uniform otherwise gaussian

  end

  properties(SetAccess=protected)
    ID_ADAPTATION = 0;
    ID_STIMULUS = 1;

    IDX_STMIDX = 1;
    IDX_T = 5;
    IDX_STATE = 6;
    IDX_MSK = 7;
    IDX_SHIXYT =8;
    IDX_SHIXYTORI = 9;
    IDX_SIZEFACTOR = 10;
    IDX_COL = 11:13;
    IDX_BBOX = 14:17;
    IDX_VELXY = 18:19;
    
    NIDX = 19;
  end
  
  
  properties(SetAccess=private,GetAccess=private)
    stmState = [];
    nextStateChange = [];
    
    stmIdx = -1;
    trackTextures = [];
    trackBbox = [];
    trackVelocity = [];
    trackOrientation = [];
  
    curStmSizeFactor = [];
    curStmCol = [];
    curStmShift = [];
    curStmShiftOri = [];
    
    lastt = 0;
  end
    
  methods(Access=private);
    function bboxes = transformBbox(self,bboxes,velocity,stmSizeFactor,stmShift,stmShiftOri);
      
      % enlarging
      bboxes(:,[1,2]) = bboxes(:,[1,2])+bsxfun(@times,bboxes(:,[3,4]),(1-stmSizeFactor(:))/2);
      bboxes(:,[3,4]) = bsxfun(@times,bboxes(:,[3,4]),stmSizeFactor(:));
      
      % get shift orientation relative to  velocity
      vel = velocity;

      ori = atan2(vel(:,1),vel(:,2));      

      ori =  ori  + stmShiftOri(:);
      
      % shift box
      bboxes(:,1) = bboxes(:,1) + sin(ori).*stmShift(:);
      bboxes(:,2) = bboxes(:,2) + cos(ori).*stmShift(:);
      
    end
  end
  
  
  methods
    
    function self = PresenterTrackTextures(varargin)

      self = self@xy.stimulus.Presenter();
      self.IDX_XY = 2:3;
      self.IDX_IDENTITYID = 4;
    end

    
    function init(self,varargin)
      init@xy.stimulus.Presenter(self,varargin{:});
    end

    function reset(self)

      reset@xy.stimulus.Presenter(self);
      
      self.stmIdx = -1;
      self.stmState = [];
      self.lastt = 0;
    end
    
    function  stmmsk = assignStm(self,x,y,t,identityIds)

      [~,ridx] = sort(identityIds,'ascend');
      sx  = x(ridx);
      sy = y(ridx);
      
      beginning = 0; 
      dt = t - self.lastt;
      if length(sx) ~= length(self.stmState)
        % should only happen once in the beginning. (nbody is constant)
        self.stmState = ones(size(x));
        self.nextStateChange = zeros(size(x));
        beginning = 1;
      end
      oldState = self.stmState;
      

      if self.stmOnInt>0
        stateChanges = find(self.nextStateChange<=t);

        if ~isempty(stateChanges)
          self.stmState(stateChanges) = ~self.stmState(stateChanges);

          offChange = ~self.stmState(stateChanges);
          onChange = ~offChange;

          % currently on 
          idx = stateChanges(onChange);
          self.nextStateChange(idx) = self.nextStateChange(idx) + ...
              gamrnd(self.stmOnInt*ones(size(idx))/self.stmOnIntCV,self.stmOnIntCV);
          idx = stateChanges(offChange);
          self.nextStateChange(idx) = self.nextStateChange(idx) + ...
              gamrnd(self.stmOffInt*ones(size(idx))/self.stmOffIntCV,self.stmOffIntCV);
        end
      end

      newOn = (oldState~=self.stmState) & (self.stmState==1);
      if beginning
        self.assignNewStmValues(sx,sy,t,1:length(sx));
      elseif any(newOn) 
        % newly on. Draw the parameters for nex stimX
        self.assignNewStmValues(sx(newOn),sy(newOn),t,find(newOn));
      end
      self.lastt = t;

      stmmsk = self.stmState;

      if self.stmBorderThres
        stmmsk = stmmsk & (sx>self.stmBorderThres & sx<1-self.stmBorderThres);
        stmmsk = stmmsk & (sy>self.stmBorderThres & sy<1-self.stmBorderThres);
      end

      if any(self.stmVelThres)
        v = sqrt(sum(self.trackVelocity.^2,2));
        stmmsk = stmmsk &  v>self.stmVelThres(1);
        if length(self.stmVelThres)>1
          stmmsk = stmmsk &  v<self.stmVelThres(2);
        end
      end
      
    end
    
    function assignNewStmValues(self,x,y,t,identityIds)

      self.curStmSizeFactor(identityIds) = self.assignStmSizeFactor(x,y,t,identityIds);
      self.curStmCol(identityIds,:) = self.assignStmCol(x,y,t,identityIds);
      self.curStmShift(identityIds) = self.assignStmShift(x,y,t,identityIds);
      self.curStmShiftOri(identityIds) = self.assignStmShiftOri(x,y,t,identityIds);

    end
    
    function [rnd] = grand(self,mus,cv)
    % [RND] = GRAND(SELF,MUS,CV) gamma-rand numbers (same size as MUS)
    % with mean MUS and coefficient of variation CV

      if cv
        rnd = gamrnd(mus/cv, cv);
      else
        rnd = mus;
      end
    end

    
    function [rnd] = orand(self,mus,sig)
    % [RND] = ORAND(SELF,MUS,SIG) normal-rand numbers (same size as MUS)
    % with mean MUS and stddev SIG. if SIG<0 uniform numbers in
    % 0,2pi are produced

      if sig>0
        rnd = randn(size(mus)).*sig+mus;
      elseif sig<0
        rnd = rand(size(mus))*(2*pi);
      else
        rnd = mus;
      end
    end

    
    function [stmSizeFactor] = assignStmSizeFactor(self,x,y,t,identityIds)
      stmSizeFactor = self.grand(self.stmSizeFactor(min(identityIds,end)),self.stmSizeFactorCV);
    end
    
    function [stmCol] = assignStmCol(self,x,y,t,identityIds)
      stmCol = self.stmCol(min(identityIds,end),:);     
    end
    
    function [stmShiftOri] = assignStmShiftOri(self,x,y,t,identityIds)
      stmShiftOri = self.orand(self.stmShiftOri(min(identityIds,end)),self.stmShiftOriSTD); 
    end
    
    function [stmShift] = assignStmShift(self,x,y,t,identityIds)
      stmShift = self.grand(self.stmShift(min(identityIds,end)),self.stmShiftCV); 
    end
    
    function [stmbboxes,stmmsk] = plotStimulus(self,x,y,t,identityIds)

      stmmsk = self.assignStm(x,y,t,identityIds);
      stmmsk = stmmsk & self.curStmSizeFactor(:)>0;

      stmbboxes = nan(length(identityIds),4);
      idx = find(stmmsk(identityIds));

      if ~isempty(idx)
        fids = identityIds(idx);
        bboxes = self.transformBbox(self.trackBbox(idx,:), ...
                                    self.trackVelocity(idx,:), ...
                                    self.curStmSizeFactor(fids), ...
                                    self.curStmShift(fids), ...
                                    self.curStmShiftOri(fids));
        
        stmbboxes(idx,:) = self.toScreenBbox(bboxes);

        for i = idx(:)'
          Screen('drawTexture',self.window,self.trackTextures(i),[],...
                 stmbboxes(i,:),[],[],[],self.curStmCol(identityIds(i),:)*255);
        end
      end      
    end

    
    
    function tracks = step(self,tracks,framesize,t)
    % overloads step to get the textures

      self.saveTrackTextures(tracks);

      % let  parent do the work
      tracks = step@xy.stimulus.Presenter(self,tracks,framesize,t);
    end
    
    
    function stmInfo = stepStimulus(self,x,y,t,identityIds)
    % stmInfo = stepStimulus(self,x,y,t,identityIds) this function will be
    % called from xy.stimulus.Presenter/step once for each frame.
    %  
      
    %% set the stmidx
      oldstmIdx = self.stmIdx;
      stmbbox = nan(length(x),4);
      stmmsk = zeros(length(x),1);

      if t < self.adaptationTime
        self.stmIdx = self.ID_ADAPTATION;
      else
        self.stmIdx = self.ID_STIMULUS;
      end

      if  self.stmIdx == self.ID_STIMULUS;
        % stmbackground
        self.plotBorder(self.borderWidth,self.colBorder,self.colBackground);
        
        %% plotting
        % plot stimulus
        [stmbbox,stmmsk] = self.plotStimulus(x,y,t,identityIds); 
      end
      
      % save stm info
      stmInfo = NaN(length(identityIds),self.NIDX);
      stmInfo(:,self.IDX_STMIDX) = self.stmIdx;
      stmInfo(:,self.IDX_XY) = [x,y];
      stmInfo(:,self.IDX_T) = t;
      stmInfo(:,self.IDX_BBOX) = stmbbox;
      stmInfo(:,self.IDX_IDENTITYID) = identityIds;
      stmInfo(:,self.IDX_VELXY) = self.trackVelocity;

      if self.stmIdx == self.ID_STIMULUS
        stmInfo(:,self.IDX_MSK) = stmmsk(identityIds);
        stmInfo(:,self.IDX_SHIXYT) = self.curStmShift(identityIds);
        stmInfo(:,self.IDX_SHIXYTORI) = self.curStmShiftOri(identityIds);
        stmInfo(:,self.IDX_SIZEFACTOR) = self.curStmSizeFactor(identityIds);
        stmInfo(:,self.IDX_COL) = self.curStmCol(identityIds,:);
        stmInfo(:,self.IDX_STATE) = self.stmState(identityIds);
      end
    end
    

    function saveTrackTextures(self,tracks);

      sbbox = self.screenBoundingBox;   

      for i = 1:length(self.trackTextures)
        Screen('Close', self.trackTextures(i));
      end
      self.trackTextures = zeros(length(tracks),1);
      
      bodyori = zeros(length(tracks),1);
      for i = 1:length(tracks)
        seg = tracks(i).segment;

        if ~isempty(seg)
          img = flipud(fliplr(seg.Image));
          bodyori(i) = seg.Orientation;
        else
          img = 0;
          bodyori(i) = 0;
        end
        self.trackTextures(i) = Screen('makeTexture',self.window,img*255,0,0);
      end

      % transform boxes
      self.trackBbox = cat(1,tracks.bbox);
      self.trackVelocity = cat(1,tracks.velocity);
      self.trackOrientation = -bodyori/180*pi+90;
      
    end

  end
  
end




