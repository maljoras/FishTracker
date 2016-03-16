classdef FishStimulusPresenter < handle;
  
  
  
  properties 
    screen = 1;
    defaultColor = [1,0.2,1];
    tmax = Inf;
    xreversed = true; % whether screen versus camera are reversed in x
    yreversed = true; % or y
    screenBoundingBox = [];
  
    muteAllFlipping = false; % CAUTION: Turns off all flipping (for debugging only)
  end
  
  
    
  properties(SetAccess=private);
  
    windowSize = [];
    ifi = []; % flip interval
  end

  properties(SetAccess=private,GetAccess=private);
    topPriorityLevel = [];
    windowRect = [];
    window = []; 
    textureIdx =  [];
  end
  
  methods 
    
    function self = FishStimulusPresenter(varargin)

      self = self@handle();

      self.setOpts(varargin{:});
      
    end

    
    function setOpts(self,varargin)
    
      if length(varargin)==1 && isstruct(varargin{1});
        for f = fieldnames(varargin{1})
          if isprop(self,f{1})
            self.(f{1}) = varargin{1}.(f{1});
          end
        end
      else
        
        if mod(length(varargin),2)
          error('expect even number of arguments or options structure');
        end
        
        for i = 1:2:length(varargin)
          if ~ischar(varargin{i}) || ~isprop(self,varargin{i});
            error('expect valid property name FishStimulusPresenter');
          end
          self.(varargin{i}) = varargin{i+1};
        end
      end
    
    end     
    
    function setScreenSize(self,frameBoundingBox)
      self.screenBoundingBox = frameBoundingBox;
    end
    
  
    function init(self,wrect)
    % inits the screen and the initial PsychToolBox stuff (needs to be called! )
      if ~exist('Screen') 
        error(['Seems that PsychoTooolbox is  cannot be accessed.']);
      end
      PsychDefaultSetup(2);
      self.release();

      Screen('Preference', 'SkipSyncTests', 1);     
      screens = Screen('Screens');
      if ~ismember(self.screen,screens)
        error(sprintf('Do not find screen %d',self.screen));
      end
      
      % open window
      if ~exist('wrect','var')
        wrect  = [];
      end
      
      black = BlackIndex(self.screen);
      [self.window, self.windowRect] = Screen('OpenWindow', self.screen, black,wrect); 

      self.ifi = Screen('GetFlipInterval', self.window);
      self.topPriorityLevel = MaxPriority(self.window);    
      Screen('Flip', self.window);
      Priority(self.topPriorityLevel);

      Screen('Preference', 'SkipSyncTests', 0);
      self.windowSize = [self.windowRect(3) - self.windowRect(1),...
                        self.windowRect(4) - self.windowRect(2)];

    end

    function release(self)
      Priority(0);       
      Screen('CloseAll');
      sca;
    end
    
    function applyInit(self,stimPresenter)
    % apply the init values from another stim presenter
      self.ifi = stimPresenter.ifi;
      [self.window,self.windowSize,self.windowRect, self.topPriorityLevel] ...
        = stimPresenter.getWindow();
    end
    
    function [window,windowSize,windowRect, topPriorityLevel] = ...
          getWindow(self);
      window = self.window;
      windowSize = self.windowSize;
      windowRect = self.windowRect;
      topPriorityLevel  = self.topPriorityLevel;
    end
  
    
    function timestamp = plotDot(self,x,y,inSize,inColor)
    % plots a dot in normalized coordinates. Width in pixel
      if isempty(x) || isempty(y)
        timestamp = NaN;
        return
      end
      
      if ~exist('inColor','var')
        inColor = self.defaultColor;
      end
      if ~exist('inSize','var')
        inSize = 20;
      end
      assert(length(x)==length(y))

      xx = self.toScreenX(x);
      yy = self.toScreenY(y);

      for i = 1:length(xx)
        s2 = inSize(min(i,end))/2;
        rect = [xx(i)-s2,yy(i)-s2,xx(i)+s2,yy(i)+s2];
        Screen('FillOval', self.window, inColor(min(i,end),:)*255, rect);
      end

      if nargout
        timestamp = self.flip();
      end
      
    end

    function timestamp = plotVLine(self,x,inWidth,inColor)
    % plots a dot in normalized coordinates. Width in pixel
      if ~exist('inColor','var')
        inColor = self.defaultColor;
      end
      if ~exist('inWidth','var')
        inWidth = 10;
      end

      xw = inWidth;
      xw = max(min(xw,10),0.5); % max supported
      
      xx = self.toScreenX(x);

      Screen('DrawLine', self.window, inColor*255,xx,self.toScreenY(0),...
             xx,self.toScreenY(1),inWidth);
      
      if nargout 
        timestamp =self.flip();
      end
      
    end
    
    function timestamp = plotVPlane(self,x,inColor1,inColor2)
    % plots a vertical half plane at x with left size color1 right color2
      if ~exist('inColor1','var')
        inColor1 = self.defaultColor;
      end
      if ~exist('inColor2','var')
        inColor2 = BlackIndex(self.screen);
      end

      self.patch(0,0,x,1,inColor1);
      self.patch(x,0,1,1,inColor2);
      
      if nargout
        timestamp = self.flip();
      end
    end

    function timestamp = plotHPlane(self,y,inColor1,inColor2)
    % plots a horizontal half plane at y with left size color1 right color2
      if ~exist('inColor1','var')
        inColor1 = self.defaultColor;
      end
      if ~exist('inColor2','var')
        inColor2 = BlackIndex(self.screen);
      end

      self.patch(0,0,1,y,inColor1);
      self.patch(0,y,1,1,inColor2);
      
      if nargout
        timestamp = self.flip();
      end
    end
    
    function bool = isFinished(self,t);
      bool = t>self.tmax;
    end
    
         
    
    function stmInfo = stepStimulus(self,x,y,t,fishIds);
    % plots the stimulus (does NOT need a flip). 

      col = zeros(length(fishIds),3);
      for i = 1:length(fishIds)
        col(i,:) = [1-(fishIds(i)==1),(fishIds(i)==2),1 - (fishId(i)-1)/length(fishIds)];
      end
      self.plotDot(x,y,50,col);     
    
      stmInfo = [x,y];
    end
    
    
    function tracks = step(self,tracks,framesize,t)
    % this function will be called from fish.Tracker after each
    % round. Calls the stepStimulus method which should be overloaded!
    
      if isempty(tracks)
        return;
      end
      
      
      if isempty(self.screenBoundingBox)
        sbbox = [1,1,framesize([2,1])-1];
      else
        sbbox = self.screenBoundingBox;
      end

      x = nan(length(tracks),1);
      y = nan(length(tracks),1);
      for i =1:length(tracks)
        x(i) = (tracks(i).centroid(1)-sbbox(1))/sbbox(3);
        y(i) = (tracks(i).centroid(2)-sbbox(2))/sbbox(4);
      end

      
      fishIds = [tracks.fishId];
      stmInfo = self.stepStimulus(x,y,t,fishIds);

      for i =1:length(tracks)
        tracks(i).stmInfo = stmInfo(i,:);
      end
    
      self.flip();
    
    end
      
    function [xx] =toScreenX(self,normx);
      if self.xreversed
        xx = (1-normx)*self.windowSize(1); 
      else
        xx = (normx)*self.windowSize(1); 
      end
    end
    
    function [yy] =toScreenY(self,normy);
      if self.yreversed
        yy = (1-normy)*self.windowSize(2); 
      else
        yy = normy*self.windowSize(2); 
      end
    end
    
    function [wyy] =toScreenHeight(self,wy);
      wyy = wy*self.windowSize(2); 
    end
    function [wxx] =toScreenWidth(self,wx);
      wxx = wx*self.windowSize(1); 
    end

    function wrect = toScreenRect(self,rect);
    % converted matlav type of rectangle in norm coordinates to the
    % expected format of PsychToolbox

      x = self.toScreenX(rect(1));
      y = self.toScreenY(rect(2));
      wx = self.toScreenWidth(rect(3));
      wy = self.toScreenHeight(rect(4));
      
      wrect = [x,y,x+wx,y+wy];

      if self.xreversed
        wrect([1,3]) = [x-wx,x];
      end
      if self.yreversed
        wrect([2,4]) = [y-wy,y];
      end
      
    end
    
    function patch(self,x,y,wx,wy,inColor)
    % plots a patch WITHOUT flipping!

      if ~exist('inColor','var')
        inColor = self.defaultColor;
      end
      xx = self.toScreenX(x);
      yy = self.toScreenY(y);
      wxx = self.toScreenWidth(wx);
      wyy = self.toScreenHeight(wy);

      if self.xreversed
        xx = xx-wxx;
      end
      if self.yreversed
        yy = yy-wyy;
      end

      Screen('FillPoly', self.window, inColor*255,...
             [xx,yy;xx+wxx-1,yy;xx+wxx-1,yy+wyy-1;xx,yy+wyy-1;xx,yy]);

    end

        
    
    function relidx = initTexture(self,mat)
      self.textureIdx(end+1) = Screen('makeTexture',self.window,mat*255);
      relidx = length(self.textureIdx);
    end
    
    function  varargout = drawTexture(self,relidx,rect)
      if relidx>length(self.textureIdx) || relidx<1
        error('Cannot find texture. Forgot to init ?');
      end
      if nargin<3
        rect = [0,0,1,1];
      end

      wrect = self.toScreenRect(rect);

      
      Screen('drawTexture',self.window,self.textureIdx(relidx),[],wrect);
      
      if nargout
        varargout{1} = self.flip();
      end
    end

    
    function timestamp = flip(self,force);
    % Flip to the screen
      if ~self.muteAllFlipping || (nargin>1 && force)
        Screen('DrawingFinished', self.window);
        timestamp = Screen('Flip', self.window);
      end
    end
    
    function timestamp = clear(self);
      timestamp = self.flip();
    end
    
  end
  
end
