classdef Presenter < handle;
  
  
  
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
    windowRect = [];
    window = []; 
    textureIdx =  [];

  end

  properties(SetAccess=private,GetAccess=private);
    topPriorityLevel = [];
  end
  
  methods 
    
    function self = Presenter(varargin)

      self = self@handle();
      self.setOpts(varargin{:});
      
    end

    
    function setOpts(self,varargin)

      if length(varargin)==1 && isstruct(varargin{1});
        for f = fieldnames(varargin{1})'
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
            error('expect valid property name fish.stimulus.Presenter');
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
  
    
    function bool = isFinished(self,t);
      bool = t>self.tmax;
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

