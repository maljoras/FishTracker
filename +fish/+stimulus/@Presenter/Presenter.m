classdef Presenter < handle;
  
  
  
  properties 
    screen = 1;
    defaultColor = [1,0.2,1];

    tmax = Inf;
    
    xreversed = true; % whether screen versus camera are reversed in x
    yreversed = true; % or y
    screenBoundingBox = [];
  
    muteAllFlipping = false; % CAUTION: Turns off all flipping (for debugging only)
  
    usePredFishId = false;  % wether to use the ID predicted by
                           % DAG. (DEPRECIATED. Results in too many
                           % switches)
    borderWidth = 0;
    colBorder= [1,1,1];
    colBackground= [0,0,0];
  
    progressBar = 1;
    progressBarUpdateInt = 50; % in frames
  end
  
  
    
  properties(SetAccess=private);
  
    windowSize = [];
    ifi = []; % flip interval
    windowRect = [];
    window = []; 
    textureIdx =  [];
  end
  
  properties(SetAccess=protected);
    IDX_FISHID = 3;
    IDX_XY = 1:2;
  end

  properties(SetAccess=private,GetAccess=private);
    topPriorityLevel = [];
    progressBarH  = [];
    istep = 0;
  end
  
  methods 
    
    function self = Presenter(varargin)

      self = self@handle();
      self.setOpts(varargin{:});
      
    end
    
    
    function reset(self);
    % RESET(SELF) called in initTracking. Can be overloaded to reset the
    % state of the stimulus presenter
    
      self.istep = 0;
      if isinf(self.tmax)
        self.progressBar = 0;
      end
      self.initProgressBar();
    end
    
    function closeProgressBar(self)

      if ishandle(self.progressBarH)
        close(self.progressBarH);
      end
      self.progressBarH = [];        
    end
    
    function initProgressBar(self)
      if self.progressBar
        self.closeProgressBar();
        global GLOBAL_CLOSE_REQ;
        GLOBAL_CLOSE_REQ = 0;
        self.progressBarH = waitbar(0,'Stimulus progress.. close to stop stimulation.',...
                                    'Name',sprintf('%s progress..',class(self)));
      
        set(self.progressBarH,'CloseRequestFcn',@closereqfun);
      end
    end
    
    
    function updateProgressBar(self,t,tmax);
      tt = datevec(seconds(tmax-t));
      waitbar(t/tmax,self.progressBarH,...
              sprintf('Estimated time: %1.0fh %1.0fm %1.1fs',tt(4),tt(5),tt(6)));
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
      
      if self.progressBar
        if bool
          self.closeProgressBar();
        else
          global GLOBAL_CLOSE_REQ
          bool =  bool || GLOBAL_CLOSE_REQ;
        end
      end
    
    end
    
    function [xx] =toScreenX(self,normx);
      if self.xreversed
        xx = (1-normx)*self.windowSize(1); 
      else
        xx = (normx)*self.windowSize(1); 
      end
    end

    function [normx] =fromScreenX(self,xx);
      if self.xreversed
        normx = -xx/self.windowSize(1)+1;
      else
        normx = xx/self.windowSize(1); 
      end
    end

    
    function [yy] =toScreenY(self,normy);
      if self.yreversed
        yy = (1-normy)*self.windowSize(2); 
      else
        yy = normy*self.windowSize(2); 
      end
    end
    
    function [normy] =fromScreenY(self,yy);
      if self.yreversed
        normy = -yy/self.windowSize(2)+1;; 
      else
        normy  = yy/self.windowSize(2) ;
      end
    end
    
    function [wyy] =toScreenHeight(self,wy);
      wyy = wy*self.windowSize(2); 
    end
    
    function [wy] =fromScreenHeight(self,wyy);
      wy  = wyy/self.windowSize(2);
    end

    function [wxx] =toScreenWidth(self,wx);
      wxx = wx*self.windowSize(1); 
    end
    
    function [wx] =fromScreenWidth(self,wxx);
      wx = wxx/self.windowSize(1); 
    end

    function wrect = toScreenRect(self,rect);
    % converted matlab type of rectangle in norm coordinates to the
    % expected format of PsychToolbox
      
      if size(rect,2)==1
        rect = rect';
      end
      
      x = self.toScreenX(rect(:,1));
      y = self.toScreenY(rect(:,2));
      wx = self.toScreenWidth(rect(:,3));
      wy = self.toScreenHeight(rect(:,4));
      
      wrect = [x,y,x+wx,y+wy];

      if self.xreversed
        wrect(:,[1,3]) = [x-wx,x];
      end
      if self.yreversed
        wrect(:,[2,4]) = [y-wy,y];
      end
      
    end

    function rect = fromScreenRect(self,wrect);
    % converted  format of PsychToolbox to  norm coordinates 
      
      if size(wrect,2)==1
        wrect = wrect';
      end

      x = wrect(:,1);
      y = wrect(:,2);
      wx = wrect(:,3)-x;
      wy = wrect(:,4)-y;
      
      if self.xreversed
        x = wrect(:,3);
        wx = -wrect(:,1) + x;
      end
      if self.yreversed
        y = wrect(:,4);
        wy = -wrect(:,2) + y;
      end
      x = self.fromScreenX(x);
      y = self.fromScreenY(y);
      wx = self.fromScreenWidth(wx);
      wy = self.fromScreenHeight(wy);
      
      rect = [x,y,wx,wy];
      
    end


    function fishIds = getFishIdsFromTracks(self,tracks)
      if self.usePredFishId
        fishIds = [tracks.predFishId];
      else
        fishIds = [tracks.fishId];
      end
    end
    
    function convertedbbox = toScreenBbox(self,bbox)
    % converts the bounding box from matlab to coordinates of
    % PsychToolbox. CAUTION: for converting from norm coordinates use
    % TOSCREENRECT
    
      sbbox = self.screenBoundingBox(:)'; % bbox of the overall screen

      nbbox = zeros(size(bbox,1),4);
      nbbox(:,1:2) =  bsxfun(@rdivide,bsxfun(@minus,bbox(:,1:2),sbbox(1,1:2)),sbbox(1,3:4));
      nbbox(:,3:4) = bsxfun(@rdivide,bbox(:,3:4),sbbox(1,3:4));
    
      convertedbbox = self.toScreenRect(nbbox);
    end

    function convertedbbox = fromScreenBbox(self,bbox)
    % converts the bounding box from PsychToolbox coordinates to
    % fishtracker coordinates. Reverses TOSCREENBBOX

      sbbox = self.screenBoundingBox(:)'; 
      nbbox = self.fromScreenRect(bbox);

      convertedbbox = zeros(size(nbbox,1),4);
      convertedbbox(:,3:4) = bsxfun(@times,nbbox(:,3:4),sbbox(1,3:4));
      convertedbbox(:,1:2) =  bsxfun(@plus,bsxfun(@times,nbbox(:,1:2),sbbox(1,3:4)),sbbox(1,1:2));

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

