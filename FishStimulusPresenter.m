classdef FishStimulusPresenter < handle;
  
  
  
  properties 
    screen = 1;
    defaultColor = [1,0.2,1];
    tmax = Inf;
  end
  
  
    
  properties(SetAccess=private);
  
    windowSize = [];
    ifi = []; % flip interval

  end
  properties(SetAccess=private,GetAccess=private);
    topPriorityLevel = [];
    windowRect = [];
    window = []; 
  end
  
  methods 
    
    function self = FishStimulusPresenter(varargin)

      self = self@handle();
      
      if mod(length(varargin),2)
        error('expect even number of arguments');
      end

      for i = 1:2:length(varargin)
        if ~ischar(varargin{i}) || ~isprop(self,varargin{i});
          error('expect valid property name FishStimulusPresenter');
        end
        self.(varargin{i}) = varargin{i+1};
      end
      
    end
    
    
    function init(self)
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
      black = BlackIndex(self.screen);
      [self.window, self.windowRect] = Screen('OpenWindow', self.screen, black); 

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
    
    
    
    function timestamp = plotDot(self,x,y,inSize,inColor)
    % plots a dot in normalized coordinates. Width in pixel
      if ~exist('inColor','var')
        inColor = self.defaultColor;
      end
      if ~exist('inSize','var')
        inSize = 20;
      end

      xx = x*self.windowSize(1);
      yy = y*self.windowSize(2);

      assert(length(xx)==length(yy))
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

      inWidth = max(min(inWidth,10),0.5); % max supported
      xx = x*self.windowSize(1); % normalize

      Screen('DrawLine', self.window, inColor*255,xx,0,xx,self.windowSize(2),inWidth);
      
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
    
    
    function tracks = step(self,tracks,framesize,t)
    % this function will be called from FishTracker after each round
    
      if isempty(tracks)
        return;
      end
      
      fishId = [tracks.fishId];
      idx = find(fishId==1);
      
      track = tracks(idx);
      sz = framesize;
      for i =1:length(tracks)

        x(i) = tracks(i).centroid(1)/sz(2);
        y(i) = tracks(i).centroid(2)/sz(1);
        col(i,:) = [1-i==1,i==2,1 - (i-1)/length(tracks)];
      end
      % only on one halfplane
      %  idx = x>0.5;
      %x(idx) = [];
      %y(idx) = [];
      %col(idx,:)= [];
      boundary = 0.08;
      x = max(min((x - boundary)/(1-2*boundary),1),0);
      y = max(min((y - boundary)/(1-2*boundary),1),0);
      self.plotDot(1-x,1-y,50,col);
      for i =1:length(tracks)
        tracks(i).stmInfo = [x(i),y(i)];
      end
      self.flip();
    end
      
    function patch(self,x,y,wx,wy,inColor)
    % plots a patch WITHOUT flipping!

      if ~exist('inColor','var')
        inColor = self.defaultColor;
      end

      xx = x*self.windowSize(1); % normalize
      yy = y*self.windowSize(2); % normalize
      wxx = wx*self.windowSize(1);
      wyy = wy*self.windowSize(2);
      Screen('FillPoly', self.window, inColor*255,...
             [xx,yy;xx+wxx-1,yy;xx+wxx-1,yy+wyy-1;xx,yy+wyy-1;xx,yy]);

    end

    function timestamp = flip(self);
    % Flip to the screen
      Screen('DrawingFinished', self.window);
      timestamp = Screen('Flip', self.window);
    end
    
    function timestamp = clear(self);
      timestamp = self.flip();
    end
    
  end
  
end

