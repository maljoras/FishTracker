classdef FishVideoReader < handle;
  
  
  
  properties 
    frameFormat = 'RGBU';
    grayFormat = 'GRAY';
    timeRange  = [];
    originalif = false;

  end
  

  properties(Access=protected);
    reader = [];
    scale = [1,1,1]/3;
    delta = 0;  
  end
  
    
  properties(SetAccess=private);

    currentFrame = 0;
    currentTime = 0;
    frame = [];
    oframe = [];
    videoFile = '';
    frameRate = [];
    duration = 0;
    nFrames = 0;
    frameSize = [];

    tframe = [];
  end
  
  
  methods (Access=protected)
    % needs to be overloaded

    function dur = a_getDuration(self);
      dur = self.a_getNFrames()/self.a_getFrameRate();
    end
    
    function nframes = a_getNFrames(self);
      nframes = get(self.reader,'FrameCount');
      if nframes==0
        % strange. cannot find how many frames. Assume Inf
        nframes = Inf;
      end
      
    end
    
    function frameRate = a_getFrameRate(self);
      frameRate = get(self.reader,'FPS');

      if isnan(frameRate) % strasnge things happen
        f = VideoReader(self.videoFile);
        frameRate = f.FrameRate;
        set(self.reader,'FPS',frameRate);
        clear('f')
      end
      
    end
    
    function frameSize= a_getFrameSize(self);
      frameSize = [get(self.reader,'FrameHeight'),get(self.reader,'FrameWidth')];
    end
    
    function [frame oframe] = a_readScaledSFrame(self,scale,delta);      
      if self.originalif
        [frame oframe]= self.reader.readScaledS(scale,delta);
      else
        [frame]= self.reader.readScaledS(scale,delta);
        oframe = [];
      end
    end
    
    function [frame oframe] = a_readScaledUFrame(self,scale,delta);      
      if self.originalif
        [frame oframe] = self.reader.readScaledU(scale,delta);
      else
        frame = self.reader.readScaledU(scale,delta);
        oframe = [];
      end
      
    end

    function [frame oframe] = a_readUFrame(self);      
      if self.originalif
        [frame oframe] = self.reader.read();
      else
        frame = self.reader.read();
        oframe = [];
      end
      
    end
    
    function [frame oframe] = a_readGraySFrame(self);      
      if self.originalif
        [frame oframe] = self.reader.readGraySingle();
      else
        frame = self.reader.readGraySingle();
        oframe = [];
      end
    end
    
    function [frame oframe] = a_readGrayUFrame(self);      
      if self.originalif
        [frame oframe] = self.reader.readGray();
      else
        frame = self.reader.readGray();
        oframe = [];
      end
    end
    
    function [frame oframe] = a_readInvertedGraySFrame(self);      
      if self.originalif
        [frame oframe] = self.reader.readInvertedGraySingle();
      else
        frame = self.reader.readInvertedGraySingle();
        oframe = [];
      end
    end
    
    function [frame oframe] = a_readInvertedGrayUFrame(self);      
      if self.originalif
        [frame oframe] = self.reader.readInvertedGray();
      else
        frame  = self.reader.readInvertedGray();
        oframe = [];
      end
    end
    
    function [frame oframe] = a_readSFrame(self);      
      if self.originalif
        [frame oframe]= self.reader.readSingle();
      else
        frame = self.reader.readSingle();
        oframe = [];
      end
    end
    
    function bool = a_hasFrame(self);
      nextFrame = get(self.reader,'PosFrames');
      bool = nextFrame<self.nFrames; % starts from 0!
    end
    
    
    function a_delete(self);
      self.reader.delete();
      self.reader = [];
    end
    
    
    function a_setCurrentTime(self,time);
    % time is given in seconds
      nextFrame = max(floor(time*self.frameRate),0);
      pos = get(self.reader,'PosFrames');
      if nextFrame~=pos
        set(self.reader,'PosFrames',nextFrame);
      end
    end

    function playMsk(self);
      vp = vision.VideoPlayer('Name',self.videoFile);
      cont = self.hasFrame();
      fgbg = cv.BackgroundSubtractorMOG2()
      %fgbg.Dist2Threshold = 500;
      fgbg.DetectShadows = 0;
      %fgbg.ShadowValue = 0.5;
      while cont
        self.readFrameFormat('GRAYU');
        msk = fgbg.apply(self.frame);
        step(vp,msk);
        cont = self.hasFrame() && isOpen(vp);
      end
      release(vp);
      self.reset();
      
    end
    
    
    function a_startReader(self)
      if ~fish.helper.hasOpenCV()
        error('Cannot find mexopencv toolbox..');
      end
      self.reader = fish.core.FishVideoCapture(self.videoFile);
    end
  
    function loadTFile(self,tfile);
    % reads the txt file with the time information from SaveVideo
      tmp = dlmread(tfile);
      self.tframe = tmp(:,3);%-tmp(1,3);
    end
  
    function startReader(self)
      self.a_startReader();
    end


  
  end
  
  
  
    
  methods
  
    
    function self = FishVideoReader(vid,trange,varargin) 
    % constructor
      self = self@handle();

      if ischar(vid) && ~exist(vid,'file')
        error('Video file not found');
      end
      self.videoFile = vid;
      
      if ischar(vid) && exist([vid '.txt'],'file')
        self.loadTFile([vid '.txt']);
      end
      
      
      if nargin>1
        self.timeRange = trange;
      end
      if nargin>2
        nargs = length(varargin);
        if nargs>0 && mod(nargs,2)
          error('expected arguments of the type ("PropName",pvalue)');
        end
        for i = 1:2:nargs
          if ~ischar(varargin{i}) 
            error('expected arguments of the type ("PropName",pvalue)');
          else
            self.(varargin{i}) = varargin{i+1};
          end
        end
      end
      
      self.startReader();
      self.init();

      if nargin==1
        trange = [];
        varargin = {};
      end

      self.timeRange = trange;

      %self.reset();
      self.verbose();

    end
    
      
    function increaseCounters(self,timeStamp);
      self.currentFrame = self.currentFrame + 1;
      if length(self.tframe) < self.currentFrame
        if nargin>1
          self.currentTime = timeStamp;
        else
          self.currentTime = self.currentTime + 1/self.frameRate;
        end
      else
        self.currentTime = self.tframe(self.currentFrame);
      end
    end

    
    function frame = getCurrentFrame(self);
      frame = self.frame;
    end

    function delete(self)
      a_delete(self);
    end
    
    function set.timeRange(self,trange)
      if isempty(trange) || diff(trange)<=0
        self.timeRange= [0,self.duration];
      else
        assert(length(trange)==2);
        if self.duration
          trange(2) = min(trange(2),self.duration);
        end
        trange(1) = max(trange(1),0);
        self.timeRange = trange;
      end
    end
    
    
    function self = init(self)

      self.frameRate = self.a_getFrameRate();
      self.duration = self.a_getDuration();
      self.nFrames = self.a_getNFrames();
      self.frameSize = self.a_getFrameSize();
    end

    
    function oneFrameBack(self)
      if self.timeRange(1)<self.currentTime
        self.setCurrentTime(self.currentTime);
        self.currentFrame = self.currentFrame -1;
      end
    end
    
    function verbose(self)
      if ~iscell(self.videoFile)
        [a,b,c] = fileparts(self.videoFile);
        fish.helper.verbose('%s using "%s" ',upper(class(self)),[b,c]);
        fish.helper.verbose('FPS: %gHz, NFrames: %d, selected range:  %1.1fs-%1.1fs',self.frameRate,self.nFrames,self.timeRange);
      else
        [a,b,c] = fileparts(self.videoFile{2});
        fish.helper.verbose('%s grabbing from camera %d (%gHz) and writing to "%s" ',upper(class(self)),self.videoFile{1},self.frameRate,[b,c]);
      end
      
    end
    

    function setCurrentTime(self,time)
      assert(time<self.timeRange(2) && time>=self.timeRange(1))
      time = max(time,0);
      
      if isempty(self.tframe)
        self.a_setCurrentTime(time);
        self.currentTime = time-1/self.frameRate;
        self.currentFrame = 0;
      else
        self.currentFrame = find(self.tframe<=time,1,'last')-1;
        self.a_setCurrentTime(self.currentFrame/self.frameRate);
        self.currentTime = self.tframe(self.currentFrame+1);
      end

      self.frame = [];
      self.oframe = [];
    end
    
    function time = getCurrentTime(self)
      time = self.currentTime;
    end
    
    function bool = hasFrame(self);
      if ~isempty(self.timeRange)
        bool = self.currentTime<=self.timeRange(2);
        bool = bool & self.a_hasFrame();
      else
        bool = false;
      end
    end
    
    function reset(self)
    % resets the counters and deleted the current frame
      self.setCurrentTime(self.timeRange(1));
      
      self.frame = [];
      self.oframe = [];
    end

    
    function varargout = readFrameFormat(self,format);
    % read increases frame counter. if ORGINALIF also returns the
    % original frame
      switch format
        case 'RGBS'
          [self.frame self.oframe]= a_readSFrame(self);
        case 'RGBU'
          [self.frame self.oframe] = a_readUFrame(self);
        case 'GRAYU'
          [self.frame self.oframe] = a_readGrayUFrame(self);
        case 'GRAYS'
          [self.frame self.oframe] = a_readGraySFrame(self);
        case 'IGRAYU'
          [self.frame self.oframe] = a_readInvertedGrayUFrame(self);
        case 'IGRAYS'
          [self.frame self.oframe] = a_readInvertedGraySFrame(self);
        case 'SCALEDU'
          [self.frame self.oframe] = a_readScaledUFrame(self,self.scale,self.delta);
        case 'SCALEDS'
          [self.frame self.oframe] = a_readScaledSFrame(self,self.scale,self.delta);
        otherwise
          error('Format not known');
      end
        
      self.increaseCounters();

      if nargout
        varargout{1} = self.frame;
        if self.originalif
          varargout{2} = self.oframe;
        else
          varargout{2} = [];
        end
      end

    end
    
    function setToScaledFormat(self,scale,delta)
      self.scale = scale;
      self.delta = delta;
      self.frameFormat = ['SCALED',self.frameFormat(end)];
    end
    
    function setToRGBFormat(self)
      self.frameFormat = ['RGB',self.frameFormat(end)];
    end
    
    
    function varargout = readFrame(self);
    % reads a frame format 
      [self.frame self.oframe] = self.readFrameFormat(self.frameFormat);
      if nargout
        varargout{1} = self.frame;
        if self.originalif
          varargout{2} = self.oframe;
        else
          varargout{2} = [];
        end
      end
    end
    
    
    
    function play(self)
      vp = vision.VideoPlayer('Name',self.videoFile);
      cont = self.hasFrame();
      while cont
        self.readFrame();
        step(vp,self.frame);
        cont = self.hasFrame() && isOpen(vp);
      end
      release(vp);
      self.reset()
    end
    
  end
  
    
  
end
