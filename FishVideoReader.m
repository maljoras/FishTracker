classdef FishVideoReader < handle;
  
  
  
  properties 
    frameFormat = 'RGBU';
    grayFormat = 'GRAY';
    timeRange  = [];
  end
  

  properties(Abstract);
    reader 
    originalif 
    scale ;
    delta ;  
    
  end
  
    
  properties(SetAccess=private);

    currentFrame = 1;
    currentTime = 0;
    frame = [];
    oframe = [];
    videoFile = '';
    frameRate = [];
    duration = 0;
    nFrames = 0;
    frameSize = [];

  end
  
  
  methods (Abstract)
    % needs to be overloaded
    
   dur       = a_getDuration(self)
   nFrames   = a_getNFrames(self);
   frameRate = a_getFrameRate(self);
   frameSize = a_getFrameSize(self);
   
   [frame,oframe] = a_readUFrame(self);      
   [frame,oframe] = a_readSFrame(self);      

   [frame,oframe] = a_readGrayUFrame(self);      
   [frame,oframe] = a_readGraySFrame(self);      

   [frame,oframe] = a_readInvertedGrayUFrame(self);      
   [frame,oframe] = a_readInvertedGraySFrame(self);      

   [frame,oframe] = a_readScaledUFrame(self,scale,delta);      
   [frame,oframe] = a_readScaledSFrame(self,scale,delta);      


   a_setCurrentTime(self,time);
   bool = a_hasFrame(self);

   a_delete(self);
  end
  
  
    
  methods
  
    
    function self = FishVideoReader(vid,trange,varargin) 
    % constructor
      self = self@handle();

      self.videoFile = vid;
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
    end
    
      
    function increaseCounters(self);
      self.currentTime = self.currentTime + 1/self.frameRate;
      self.currentFrame = self.currentFrame + 1;
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
      [a,b,c] = fileparts(self.videoFile);
      verbose('%s using "%s" ',upper(class(self)),[b,c]);
      verbose('FPS: %gHz, NFrames: %d, selected range:  %1.1fs-%1.1fs',self.frameRate,self.nFrames,self.timeRange);
    end
      
    

    function setCurrentTime(self,time)
      assert(time<self.timeRange(2) && time>=self.timeRange(1))
      time = max(time,0);
      self.a_setCurrentTime(time);
      self.currentTime = time-1/self.frameRate;
      self.frame = [];
      self.oframe = [];
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
      self.currentFrame = 0;
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
