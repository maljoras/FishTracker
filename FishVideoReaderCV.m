classdef FishVideoReaderCV < FishVideoReader;
% implements the reader using the matlab included video reader
  
  
  properties 
    reader = [];
    originalif = false;
    scale = [1,1,1]/3;
    delta = 0;  

  end
  
  methods
    
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
      nextFrame = floor(time*self.frameRate);
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
    

    
    
    function self = FishVideoReaderCV(vid,trange,varargin) 
    % SELF = FISHVIDEOREADERCV(VID,TRANGE,...) 

      if ~hasOpenCV()
        error('Cannot find mexopencv toolbox..');
      end

      if nargin==1
        trange = [];
        varargin = {};
      end

      self@FishVideoReader(vid,varargin{:});
      
      self.reader = FishVideoCapture(vid);
      self.init();
      self.timeRange = trange;

      %self.reset();
      self.verbose();

    end
    
  end
end
