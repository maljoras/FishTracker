classdef FishVideoHandlerMex < handle & FishBlobAnalysis & FishVideoReader
%FISHVIDEOHANDER  wrapper class
%
% Class for video reading of the FishTracker
%
%
  
  
  
  properties (SetAccess=private)
    % Object ID
    id;

  end
  properties;
    reader;

  end
  
  
  properties (Dependent)
    plotting
    useScaled
    delta 
    scale
    
    ploting

    minArea 
    maxArea 
    minextent
    maxextent
    minWidth 
    featureheight
    featurewidth
    colorfeature
    originalif  
    
    minContourSize
    
    
% $$$     PosMsec      % Current position of the video file in milliseconds or video capture timestamp.
% $$$     PosFrames    % 0-based index of the frame to be decoded/captured next.
% $$$     AVIRatio     % Relative position of the video file: 0 - start of the film, 1 - end of the film.
% $$$     FrameWidth   % Width of the frames in the video stream.
% $$$     FrameHeight  % Height of the frames in the video stream.
% $$$     FPS          % Frame rate.
% $$$     FourCC       %  4-character code of codec.
% $$$     FrameCount   %  Number of frames in the video file.
% $$$     Format       %  Format of the Mat objects returned by retrieve() .
% $$$     ConvertRGB   %  Boolean flags indicating whether images should be converted to RGB.
% $$$     Mode            %   Backend-specific value indicating the current capture mode.
    
    DetectShadows
    history
% $$$     Dist2Threshold
% $$$     KNNSamples
% $$$     NSamples
% $$$     ShadowThreshold
% $$$     ShadowValue
    
  end

  methods 
    function self = FishVideoHandlerMex(vidname,timerange,opts)
    %VIDEOCAPTURE  Create a new FishVideoHandler object
    %
      if nargin < 1, filename = getVideoFile(); end


      self@FishBlobAnalysis(); 
      self@FishVideoReader(vidname); 
      self.id = FishVideoHandler_(vidname);      
      self.setFishSize([], []);
      self.init();
      self.timeRange = [];
      self.verbose();
      
      if nargin >1
        self.timeRange = timerange;
      end
      
      if nargin >2
        setOpts(self,opts);
      end
      
    end
    
    function setOpts(self,opts)
      
      for f1 = fieldnames(opts)'
        if any(strcmp(f1{1},{'blobs','reader','detector'}))
          for f2 = fieldnames(opts.(f1{1}))'
            if isprop(self,f2{1})
              self.(f2{1}) = opts.(f1{1}).(f2{1});
            else
              warning(sprintf('unknown property: %s',f2{1}));
            end
          end
        end
      end
      
      
    end
    
    
    function delete(self)
    %DELETE  Destructor of VideoCapture object
      FishVideoHandler_(self.id, 'delete');
    end
    
    function [seg,varargout] = step(self)
    % STEP one frame
    %
    % [segments [contours, bwimg, frame]] = vh.step();
    %
    %
      if nargout==1
        [seg] = FishVideoHandler_(self.id, 'step');
      else
        [seg,varargout{1:nargout-1}] = FishVideoHandler_(self.id, 'step');
      end
      self.increaseCounters(); % implicit read frame, so increase
                               % the counters

      % split regions / mser not yet implemented

    end
    
    %%%%%%%%%%%% ABSTRCT CLASSES 
    
    function dur = a_getDuration(self);
      dur = self.a_getNFrames()/self.a_getFrameRate();
    end
    
    function nframes = a_getNFrames(self);
      nframes = FishVideoHandler_(self.id,'capget','FrameCount');
      if nframes==0
        % strange. cannot find how many frames. Assume Inf
        nframes = Inf;
      end
      
    end
    
    function frameRate = a_getFrameRate(self);
      frameRate = FishVideoHandler_(self.id,'capget','FPS');

      if isnan(frameRate) % strasnge things happen
        f = VideoReader(self.videoFile);
        frameRate = f.FrameRate;
        FishVideoHandler_(self.id,'capset','FPS',frameRate);
        clear('f')
      end
      
    end
    
    function frameSize= a_getFrameSize(self);
      
      frameSize = [ FishVideoHandler_(self.id,'capget','FrameHeight'),...
                    FishVideoHandler_(self.id,'capget','FrameWidth')];
    end
    
    function bool = a_hasFrame(self);
      nextFrame = FishVideoHandler_(self.id,'capget','PosFrames');
      bool = nextFrame<self.nFrames; % starts from 0!
    end
    
    
    function a_delete(self);
      FishVideoHandler_(self.id,'delete');
    end
    
    function a_setCurrentTime(self,time);
    % time is given in seconds
      nextFrame = floor(time*self.frameRate);
      pos = FishVideoHandler_(self.id,'capget','PosFrames');
      if nextFrame~=pos
        FishVideoHandler_(self.id,'capset','PosFrames',nextFrame);
      end
    end

    
    function  [frame,oframe] = a_readUFrame(self);      

      bool = self.useScaled;
      self.useScaled = false;

      if self.originalif
        [~,~,frame,oframe] = self.step();
      else
        [~,~,frame] = self.step();
        oframe = [];
      end
      self.useScaled = bool;
      
    end
    

    function [frame,oframe] = a_readSFrame(self);
      [frame,oframe] = a_readUFrame(self);
      frame = single(frame)/255.
      oframe = single(oframe)/255.;
    end

    function  [frame,oframe] = a_readGrayUFrame(self);      
      [frame,oframe] = a_readUFrame(self); 
    end


    function [frame,oframe] = a_readGraySFrame(self);  
      [frame,oframe] = a_readSFrame(self);
    end
    
    function [frame,oframe] = a_readInvertedGrayUFrame(self); 
      [frame,oframe] = a_readGrayUFrame(self); 
      frame = 255-frame;
    end
    
    function [frame,oframe] = a_readInvertedGraySFrame(self);      
      [frame,oframe] = a_readGraySFrame(self); 
      frame = 1-frame;
    end
    
    function  [frame,oframe] = a_readScaledUFrame(self);      
      
      bool = self.useScaled;
      self.useScaled = true;
      if self.originalif
        [~,~,frame,oframe] = self.step();
      else
        [~,~,frame] = self.step();
        oframe = [];
      end
      self.useScaled = bool;
    end

    function  [frame,oframe] = a_readScaledSFrame(self);      
      [frame,oframe] = a_readScaledUFrame(self);
      frame = single(frame)/255.;
    end

    
    %%%%%%
    function value = get.minArea(self)
      value =  FishVideoHandler_(self.id, 'get','minArea');
    end
    function set.minArea(self,value)
      FishVideoHandler_(self.id, 'set','minArea',value);
    end
    function value = get.maxArea(self)
      value =  FishVideoHandler_(self.id, 'get','maxArea');
    end
    function set.maxArea(self,value)
      FishVideoHandler_(self.id, 'set','maxArea',value);
    end

    function value = get.minextent(self)
      value =  FishVideoHandler_(self.id, 'get','minExtent');
    end
    function set.minextent(self,value)
      FishVideoHandler_(self.id, 'set','minExtent',value);
    end
    function value = get.originalif(self)
      value =  FishVideoHandler_(self.id, 'get','colorfeature');
    end
    function set.originalif(self,value)
      FishVideoHandler_(self.id, 'set','colorfeature',value);
    end
    function value = get.colorfeature(self)
      value =  FishVideoHandler_(self.id, 'get','colorfeature');
    end
    function set.colorfeature(self,value)
      FishVideoHandler_(self.id, 'set','colorfeature',value);
    end

    function value = get.maxextent(self)
      value =  FishVideoHandler_(self.id, 'get','maxExtent');
    end
    function set.maxextent(self,value)
      FishVideoHandler_(self.id, 'set','maxExtent',value);
    end
    function value = get.minWidth(self)
      value =  FishVideoHandler_(self.id, 'get','minWidth');
    end
    function set.minWidth(self,value)
      FishVideoHandler_(self.id, 'set','minWidth',value);
    end

    function set.featurewidth(self,value)
      FishVideoHandler_(self.id, 'set','featurewidth',value);
    end
    function set.featureheight(self,value)
      FishVideoHandler_(self.id, 'set','featureheight',value);
    end
    function value = get.featurewidth(self)
      value = FishVideoHandler_(self.id, 'get','featurewidth');
    end
    function value = get.featureheight(self)
      value = FishVideoHandler_(self.id, 'get','featureheight');
    end
    
    function value = get.plotting(self)
      value = FishVideoHandler_(self.id, 'get','plotif');
    end
    
    function set.plotting(self,value)
      FishVideoHandler_(self.id, 'set','plotif',value);
    end

    function setScaledFormat(self,rgbchannel,delta)
      self.useScaled = 1;
      self.delta = delta;
      self.scale = rgbchannel;
    end
      
    function value = get.scale(self)
      value= FishVideoHandler_(self.id, 'getScale');
    end
    function set.scale(self,value)
      FishVideoHandler_(self.id, 'setScale',value(1),value(2),value(3));
    end
    
    function value = get.delta(self)
      value= FishVideoHandler_(self.id, 'get','delta');
    end
    function set.delta(self,value)
      FishVideoHandler_(self.id, 'set','delta',sum(value));
    end
    
    
    
    function set.useScaled(self,value)
      if value
        warning('Scaled is yet unsupported');
      end
      FishVideoHandler_(self.id, 'set', 'scaled',value);
    end
    
    function value=get.useScaled(self)
      value = FishVideoHandler_(self.id, 'get','scaled');
    end
    
    
    function msec = getTimePos(self,onoff)
    % MSEC = GETTIMEPOS() gets the current time position

      msec =FishVideoHandler_(self.id, 'getTimePos')-1/self.FPS;
    end
    
    function setTimePos(self,msec)
    % GETTIMEPOS(MSEC) sets the current time position

      FishVideoHandler_(self.id, 'setTimePos',msec);
    end
    
    function frame = getCurrentFrame(self);
      frame = FishVideoHandler_(self.id,'getframe');
    end

    function value = get(self, key)
    %GET  Returns the specified VideoCapture property
    %
    %   value = cap.get(PropertyName)
    %
    % The method returns a Property value of VideoCapture.
    % PropertyName can be one of the followings:
    %
    %    'PosMsec'       Current position of the video file in milliseconds or video capture timestamp.
    %    'PosFrames'     0-based index of the frame to be decoded/captured next.
    %    'AVIRatio'      Relative position of the video file: 0 - start of the film, 1 - end of the film.
    %    'FrameWidth'    Width of the frames in the video stream.
    %    'FrameHeight'   Height of the frames in the video stream.
    %    'FPS'           Frame rate.
    %    'FourCC'        4-character code of codec.
    %    'FrameCount'    Number of frames in the video file.
    %    'Format'        Format of the Mat objects returned by retrieve() .
    %    'Mode'          Backend-specific value indicating the current capture mode.
    %    'Brightness'    Brightness of the image (only for cameras).
    %    'Contrast'      Contrast of the image (only for cameras).
    %    'Saturation'    Saturation of the image (only for cameras).
    %    'Hue'           Hue of the image (only for cameras).
    %    'Gain'          Gain of the image (only for cameras).
    %    'Exposure'      Exposure (only for cameras).
    %    'ConvertRGB'    Boolean flags indicating whether images should be converted to RGB.
    %    'Rectification' Rectification flag for stereo cameras (note: only supported by DC1394 v 2.x backend currently)
    %
    % See also cv.VideoCapture
    %
      value = FishVideoHandler_(self.id, 'capget', key);
    end
    
    function set(self, key, value)
    %SET  Sets a property in the VideoCapture.
    %
    %    cap.set(PropertyName, value)
    %
    % The method set a Property value of VideoCapture.
    % PropertyName can be one of the followings:
    %
    %    'PosMsec'       Current position of the video file in milliseconds or video capture timestamp.
    %    'PosFrames'     0-based index of the frame to be decoded/captured next.
    %    'AVIRatio'      Relative position of the video file: 0 - start of the film, 1 - end of the film.
    %    'FrameWidth'    Width of the frames in the video stream.
    %    'FrameHeight'   Height of the frames in the video stream.
    %    'FPS'           Frame rate.
    %    'FourCC'        4-character code of codec.
    %    'FrameCount'    Number of frames in the video file.
    %    'Format'        Format of the Mat objects returned by retrieve() .
    %    'Mode'          Backend-specific value indicating the current capture mode.
    %    'Brightness'    Brightness of the image (only for cameras).
    %    'Contrast'      Contrast of the image (only for cameras).
    %    'Saturation'    Saturation of the image (only for cameras).
    %    'Hue'           Hue of the image (only for cameras).
    %    'Gain'          Gain of the image (only for cameras).
    %    'Exposure'      Exposure (only for cameras).
    %    'ConvertRGB'    Boolean flags indicating whether images should be converted to RGB.
    %    'Rectification' Rectification flag for stereo cameras (note: only supported by DC1394 v 2.x backend currently)
    %
    % See also cv.VideoCapture
    %
      FishVideoHandler_(self.id, 'capset', key, value);
    end
    
    function value = get.history(self)
      value = FishVideoHandler_(self.id, 'bsget', 'History');
    end
    function set.history(self, value)
      FishVideoHandler_(self.id, 'bsset', 'History', value);
    end


% $$$     function value = get.PosMsec(self)
% $$$       value =  getTimePos(self);
% $$$     end
% $$$     function value = get.PosFrames(self)
% $$$       value =  FishVideoHandler_(self.id, 'get', 'PosFrames');
% $$$     end
% $$$     function value = get.AVIRatio(self)
% $$$       value =  FishVideoHandler_(self.id, 'get', 'AVIRatio');
% $$$     end
% $$$     function value = get.FrameWidth(self)
% $$$       value =  FishVideoHandler_(self.id, 'get', 'FrameWidth');
% $$$     end
% $$$     function value = get.FrameHeight(self)
% $$$       value =  FishVideoHandler_(self.id, 'get', 'FrameHeight');
% $$$     end
% $$$     function value = get.FPS(self)
% $$$       value =  FishVideoHandler_(self.id, 'get', 'FPS');
% $$$     end
% $$$     
% $$$     function value = get.FourCC(self)
% $$$       value =  FishVideoHandler_(self.id, 'get', 'FourCC');
% $$$     end
% $$$     function value = get.FrameCount(self)
% $$$       value =  FishVideoHandler_(self.id, 'get', 'FrameCount');
% $$$     end
% $$$     function value = get.Format(self)
% $$$       value =  FishVideoHandler_(self.id, 'get', 'Format');
% $$$     end
% $$$     function value = get.Mode(self)
% $$$       value =  FishVideoHandler_(self.id, 'get', 'Mode');
% $$$     end
% $$$     function value = get.ConvertRGB(self)
% $$$       value =  FishVideoHandler_(self.id, 'get', 'ConvertRGB');
% $$$     end
% $$$ 
% $$$     function set.PosMsec(self,value)
% $$$       setTimePos(self,value);
% $$$     end
% $$$     function set.PosFrames(self,value)
% $$$       error('not supported. Use PosMsec instead');
% $$$     end
% $$$     
    
    
% $$$     function value = get.DetectShadows(self)
% $$$       value = FishVideoHandler_(self.id, 'bsget', 'DetectShadows');
% $$$     end
% $$$     function set.DetectShadows(self, value)
% $$$       FishVideoHandler_(self.id, 'bsset', 'DetectShadows', value);
% $$$     end
% $$$ 
% $$$     function value = get.Dist2Threshold(self)
% $$$       value = FishVideoHandler_(self.id, 'bsget', 'Dist2Threshold');
% $$$     end
% $$$     function set.Dist2Threshold(self, value)
% $$$       FishVideoHandler_(self.id, 'bsset', 'Dist2Threshold', value);
% $$$     end
% $$$     function value = get.KNNSamples(self)
% $$$       value = FishVideoHandler_(self.id, 'bsget', 'kNNSamples');
% $$$     end
% $$$     function set.KNNSamples(self, value)
% $$$       FishVideoHandler_(self.id, 'bsset', 'kNNSamples', value);
% $$$     end
% $$$ 
% $$$     function value = get.NSamples(self)
% $$$       value = FishVideoHandler_(self.id, 'bsget', 'NSamples');
% $$$     end
% $$$     function set.NSamples(self, value)
% $$$       FishVideoHandler_(self.id, 'bsset', 'NSamples', value);
% $$$     end
% $$$ 
% $$$     function value = get.ShadowThreshold(self)
% $$$       value = FishVideoHandler_(self.id, 'bsget', 'ShadowThreshold');
% $$$     end
% $$$     function set.ShadowThreshold(self, value)
% $$$       FishVideoHandler_(self.id, 'bsset', 'ShadowThreshold', value);
% $$$     end
% $$$ 
% $$$     function value = get.ShadowValue(self)
% $$$       value = FishVideoHandler_(self.id, 'bsget', 'ShadowValue');
% $$$     end
% $$$     function set.ShadowValue(self, value)
% $$$       FishVideoHandler_(self.id, 'bsset', 'ShadowValue', value);
% $$$     end
  
    
    function [oimg,msk,oimg_col] = a_imrotate(self,ori,oimg,mback,omsk,oimg_col,mback_col);
    end
    
    function  [outregion] = a_computeMSERregions(self,inregion,bb2);
    end
    
    function regions = a_getRegions(self,bwimg,Iframe,rprops);
    end
    function  [boundingBox,centroid] = a_getMaxAreaRegion(self,bwimg);
    end
    function  doimg = a_interp(self,foimg,xshift,mback,type);
    end
    function  msk = a_closeHoles(self,msk);
    end
    function a_init(self);
    end
    

  end
  

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Functions for FishBlobAnalysis
  methods(Access=protected)  
    function segm = detect(self, bwimg, Iframe, Cframe)
    %dummy overload
    end
    function rp = splitRegions(self,rp,Iframe, Cframe)
    %dummy overload
    end

    function segm = getFishFeatures(self,segments);
    %dummy overload
    end
  
    
  end
  
  methods (Hidden = true)

    function release(self)
    %RELEASE  Closes video file or capturing device.
    %
    % The methods are automatically called by subsequent open() and by destructor.
    %
    % See also cv.VideoCapture.open
    %
      FishVideoHandler_(self.id, 'release');
    end
  end
  
  
end
