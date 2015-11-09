classdef FishVideoHandlerMex < handle & FishBlobAnalysis & FishVideoReader
%FISHVIDEOHANDER  wrapper class
%
% Class for video reading of the FishTracker
%
%
  
  
  
  properties (SetAccess=private)
    % Object ID
    id = [];

  end
  
  properties (Dependent)

    useScaled
    
    %ploting
    %minContourSize
    DetectShadows
    history
    nprobe
  end

  methods 
    function self = FishVideoHandlerMex(vidname,timerange,opts)
    %VIDEOCAPTURE  Create a new FishVideoHandlerMex object. VIDNAME
    %can be either a videofile or a cell like {CAMIDX, WRITEFILE}
    %for online capture from point grey devices;
    %
      
      if nargin < 1, vidname = getVideoFile(); end
      
      % make sure the image librabry of matlkab is loaded
      f = figure;
      image();
      close(f);
      vision.VideoPlayer();
      
      self@FishVideoReader(vidname);       
      self@FishBlobAnalysis(); 


      if nargin >1
        self.timeRange = timerange;
      end
      
      if nargin >2
        setOpts(self,opts);
      end
      
    end
    
    function setOpts(self,opts)

      for f1 = fieldnames(opts)'
        if any(strcmp(f1{1},{'blob','reader','detector'}))
          for f2 = fieldnames(opts.(f1{1}))'
            if isprop(self,f2{1})
              verbose('Set option "%s" to "%s"',f2{1},num2str(opts.(f1{1}).(f2{1})));
              self.(f2{1}) = opts.(f1{1}).(f2{1});
            else
              warning(sprintf('unknown property: %s',f2{1}));
            end
          end
        end
      end
    end
    
    
    
    function [seg,bwimg,frame,cframe] = step(self)
    % STEP one frame
    %
    % [segments [contours, bwimg, frame]] = vh.step();
    %
    %
% $$$       if nargout==1
% $$$         [seg] = FishVideoHandler_(self.id, 'step');
% $$$       else
% $$$         [seg,varargout{1:nargout-1}] = FishVideoHandler_(self.id, 'step');
% $$$       end

      [seg,bwimg,frame,cframe] = FishVideoHandler_(self.id, 'step');
      self.increaseCounters(); % implicit read frame, so increase
                               % the counters

      % split regions / mser not yet implemented
      self.segm = seg; % for a_getRegion dummy
      seg = self.stepBlob(bwimg,frame,cframe);
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
      frame = FishVideoHandler_(self.id,'getFrame');
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
      value = FishVideoHandler_(self.id, 'get', 'History');
    end
    function set.history(self, value)
      FishVideoHandler_(self.id, 'set', 'History', value);
    end

    function value = get.DetectShadows(self)
      value = FishVideoHandler_(self.id, 'get', 'DetectShadows');
    end
    function set.DetectShadows(self, value)
      FishVideoHandler_(self.id, 'set', 'DetectShadows', value);
    end
    
    function value = get.nprobe(self)
      value = FishVideoHandler_(self.id, 'get', 'nprobe');
    end
    function set.nprobe(self, value)
      FishVideoHandler_(self.id, 'set', 'nprobe', value);
    end
    
    
    % OVERLOAD FEATURE DETECTION (COMMENT OUT TO GET MSER)
    function segm = detect(self, bwimg, Iframe, Cframe)
      
    % get the spots from the binary image
      rp = self.a_getRegions(bwimg,Iframe,[self.rprops,{'Image'}]);
      
      for i = 1:length(rp)
        rp(i).FishFeature = rp(i).FishFeatureCRemap;
      end
      
      segm = rp;
    end
    

  end
  
  %%%%%%%%%%%% Overloaded methods (READER) 
  methods(Access=protected);
    
    function a_startReader(self);

      if isempty(self.id)
        if iscell(self.videoFile)
          camidx = self.videoFile{1};
          self.id = FishVideoHandler_('camera',camidx,self.videoFile{2});      
        else
          self.id = FishVideoHandler_(self.videoFile);      
        end
        FishVideoHandler_(self.id,'start');      
        self.reader = self.id;
      end
      
    end
    
    
    function dur = a_getDuration(self);
      dur = self.a_getNFrames()/self.a_getFrameRate();
    end
    
    function nframes = a_getNFrames(self);
      nframes = FishVideoHandler_(self.id,'get','FrameCount');
      if nframes==0
        % strange. cannot find how many frames. Assume Inf. Maybe camera
        nframes = Inf;
      end
      
    end
    
    function frameRate = a_getFrameRate(self);
      frameRate = FishVideoHandler_(self.id,'get','FPS');

      if isnan(frameRate) % strange thing that this sometimes happens
        f = VideoReader(self.videoFile);
        frameRate = f.FrameRate;
        FishVideoHandler_(self.id,'set','FPS',frameRate);
        clear('f')
      end
      
    end
    
    function frameSize= a_getFrameSize(self);
      
      frameSize = [ FishVideoHandler_(self.id,'get','FrameHeight'),...
                    FishVideoHandler_(self.id,'get','FrameWidth')];
    end
    
    function bool = a_hasFrame(self);
      nextFrame = FishVideoHandler_(self.id,'get','PosFrames');
      bool = nextFrame<self.nFrames; % starts from 0!
    end
    
    
    function a_delete(self); %% overwrite VideoCapture.delete: do nothing
      FishVideoHandler_(self.id,'delete');
    end
    
    
    function a_setCurrentTime(self,time);
    % time is given in seconds
      nextFrame = floor(time*self.frameRate);
      pos = FishVideoHandler_(self.id,'get','PosFrames');
      if nextFrame~=pos 
        FishVideoHandler_(self.id,'set','PosFrames',nextFrame);
      end
    end

    
    function  [frame,oframe] = a_readUFrame(self);      


      if self.useScaled
        self.useScaled = false;
      end
      
      if self.originalif
        [~,~,frame,oframe] = self.step();
      else
        [~,~,frame] = self.step();
        oframe = [];
      end

      
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
      if ~self.useScaled 
        self.useScaled = true;
        FishVideoHandler_(self.id, 'setScale',self.scale(1),self.scale(2),self.scalke(3));
        FishVideoHandler_(self.id, 'set','delta',sum(self.delta));
      end
      
      if self.originalif
        [~,~,frame,oframe] = self.step();
      else
        [~,~,frame] = self.step();
        oframe = [];
      end

    end

    function  [frame,oframe] = a_readScaledSFrame(self);      
      [frame,oframe] = a_readScaledUFrame(self);
      frame = single(frame)/255.;
    end

% $$$ 
% $$$     function [oimg,msk,oimg_col] = a_imrotate(self,ori,oimg,mback,omsk,oimg_col,mback_col);
% $$$     end
% $$$     
% $$$     function  [outregion] = a_computeMSERregions(self,inregion,bb2);
% $$$     end
    
    function regions = a_getRegions(self,bwimg,Iframe,rprops);
    % already done. 
      regions = self.segm;
    end
% $$$     function  [boundingBox,centroid] = a_getMaxAreaRegion(self,bwimg);
% $$$     end
% $$$     function  doimg = a_interp(self,foimg,xshift,mback,type);
% $$$     end
% $$$     function  msk = a_closeHoles(self,msk);
% $$$     end
    
    function a_init(self);
    % pass all  the properties to the C-core
          
      a_init@FishBlobAnalysis(self);
      
      FishVideoHandler_(self.id, 'set','minArea',self.minArea);
      FishVideoHandler_(self.id, 'set','maxArea',self.maxArea);
      FishVideoHandler_(self.id, 'set','maxExtent',self.maxextent);
      FishVideoHandler_(self.id, 'set','minExtent',self.minextent);
      FishVideoHandler_(self.id, 'set','colorfeature',self.originalif);
      FishVideoHandler_(self.id, 'set','minWidth',self.minWidth);
      FishVideoHandler_(self.id, 'set','featureheight',self.featureheight);
      FishVideoHandler_(self.id, 'set','featurewidth',self.featurewidth);
    
    end
    

  end
  

% $$$   
% $$$   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$   
% $$$   methods (Hidden = true)
% $$$ 
% $$$     function release(self)
% $$$     %RELEASE  Closes video file or capturing device.
% $$$     %
% $$$     % The methods are automatically called by subsequent open() and by destructor.
% $$$     %
% $$$     % See also cv.VideoCapture.open
% $$$     %
% $$$       FishVideoHandler_(self.id, 'release');
% $$$     end
% $$$   end
  
  
end
