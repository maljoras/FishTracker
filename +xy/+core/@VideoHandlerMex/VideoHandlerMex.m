classdef VideoHandlerMex < handle & xy.core.BlobAnalysis & xy.core.VideoReader
%FISHVIDEOHANDLER  wrapper class
%
% Class for video handling of the xy.Tracker. Inherits from
% BlobAnalisys and VideoReader, since both video
% capturing/reading and blob analysis is done in a multi-threading
% C++ code. 
%
% NOTE: if "unknown error" occurs in MEX file then try to set in
% linux terminal:
%   sudo echo "vm.max_map_count=16777216" | tee -a /etc/sysctl.conf
%   sudo sysctl -w vm.max_map_count=16777216

  
  properties (SetAccess=private)
    id = [];
    knnMethod = true;
    useScaled = false;
  end
  
  properties (Dependent)
    DetectShadows
    history
    nprobe

    computeSegments
    resizeif
    resizescale
    fixedSize 
    difffeature
    
    inverted
    nskip
    adjustThresScale
  end

  properties 
    codec = 'X264';
  end

  
  methods(Static)
    function bool = installed()
      bool = ~~exist('xyVideoHandler_');
    end
  end
    
  
  methods 
    
    
    function self = VideoHandlerMex(vidname,timerange,knnMethod,opts)
    %VIDEOCAPTURE  Create a new VideoHandlerMex object. VIDNAME
    %can be either a videofile or a cell like {CAMIDX, WRITEFILE}
    %for online capture from point grey devices;
    %
      
      if nargin < 1, vidname = xy.helper.getVideoFile(); end

      % make sure the image library of matlab is loaded (seems not
      % necessary anymore, was a bug in FlyCapture SDK)
      %f = figure;
      %image();
      %close(f);
      %vision.VideoPlayer();
      global global_knnMethod; % somewaht a hack... since foreground detector is now
                               % implicitely included in the reader...     
      if exist('knnMethod','var') && ~isempty(knnMethod)
        global_knnMethod = knnMethod;
      else
        global_knnMethod = true;
      end
      if ~exist('xyVideoHandler_')
        error('Cannot find mex file. Forgot to run make ?');
      end


      self@xy.core.VideoReader(vidname);       
      self@xy.core.BlobAnalysis(opts);       
      if nargin >1
        self.timeRange = timerange;
      end
      
      self.frameFormat = 'RGBU'; % default
      self.knnMethod = global_knnMethod;
      clear global global_knnMethod;
      

      if nargin >3
        setOpts(self,opts);
      else
        self.initialize(); % initialize blob
      end
      
    end
     
    
    function release(self)
    %RELEASE  Closes video file or capturing device.
    %
    % The methods are automatically called by subsequent open() and by destructor.
    %
    %
      xyVideoHandler_(self.id, 'stop');
    end

    
    
    function setOpts(self,opts)
      
      for f1 = fieldnames(opts)'
        if any(strcmp(f1{1},{'blob','reader','detector'}))
          for f2 = fieldnames(opts.(f1{1}))'
            if isprop(self,f2{1})
              if any(self.(f2{1}) ~= opts.(f1{1}).(f2{1})) 
                xy.helper.verbose('Set "%s.%s" to "%s"',f1{1},f2{1},num2str(opts.(f1{1}).(f2{1})));
                self.(f2{1}) = opts.(f1{1}).(f2{1});
              end
            else
              warning(sprintf('unknown property: %s',f2{1}));
            end
          end
        end
      end
      self.initialize(0); % update Blob
    end
    
    
    
    function [seg,timeStamp,frame] = step(self)
    % STEP one frame
    %
    % [segments, timestamp, [ frame]] = vh.step();
    %
    %
      
      if nargout<3 
        [seg,timeStamp] = xyVideoHandler_(self.id, 'step');
      else
        [seg,timeStamp,frame] = xyVideoHandler_(self.id, 'step');
      end
      if timeStamp>=0
        self.increaseCounters(timeStamp); % implicit read frame, so increase
                                          % the counters
      else
        self.increaseCounters(Inf); 
      end
      
% $$$       if self.computeSegments
% $$$         self.segm = seg;
% $$$         bwimg = self.getCurrentBWImg(); 
% $$$         % split regions / mser not yet implemented
% $$$         seg = self.stepBlob(bwimg,frame,[]);
% $$$       else
% $$$         seg = [];
% $$$       end
      
    end
  
    
    function set.computeSegments(self,value)
      xyVideoHandler_(self.id, 'set', 'computeSegments',value);
    end
    
    function value = get.computeSegments(self)
      value =xyVideoHandler_(self.id, 'get', 'computeSegments');
    end
    
    function value = get.resizeif(self)
      value = xyVideoHandler_(self.id, 'get', 'resizeif');
    end
    function set.resizeif(self,value)
      xyVideoHandler_(self.id, 'set', 'resizeif',value);
    end
    
    function value = get.difffeature(self)
      value = xyVideoHandler_(self.id, 'get', 'difffeature');
    end
    function set.difffeature(self,value)
      xyVideoHandler_(self.id, 'set', 'difffeature',value);
    end

    function value = get.fixedSize(self)
      value = xyVideoHandler_(self.id, 'get', 'fixedSize');
    end
    
    function set.fixedSize(self,value)
      xyVideoHandler_(self.id, 'set', 'fixedSize',value);
    end

    function value = get.resizescale(self)
      value = xyVideoHandler_(self.id, 'get', 'resizescale');
    end
    
    function set.resizescale(self,value)
      xyVideoHandler_(self.id, 'set', 'resizescale',value);
    end
    
    
    function setToScaledFormat(self,scale,delta)
      
      if nargin<2
        scale = self.scale;
      end
      if nargin<3
        delta = self.delta;
      end
      setToScaledFormat@xy.core.VideoReader(self,scale,delta);
      
      xyVideoHandler_(self.id, 'setScale', self.scale(1),self.scale(2),self.scale(3));
      xyVideoHandler_(self.id, 'set','delta', sum(self.delta));
      xyVideoHandler_(self.id, 'set', 'scaled',true);
      self.useScaled = true;
    end
    
    function setToRGBFormat(self)
      setToRGBFormat@xy.core.VideoReader(self);
      self.useScaled = false;
      xyVideoHandler_(self.id, 'set', 'scaled',false);
    end
    
    function frame = getCurrentFrame(self)
      frame = xyVideoHandler_(self.id,'getFrame');
    end
    
    function bwimg = getCurrentBWImg(self)
      bwimg = xyVideoHandler_(self.id,'getBWImg');
    end

    function value = get(self, key)
      value = xyVideoHandler_(self.id, 'capget', key);
    end
    
    function set(self, key, value)
      xyVideoHandler_(self.id, 'set', key, value);
    end
    
    function value = get.history(self)
      value = xyVideoHandler_(self.id, 'get', 'History');
    end
    function set.history(self, value)
      xyVideoHandler_(self.id, 'set', 'History', value);
    end

    function value = get.DetectShadows(self)
      value = xyVideoHandler_(self.id, 'get', 'DetectShadows');
    end
    function set.DetectShadows(self, value)
      xyVideoHandler_(self.id, 'set', 'DetectShadows', value);
    end
    
    function value = get.nprobe(self)
      value = xyVideoHandler_(self.id, 'get', 'nprobe');
    end
    function set.nprobe(self, value)
      xyVideoHandler_(self.id, 'set', 'nprobe', value);
    end
    
    function value = get.nskip(self)
      value = xyVideoHandler_(self.id, 'get', 'nskip');
    end
    function set.nskip(self, value)
      xyVideoHandler_(self.id, 'set', 'nskip', value);
    end
    function value = get.adjustThresScale(self)
      value = xyVideoHandler_(self.id, 'get', 'adjustThresScale');
    end
    function set.adjustThresScale(self, value)
      xyVideoHandler_(self.id, 'set', 'adjustThresScale', value);
    end

    function value = get.knnMethod(self)
      value = xyVideoHandler_(self.id, 'get', 'knnMethod');
    end
    
    function value = get.inverted(self)
      value = xyVideoHandler_(self.id, 'get', 'inverted');
    end
    function set.inverted(self, value)
      xyVideoHandler_(self.id, 'set', 'inverted', value);
    end

    
    % OVERLOAD FEATURE DETECTION (COMMENT OUT TO GET MSER)
    function segm = detect(self, bwimg, Iframe, Cframe)
      
    % get the spots from the binary image
      rp = self.a_getRegions(bwimg,Iframe,[self.rprops,{'Image'}]);

% $$$ 
% $$$       for i = 1:length(rp)
% $$$ % $$$         img = rp(i).IdentityFeatureCRemap;
% $$$ % $$$         sz = size(img);
% $$$ % $$$         sz(1) = ceil(sz(1)/3);
% $$$ % $$$         sz(2) = ceil(sz(2)/3);
% $$$ % $$$         
% $$$ % $$$         dimg = zeros(sz);      
% $$$ % $$$         for j = 1:size(img,3)
% $$$ % $$$           tmp = dct2(img(:,:,j));
% $$$ % $$$           dimg(:,:,j) = tmp(1:sz(1),1:sz(2));
% $$$ % $$$         end
% $$$ % $$$         
% $$$ % $$$         rp(i).IdentityFeature = dimg;
% $$$          rp(i).IdentityFeature = rp(i).IdentityFeatureCRemap;
% $$$       end
      
      segm = rp;
    end
    
    
    function resetBkg(self)
      xyVideoHandler_(self.id, 'resetBkg');
    end
    
    function plotting(self,bool)
      xyVideoHandler_(self.id,'set','plotif',bool);      
    end
    
    function bool = isGrabbing(self)
      bool = xyVideoHandler_(self.id,'get','camera') ; 
    end

    function verbose(self)
      verbose@xy.core.VideoReader(self);
      if self.knnMethod 
        xy.helper.verbose('Using KNN for foreground subtraction...');
      else
        xy.helper.verbose('Using Thresholder for foreground subtraction...');
      end
    end
    
  end
  
  %%%%%%%%%%%% Overloaded methods (READER) 
  methods(Access=protected);
    
    function a_startReader(self)

      global global_knnMethod; % somewaht a hack...
      knnMethod = global_knnMethod;
      if isempty(self.id)
        if iscell(self.videoFile)
          camidx = self.videoFile{1};
          if exist(self.videoFile{2},'file')
            fprintf(['\n\nREALLY OVERWRITE EXISTING VIDEO FILE ? CTRL-C ' ...
                     'to abort..ENTER to continue:']);
            pause;
          end
          
          self.id = xyVideoHandler_('camera',camidx,self.videoFile{2},upper(self.codec),knnMethod);      
        else
          self.id = xyVideoHandler_(self.videoFile,knnMethod);      
        end
        xyVideoHandler_(self.id,'start');      
        self.reader = self.id;
      end
    end
    

    
    function dur = a_getDuration(self)
      dur = self.a_getNFrames()/self.a_getFrameRate();
    end
    
    function nframes = a_getNFrames(self)
      nframes = xyVideoHandler_(self.id,'get','FrameCount');
      
      if nframes<=0
        % strange. cannot find how many frames. Assume Inf. Maybe camera
        nframes = Inf;
      end
      
    end
    
    function frameRate = a_getFrameRate(self)
      frameRate = xyVideoHandler_(self.id,'get','FPS');

      if isnan(frameRate) % strange thing that this sometimes happens
        f = VideoReader(self.videoFile);
        frameRate = f.FrameRate;
        xyVideoHandler_(self.id,'set','FPS',frameRate);
        clear('f')
      end
      
    end
    
    function frameSize= a_getFrameSize(self)
      
      frameSize = [ xyVideoHandler_(self.id,'get','FrameHeight'),...
                    xyVideoHandler_(self.id,'get','FrameWidth')];
    end
    
    function bool = a_hasFrame(self)
      nextFrame = xyVideoHandler_(self.id,'get','PosFrames');
      bool = nextFrame<self.nFrames && ~isinf(self.currentTime); % starts from 0!
    end
    
    
    function a_delete(self) %% overwrite VideoCapture.delete: do nothing
      xyVideoHandler_(self.id,'delete');
    end
    
    
    function a_setCurrentTime(self,time)
    % time is given in seconds
      nextFrame = max(floor(time*self.frameRate),0);
      pos = xyVideoHandler_(self.id,'get','PosFrames');
      if nextFrame~=pos 
        xyVideoHandler_(self.id,'set','PosFrames',nextFrame);
      end
    end

    
    function  [frame,oframe] = a_readUFrame(self)


      if self.useScaled
        xy.helper.verbose('WARNING: Set permanently to non-scaled format')
        self.setToRGBFormat();
      end  
      [~,~,frame] = self.step();
      oframe = frame;
    end
    

    function [frame,oframe] = a_readSFrame(self)
      [frame,oframe] = a_readUFrame(self);
      frame = single(frame)/255.;
      oframe = single(oframe)/255.;
    end

    function  [frame,oframe] = a_readGrayUFrame(self) 
      [frame,oframe] = a_readUFrame(self); 
    end


    function [frame,oframe] = a_readGraySFrame(self)
      [frame,oframe] = a_readSFrame(self);
    end
    
    function [frame,oframe] = a_readInvertedGrayUFrame(self)
      [frame,oframe] = a_readGrayUFrame(self); 
      frame = 255-frame;
    end
    
    function [frame,oframe] = a_readInvertedGraySFrame(self)
      [frame,oframe] = a_readGraySFrame(self); 
      frame = 1-frame;
    end
    
    function  [frame,oframe] = a_readScaledUFrame(self)     

      if ~self.useScaled
        xy.helper.verbose('WARNING: Set permanently to scaled format')
        self.setToScaledFormat();
      end
      
      [~,~,frame] = self.step();
      oframe = frame;
      
    end

    function  [frame,oframe] = a_readScaledSFrame(self)   
      [frame,oframe] = a_readScaledUFrame(self);
      frame = single(frame)/255.;
    end

% $$$ 
% $$$     function [oimg,msk,oimg_col] = a_imrotate(self,ori,oimg,mback,omsk,oimg_col,mback_col);
% $$$     end
% $$$     
% $$$     function  [outregion] = a_computeMSERregions(self,inregion,bb2);
% $$$     end
    
    function regions = a_getRegions(self,bwimg,Iframe,rprops)
    % already done. 
      regions = self.segm;
    end
% $$$     function  [boundingBox,centroid] = a_getMaxAreaRegion(self,bwimg);
% $$$     end
% $$$     function  doimg = a_interp(self,foimg,xshift,mback,type);
% $$$     end
% $$$     function  msk = a_closeHoles(self,msk);
% $$$     end
    
    function a_init(self)
    % pass all  the properties to the C-core
          
    %a_init@xy.core.BlobAnalysis(self);
      
      xyVideoHandler_(self.id, 'set','minArea',self.minArea);
      xyVideoHandler_(self.id, 'set','maxArea',self.maxArea);
      xyVideoHandler_(self.id, 'set','maxExtent',self.maxextent);
      xyVideoHandler_(self.id, 'set','minExtent',self.minextent);
      xyVideoHandler_(self.id, 'set','colorfeature',self.colorfeature || self.originalif);
      xyVideoHandler_(self.id, 'set','minWidth',self.minWidth);
      xyVideoHandler_(self.id, 'set','featureheight',self.featureheight);
      xyVideoHandler_(self.id, 'set','featurewidth',self.featurewidth);

    end
    

  end
  

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  
end
