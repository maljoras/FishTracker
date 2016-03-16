classdef FishVideoHandlerMex < handle & fish.core.FishBlobAnalysis & fish.core.FishVideoReader
%FISHVIDEOHANDLER  wrapper class
%
% Class for video handling of the fish.Tracker. Inherits from
% FishBlobAnalisys and FishVideoReader, since both video
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

  methods(Static)
    function bool = installed();
      bool = ~~exist('FishVideoHandler_');
    end
  end
    
  
  methods 
    
    function self = FishVideoHandlerMex(vidname,timerange,knnMethod,opts)
    %VIDEOCAPTURE  Create a new FishVideoHandlerMex object. VIDNAME
    %can be either a videofile or a cell like {CAMIDX, WRITEFILE}
    %for online capture from point grey devices;
    %
      
      if nargin < 1, vidname = fish.helper.getVideoFile(); end

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
      if ~exist('FishVideoHandler_')
        error('Cannot find mex file. Forgot to run make ?');
      end
      self@fish.core.FishVideoReader(vidname);       
      self@fish.core.FishBlobAnalysis();       
      
      self.frameFormat = 'RGBU'; % default
      self.knnMethod = global_knnMethod;
      clear global global_knnMethod;
      
      if nargin >1
        self.timeRange = timerange;
      end

      if nargin >3
        setOpts(self,opts);
      else
        self.initialize(); % initialize blob
      end
      
    end
    
    function setOpts(self,opts)
      
      for f1 = fieldnames(opts)'
        if any(strcmp(f1{1},{'blob','reader','detector'}))
          for f2 = fieldnames(opts.(f1{1}))'
            if isprop(self,f2{1})
              if any(self.(f2{1}) ~= opts.(f1{1}).(f2{1})) 
                fish.helper.verbose('Set "%s.%s" to "%s"',f1{1},f2{1},num2str(opts.(f1{1}).(f2{1})));
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
    
    
    
    function [seg,frame] = step(self)
    % STEP one frame
    %
    % [segments [ frame]] = vh.step();
    %
    %
      
      if nargout<2 
        [seg,timeStamp] = FishVideoHandler_(self.id, 'step');
      else
        [seg,timeStamp,frame] = FishVideoHandler_(self.id, 'step');
      end
      
      self.increaseCounters(timeStamp); % implicit read frame, so increase
                                        % the counters
            

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
      FishVideoHandler_(self.id, 'set', 'computeSegments',value);
    end
    
    function value = get.computeSegments(self)
      value =FishVideoHandler_(self.id, 'get', 'computeSegments');
    end
    
    function value = get.resizeif(self);
      value = FishVideoHandler_(self.id, 'get', 'resizeif');
    end
    function set.resizeif(self,value);
      FishVideoHandler_(self.id, 'set', 'resizeif',value);
    end
    
    function value = get.difffeature(self);
      value = FishVideoHandler_(self.id, 'get', 'difffeature');
    end
    function set.difffeature(self,value);
      FishVideoHandler_(self.id, 'set', 'difffeature',value);
    end

    function value = get.fixedSize(self);
      value = FishVideoHandler_(self.id, 'get', 'fixedSize');
    end
    
    function set.fixedSize(self,value);
      FishVideoHandler_(self.id, 'set', 'fixedSize',value);
    end


    function value = get.resizescale(self);
      value = FishVideoHandler_(self.id, 'get', 'resizescale');
    end
    
    function set.resizescale(self,value);
      FishVideoHandler_(self.id, 'set', 'resizescale',value);
    end
    
    
    function setToScaledFormat(self,scale,delta)
      
      if nargin<2
        scale = self.scale;
      end
      if nargin<3
        delta = self.delta;
      end
      setToScaledFormat@fish.core.FishVideoReader(self,scale,delta);
      
      FishVideoHandler_(self.id, 'setScale', self.scale(1),self.scale(2),self.scale(3));
      FishVideoHandler_(self.id, 'set','delta', sum(self.delta));
      FishVideoHandler_(self.id, 'set', 'scaled',true);
      self.useScaled = true;
    end
    
    function setToRGBFormat(self)
      setToRGBFormat@fish.core.FishVideoReader(self);
      self.useScaled = false;
      FishVideoHandler_(self.id, 'set', 'scaled',false);
    end
    
    function frame = getCurrentFrame(self);
      frame = FishVideoHandler_(self.id,'getFrame');
    end
    
    function bwimg = getCurrentBWImg(self);
      bwimg = FishVideoHandler_(self.id,'getBWImg');
    end

    function value = get(self, key)
      value = FishVideoHandler_(self.id, 'capget', key);
    end
    
    function set(self, key, value)
      FishVideoHandler_(self.id, 'set', key, value);
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
    
    function value = get.nskip(self)
      value = FishVideoHandler_(self.id, 'get', 'nskip');
    end
    function set.nskip(self, value)
      FishVideoHandler_(self.id, 'set', 'nskip', value);
    end
    function value = get.adjustThresScale(self)
      value = FishVideoHandler_(self.id, 'get', 'adjustThresScale');
    end
    function set.adjustThresScale(self, value)
      FishVideoHandler_(self.id, 'set', 'adjustThresScale', value);
    end

    function value = get.knnMethod(self)
      value = FishVideoHandler_(self.id, 'get', 'knnMethod');
    end
    
    function value = get.inverted(self)
      value = FishVideoHandler_(self.id, 'get', 'inverted');
    end
    function set.inverted(self, value)
      FishVideoHandler_(self.id, 'set', 'inverted', value);
    end

    
    % OVERLOAD FEATURE DETECTION (COMMENT OUT TO GET MSER)
    function segm = detect(self, bwimg, Iframe, Cframe)
      
    % get the spots from the binary image
      rp = self.a_getRegions(bwimg,Iframe,[self.rprops,{'Image'}]);

% $$$ 
% $$$       for i = 1:length(rp)
% $$$ % $$$         img = rp(i).FishFeatureCRemap;
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
% $$$ % $$$         rp(i).FishFeature = dimg;
% $$$          rp(i).FishFeature = rp(i).FishFeatureCRemap;
% $$$       end
      
      segm = rp;
    end
    
    
    function resetBkg(self)
      FishVideoHandler_(self.id, 'resetBkg');
    end
    
    function plotting(self,bool)
      FishVideoHandler_(self.id,'set','plotif',bool);      
    end
    
    function bool = isGrabbing(self);
      bool = FishVideoHandler_(self.id,'get','camera') ; 
    end

    function verbose(self);
      verbose@fish.core.FishVideoReader(self);
      if self.knnMethod 
        fish.helper.verbose('Using KNN for foreground subtraction...');
      else
        fish.helper.verbose('Using Thresholder for foreground subtraction...');
      end
    end
    
  end
  
  %%%%%%%%%%%% Overloaded methods (READER) 
  methods(Access=protected);
    
    function a_startReader(self);

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
          
          self.id = FishVideoHandler_('camera',camidx,self.videoFile{2},knnMethod);      
        else
          self.id = FishVideoHandler_(self.videoFile,knnMethod);      
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
      nextFrame = max(floor(time*self.frameRate),0);
      pos = FishVideoHandler_(self.id,'get','PosFrames');
      if nextFrame~=pos 
        FishVideoHandler_(self.id,'set','PosFrames',nextFrame);
      end
    end

    
    function  [frame,oframe] = a_readUFrame(self);      


      if self.useScaled
        warning('Set permanently to non-scaled format')
        self.setToRGBFormat();
      end  
      [~,~,frame] = self.step();
      oframe = frame;
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

      if ~self.useScaled
        warning('Set permanently to scaled format')
        self.setToScaledFormat();
      end
      
      [~,~,frame] = self.step();
      oframe = frame;
      
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
          
      %a_init@fish.core.FishBlobAnalysis(self);
      
      FishVideoHandler_(self.id, 'set','minArea',self.minArea);
      FishVideoHandler_(self.id, 'set','maxArea',self.maxArea);
      FishVideoHandler_(self.id, 'set','maxExtent',self.maxextent);
      FishVideoHandler_(self.id, 'set','minExtent',self.minextent);
      FishVideoHandler_(self.id, 'set','colorfeature',self.colorfeature || self.originalif);
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
