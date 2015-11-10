classdef FishVideoHandler < handle & FishVideoReader & FishBlobAnalysis
%FISHVIDEOHANDER  wrapper class
%
% Class for video reading of the FishTracker
%
%
  
  properties
    computeSegments = true;
  end
  
  
  properties (SetAccess = private)
    detector;
  end
  
  properties (Dependent)
    history

  end

  methods
    
    function self = FishVideoHandler(vidname,timerange,opts)
    %VIDEOCAPTURE  Create a new FishVideoHandler object
    %
      if iscell(vidname)
        error('Capture not supported');
        vidname = vidname{1}; % capture device might be suported by
                              % open cv
      end
      
      self@FishBlobAnalysis(); 
      self@FishVideoReader(vidname,timerange); 
      self.detector = FishForegroundDetector();  

      if exist('opts','var')
        self.setOpts(opts);
      end
      
    end
    
    
    function self = setOpts(self,opts)
    % to set the OPTIONS use the keywords "blob" "reader" and "detector"
      setif = 0;
      if isfield(opts,'blob')
        for f = fieldnames(opts.blob)'
          if isprop(self,f{1})
            self.(f{1}) = opts.blob.(f{1});
            setif = 1;
            verbose('Set BLOB option "%s" to "%s"',f{1},num2str(opts.blob.(f{1})));
          end
        end
      end
      
      if isfield(opts,'reader')
        for f = fieldnames(opts.reader)'
          if isprop(self,f{1})
            self.(f{1}) = opts.reader.(f{1});
            setif = 1;
            verbose('Set READER option "%s" to "%s"',f{1},num2str(opts.reader.(f{1})));
          end
        end
      end

      if isfield(opts,'detector')
        for f = fieldnames(opts.detector)'
          if isprop(self.detector,f{1})
            self.detector.(f{1}) = opts.detector.(f{1});
            setif = 1;
            verbose('Set DETECTOR option "%s" to "%s"',f{1},num2str(opts.detector.(f{1})));
          end
        end
      end
      
      if ~setif
        warning('Nothing set in VideoHandler');
      end

      self.frameFormat = [self.grayFormat,self.detector.expectedFrameFormat];
    end

    
    function [segm,bwmsk,frame,varargout] = step(self)
    % STEP one frame
    %
    % [segments,bwimg, frame, [oframe]] = vh.step();
    %
    %

      if nargout>3 || self.colorfeature
        self.originalif = true; 
        [frame,oframe] = self.readFrame();            
        varargout{1} = oframe;
      else
        self.originalif = false;
        frame = self.readFrame();
        oframe = frame;
      end
      bwmsk = self.detector.step(frame);

      if self.computeSegments
        segm = self.stepBlob(bwmsk,frame,oframe);
      else
        segm = [];
      end
      
    end        

    %% detector methods
    function set.history(self,value)
      self.detector.history = value;
    end
    function value = get.history(self)
      value = self.detector.history;
    end
    
    function resetBkg(self)
      self.detector.reset();
    end
    
    function plotting(self,bool);
    % plotting not implemented...
    end
    

  end
end
