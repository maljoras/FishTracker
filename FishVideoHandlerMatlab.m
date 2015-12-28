classdef FishVideoHandlerMatlab < handle & FishVideoReader & FishBlobAnalysisMatlab
%FISHVIDEOHANDER  wrapper class
%
% Class for video reading of the FishTracker
%
%

  properties
    computeSegments = true; % NOT SUPORTED
    resizeif = 0; % % NOT SUPORTED
    resizescale = 1; % NOT SUPORTED
    knnMethod = false;% NOT SUPORTED
    fixedSize = 0;  % NOT SUPORTED WITHOUT MEX
    difffeature = false; % no supported
  end

  properties (SetAccess = private)
    detector;
  end
  
  properties (Dependent)
    history
  end

  methods
    
    function self = FishVideoHandlerMatlab(vidname,timerange,opts)
    %VIDEOCAPTURE  Create a new FishVideoHandler object
    %
      self@FishBlobAnalysisMatlab(); 
      self@FishVideoReader(vidname,timerange);  %%% SOMEHOW
                                                  %%% MATALB READER
                                                  %%% DOES
                                                  %%% NOT WORK ?!?
      self.detector = FishForegroundDetectorMatlab();  

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
          end
        end
      end
      
      if isfield(opts,'reader')
        for f = fieldnames(opts.reader)'
          if isprop(self,f{1})
            self.(f{1}) = opts.reader.(f{1});
            setif = 1;
          end
        end
      end

      if isfield(opts,'detector')
        for f = fieldnames(opts.detector)'
          if isprop(self.detector,f{1})
            self.detector.(f{1}) = opts.detector.(f{1});
            setif = 1;
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
      segm = self.stepBlob(bwmsk,frame,oframe);
    end        

    %% detector methods
    function set.history(self,value)
      self.detector.history = value;
    end
    function value = get.history(self)
      value = self.detector.history;
    end


  end
end
