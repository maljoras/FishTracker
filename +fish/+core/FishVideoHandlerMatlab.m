classdef FishVideoHandlerMatlab < handle & fish.core.FishVideoReaderMatlab & fish.core.FishBlobAnalysisMatlab
%FISHVIDEOHANDER  wrapper class
%
% Class for video reading of the fish.Tracker
%
%

  properties
    computeSegments = true; 
    resizeif = 0; 
    resizescale = 1;
    
    knnMethod = false;% NOT SUPORTED
    fixedSize = 0;  % NOT SUPORTED WITHOUT MEX
    difffeature = false; % not supported
  end

  properties (SetAccess = private)
    detector;
  end
  
  properties (Dependent)
    history
  end

  methods
    
    function self = FishVideoHandlerMatlab(vidname,timerange,useKnn,opts)
    %VIDEOCAPTURE  Create a new FishVideoHandler object
    %
      if useKnn
        error('knn method not supported in Matlab');
      end
      
      self@fish.core.FishBlobAnalysisMatlab(); 
      self@fish.core.FishVideoReaderMatlab(vidname,timerange);  %%% SOMEHOW
                                                  %%% MATALB READER
                                                  %%% DOES
                                                  %%% NOT WORK ?!?
      self.detector = fish.core.FishForegroundDetectorMatlab();  

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
        fish.helper.verbose('WARNING: Nothing set in VideoHandler');
      end

      self.frameFormat = [self.grayFormat,self.detector.expectedFrameFormat];
    end
    
    function plotting(self,bool);
    % plotting not implemented...
    end

    
    function [segm,timeStamp,frame,varargout] = step(self)
    % STEP one frame
    %
    % [segments,timeStamp, frame, [oframe]] = vh.step();
    %
    %

      if nargout>3 || self.colorfeature
        self.originalif = true; 
        [frame,oframe] = self.readFrame();            
        if self.resizeif
          oframe = imresize(oframe,self.resizescale);
          frame = imresize(frame,self.resizescale);
        end
        varargout{1} = oframe;
      else
        self.originalif = false;
        frame = self.readFrame();
        if self.resizeif
          frame = imresize(frame,self.resizescale);
        end
        oframe = frame;
      end
      timeStamp = self.currentTime;
      
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


  end
end
