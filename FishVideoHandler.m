classdef FishVideoHandler < handle
    %FISHVIDEOHANDER  wrapper class
    %
    % Class for video reading of the FishTracker
    %
    %
    
    properties (SetAccess = private)
        % Object ID
        id;
        fname;
    end
    properties (Dependent)
      Plotting
      UseScaled
      
      MinContourSize
      
      PosMsec      % Current position of the video file in milliseconds or video capture timestamp.
      PosFrames    % 0-based index of the frame to be decoded/captured next.
      AVIRatio     % Relative position of the video file: 0 - start of the film, 1 - end of the film.
      FrameWidth   % Width of the frames in the video stream.
      FrameHeight  % Height of the frames in the video stream.
      FPS          % Frame rate.
      FourCC       %  4-character code of codec.
      FrameCount   %  Number of frames in the video file.
      Format       %  Format of the Mat objects returned by retrieve() .
      ConvertRGB   %  Boolean flags indicating whether images should be converted to RGB.
      Rectification % Rectification flag for stereo cameras (note: only supported by DC1394 v 2.x backend currently)
      Mode            %   Backend-specific value indicating the current capture mode.
      
      DetectShadows
      Dist2Threshold
      History
      KNNSamples
      NSamples
      ShadowThreshold
      ShadowValue
    
    end

    methods
        function this = FishVideoHandler(filename)
            %VIDEOCAPTURE  Create a new FishVideoHandler object
            %
            if nargin < 1, filename = getVideoFile(); end
            this.id = FishVideoHandler_(filename);
            this.fname = filename;
          end
        
        function delete(this)
            %DELETE  Destructor of VideoCapture object
            FishVideoHandler_(this.id, 'delete');
        end
        
        function varargout = step(this)
            % STEP one frame
            %
            % [segments [contours, bwimg, frame]] = vh.step();
            %
            %
            [varargout{1:nargout}] = FishVideoHandler_(this.id, 'step');
        end
        
        function frame = getFrame(this)
            %GETDATA gets the segments 
            %
            %    [segments [contours, bwimg, frame]] = vh.getData();
            %
            %
            frame = FishVideoHandler_(this.id, 'getframe');
          end
          
        function setScale(this,rgbchannel,delta)
            %SETSCALE(RGBSCALE,DELTA) for the scaled image. Have to turn on!
            %
            %   >> vh.setScale([3,1,2],[2,1,1]);
            %   >> vf.setScale('on');
            if ~nargin==2 && ischar(rgbchannel)
              FishVideoHandler_(this.id, 'setScale',rgbchannel);
            else
              FishVideoHandler_(this.id, 'setScale',rgbchannel(1),rgbchannel(2),rgbchannel(3),sum(delta));
            end
          end
          
          function set.Plotting(this,value)
            if ~value
              FishVideoHandler_(this.id, 'setPlotting','off');
            else
              FishVideoHandler_(this.id, 'setPlotting','on');
            end
          end
          
          function value = get.Plotting(this)
           value = FishVideoHandler_(this.id, 'isPlotting');
          end
          
          function set.UseScaled(this,value)
            if ~value
              FishVideoHandler_(this.id, 'setScaled','off');
            else
              FishVideoHandler_(this.id, 'setScaled','on');
            end
          end
          
          function value=get.UseScaled(this)
           value = FishVideoHandler_(this.id, 'isScaled');
          end
          
          
          function msec = getTimePos(this,onoff)
          % MSEC = GETTIMEPOS() gets the current time position

           msec =FishVideoHandler_(this.id, 'getTimePos')-1/this.FPS;
          end
          
          function setTimePos(this,msec)
          % GETTIMEPOS(MSEC) sets the current time position

           FishVideoHandler_(this.id, 'setTimePos',msec);
          end

          
        function value = get(this, key)
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
            value = FishVideoHandler_(this.id, 'get', key);
        end
        
        function set(this, key, value)
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
            FishVideoHandler_(this.id, 'set', key, value);
        end
        function value = get.MinContourSize(this)
          value =  FishVideoHandler_(this.id, 'getMinContourSize');
        end
        function set.MinContourSize(this,value)
          FishVideoHandler_(this.id, 'setMinContourSize',value);
        end
        
        function value = get.PosMsec(this)
          value =  getTimePos(this);
        end
        function value = get.PosFrames(this)
          value =  FishVideoHandler_(this.id, 'get', 'PosFrames');
        end
        function value = get.AVIRatio(this)
          value =  FishVideoHandler_(this.id, 'get', 'AVIRatio');
        end
        function value = get.FrameWidth(this)
          value =  FishVideoHandler_(this.id, 'get', 'FrameWidth');
        end
        function value = get.FrameHeight(this)
          value =  FishVideoHandler_(this.id, 'get', 'FrameHeight');
        end
        function value = get.FPS(this)
          value =  FishVideoHandler_(this.id, 'get', 'FPS');
        end
        
        function value = get.FourCC(this)
          value =  FishVideoHandler_(this.id, 'get', 'FourCC');
        end
        function value = get.FrameCount(this)
          value =  FishVideoHandler_(this.id, 'get', 'FrameCount');
        end
        function value = get.Format(this)
          value =  FishVideoHandler_(this.id, 'get', 'Format');
        end
        function value = get.Mode(this)
          value =  FishVideoHandler_(this.id, 'get', 'Mode');
        end
        function value = get.ConvertRGB(this)
          value =  FishVideoHandler_(this.id, 'get', 'ConvertRGB');
        end

        function set.PosMsec(this,value)
          setTimePos(this,value);
        end
        function set.PosFrames(this,value)
           error('not supported. Use PosMsec instead');
        end
        function set.AVIRatio(this,value)
           FishVideoHandler_(this.id, 'set', 'AVIRatio');
        end
        function set.FrameWidth(this,value)
           FishVideoHandler_(this.id, 'set', 'FrameWidth');
        end
        function set.FrameHeight(this,value)
           FishVideoHandler_(this.id, 'set', 'FrameHeight');
        end
        function set.FPS(this,value)
           FishVideoHandler_(this.id, 'set', 'FPS');
        end
        function set.FourCC(this,value)
           FishVideoHandler_(this.id, 'set', 'FourCC');
        end
        function set.FrameCount(this,value)
          FishVideoHandler_(this.id, 'set', 'FrameCount');
        end
        function set.Format(this,value)
          FishVideoHandler_(this.id, 'set', 'Format');
        end
        function set.Mode(this,value)
          FishVideoHandler_(this.id, 'set', 'Mode');
        end
        function set.ConvertRGB(this,value)
          FishVideoHandler_(this.id, 'set', 'ConvertRGB');
        end
        
        
        function value = get.DetectShadows(this)
            value = FishVideoHandler_(this.id, 'bsget', 'DetectShadows');
        end
        function set.DetectShadows(this, value)
            FishVideoHandler_(this.id, 'bsset', 'DetectShadows', value);
        end

        function value = get.Dist2Threshold(this)
            value = FishVideoHandler_(this.id, 'bsget', 'Dist2Threshold');
        end
        function set.Dist2Threshold(this, value)
            FishVideoHandler_(this.id, 'bsset', 'Dist2Threshold', value);
        end

        function value = get.History(this)
            value = FishVideoHandler_(this.id, 'bsget', 'History');
        end
        function set.History(this, value)
            FishVideoHandler_(this.id, 'bsset', 'History', value);
        end

        function value = get.KNNSamples(this)
            value = FishVideoHandler_(this.id, 'bsget', 'kNNSamples');
        end
        function set.KNNSamples(this, value)
            FishVideoHandler_(this.id, 'bsset', 'kNNSamples', value);
        end

        function value = get.NSamples(this)
            value = FishVideoHandler_(this.id, 'bsget', 'NSamples');
        end
        function set.NSamples(this, value)
            FishVideoHandler_(this.id, 'bsset', 'NSamples', value);
        end

        function value = get.ShadowThreshold(this)
            value = FishVideoHandler_(this.id, 'bsget', 'ShadowThreshold');
        end
        function set.ShadowThreshold(this, value)
            FishVideoHandler_(this.id, 'bsset', 'ShadowThreshold', value);
        end

        function value = get.ShadowValue(this)
            value = FishVideoHandler_(this.id, 'bsget', 'ShadowValue');
        end
        function set.ShadowValue(this, value)
            FishVideoHandler_(this.id, 'bsset', 'ShadowValue', value);
        end
    end


    methods (Hidden = true)

        function release(this)
            %RELEASE  Closes video file or capturing device.
            %
            % The methods are automatically called by subsequent open() and by destructor.
            %
            % See also cv.VideoCapture.open
            %
            FishVideoHandler_(this.id, 'release');
        end

    end

end