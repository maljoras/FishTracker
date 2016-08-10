classdef VideoCapture < handle
    % FISHVIDEOCAPTURE VideoCapture class. Slightly modfied from
    % cv.VideoCapture wrapper taken from the MexOpenCV toolbox
    % https://github.com/kyamagu/mexopencv
    %
    % Class for video capturing from video files or cameras. The class
    % provides Matlab API for capturing video from cameras or for reading
    % video files. Here is how the class can be used:
    %
    %    cap = VideoCapture(source);
    %    pause(3); % Note below
    %    for t = 1:30
    %       imshow(cap.read);
    %       pause(0.1);
    %    end
    %
    % ## Note
    % In some environment, there is a concurrency issue during camera
    % initialization. To avoid unexpected crash, pause for seconds after
    % the initialization of VideoCapture object.
    %

    properties (SetAccess = private)
        id
    end

    methods(Static)
      function bool = installed();
        bool = ~~exist('xyVideoCapture_');
      end
      
    end
    
    
    
    methods
        function this = VideoCapture(filename)
            if nargin < 1, filename = 0; end
            this.id = xyVideoCapture_(filename);
        end
        
        function delete(this)
            xyVideoCapture_(this.id, 'delete');
        end
        
        function varargout = read(this)
            %READ  Grabs, decodes and returns the next video frame
            %
            %    frame = cap.read()
            %
            % The method captures the next video frame and return it. If
            % capturing fails, empty array will be returned instead.
            %
            %
            if nargout==1
              varargout{1}  = xyVideoCapture_(this.id, 'read');
            else
              [varargout{1},varargout{2}]  = xyVideoCapture_(this.id, 'read');
            end
          end
          
        function varargout = readScaledS(this,rgbvec_scale,rgbvec_delta)
            %READ  Grabs, decodes and returns the next video frame
            %and scale each channel according to RGB values in
            %SCALE and shifts channels according to DELTA
            %
            %    frame = cap.read(SCALE,DELTA)
            %
            % The method captures the next video frame and return it. If
            % capturing fails, empty array will be returned instead.
            %
            % See also xyVideoCapture_
            %
            if nargout==1
              varargout{1} = xyVideoCapture_(this.id, 'readScaledS', rgbvec_scale,rgbvec_delta);
            else
              [varargout{1},varargout{2}] = xyVideoCapture_(this.id, 'readScaledS', rgbvec_scale,rgbvec_delta);
            end
            
          end
          
        function varargout = readScaledU(this,rgbvec_scale,rgbvec_delta)
            %READ  Grabs, decodes and returns the next video frame
            %and scale each channel according to RGB values in
            %SCALE and shifts channels according to DELTA. Returns UINT8. 
            %
            %    frame = cap.read(SCALE,DELTA)
            %
            % The method captures the next video frame and return it. If
            % capturing fails, empty array will be returned instead.
            %
            % See also VideoCapture
            %
            if nargout==1
              varargout{1} = xyVideoCapture_(this.id, 'readScaledU', rgbvec_scale,rgbvec_delta);
            else
              [varargout{1},varargout{2}] = xyVideoCapture_(this.id, 'readScaledU', rgbvec_scale,rgbvec_delta);
            end
            
          end
          

        function varargout = readSingle(this)
            %READ  Grabs, decodes and returns the next video
            %frame. Outputs in single format 
            %
            %    frame = cap.readSingle()
            %
            % The method captures the next video frame and return it. If
            % capturing fails, empty array will be returned instead.
            %
            % See also VideoCapture
            %
            if nargout==1
              varargout{1} = xyVideoCapture_(this.id, 'readSingle');
            else
              [varargout{1},varargout{2}] = xyVideoCapture_(this.id, 'readSingle');
            end
          end
          
        function varargout = readGraySingle(this)
            %READ  Grabs, decodes and returns the next video frame
            %in gray and single format
            %
            %    frame = cap.readGraySingle()
            %
            % The method captures the next video frame and return it. If
            % capturing fails, empty array will be returned instead.
            %
            % See also VideoCapture
            %
            if nargout==1
               varargout{1} = xyVideoCapture_(this.id, 'readGraySingle');
            else
              [varargout{1},varargout{2}] = xyVideoCapture_(this.id, 'readGraySingle');
            end
            
        end
        
        function varargout = readInvertedGraySingle(this)
            %READ  Grabs, decodes and returns the next video frame
            %in gray and single format
            %
            %    frame = cap.readGraySingle()
            %
            % The method captures the next video frame and return it. If
            % capturing fails, empty array will be returned instead.
            %
            % See also VideoCapture
            %
            if nargout==1
              varargout{1} = xyVideoCapture_(this.id, 'readInvertedGraySingle');
            else
              [varargout{1},varargout{2}] = xyVideoCapture_(this.id, 'readInvertedGraySingle');
            end
          end
        
        function varargout = readGray(this)
            %READ  Grabs, decodes and returns the next video frame
            %in gray and uint8 format
            %
            %    frame = cap.readGray()
            %
            % The method captures the next video frame and return it. If
            % capturing fails, empty array will be returned instead.
            %
            % See also VideoCapture
            %
            if nargout==1
              varargout{1} = xyVideoCapture_(this.id, 'readGray');
            else
              [varargout{1},varargout{2}] = xyVideoCapture_(this.id, 'readGray');
            end
          end
          
        function varargout = readInvertedGray(this)
            %READ  Grabs, decodes and returns the next video frame
            %in gray and uint8 format
            %
            %    frame = cap.readGray()
            %
            % The method captures the next video frame and return it. If
            % capturing fails, empty array will be returned instead.
            %
            % See also VideoCapture
            %
            if nargout==1
              varargout{1} = xyVideoCapture_(this.id, 'readInvertedGray');
            else
              [varargout{1},varargout{2}] = xyVideoCapture_(this.id, 'readInvertedGray');
            end

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
          
            value = xyVideoCapture_(this.id, 'get', key);
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
            % See also VideoCapture
            %
            xyVideoCapture_(this.id, 'set', key, value);
        end
    end

    methods (Hidden = true)
        function successFlag = open(this, filename)
          if nargin < 1, filename = 0; end
          successFlag = xyVideoCapture_(this.id, 'open', filename);
        end

        function flag = isOpened(this)
          flag = xyVideoCapture_(this.id, 'isOpened');
        end

        function release(this)
          xyVideoCapture_(this.id, 'release');
        end

        function successFlag = grab(this)
          successFlag = xyVideoCapture_(this.id, 'grab');
        end

        function frame = retrieve(this)
          frame = xyVideoCapture_(this.id, 'retrieve');
        end
    end

end
