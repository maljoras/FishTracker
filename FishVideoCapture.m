classdef FishVideoCapture < handle
    %FISHVIDEOCAPTURE  FishVideoCapture class. Modfied from
    %cv.VideoCapture wrapper. 
    %
    % Class for video capturing from video files or cameras. The class
    % provides Matlab API for capturing video from cameras or for reading
    % video files. Here is how the class can be used:
    %
    %    cap = FishVideoCapture;
    %    pause(3); % Note below
    %    for t = 1:30
    %       imshow(cap.read);
    %       pause(0.1);
    %    end
    %
    % ## Note
    % In some environment, there is a concurrency issue during camera
    % initialization. To avoid unexpected crash, pause for seconds after
    % the initialization of FishVideoCapture object.
    %
    % See also FishVideoCapture.FishVideoCapture FishVideoCapture.read
    % FishVideoCapture.get FishVideoCapture.set
    %
    
    properties (SetAccess = private)
        % Object ID
        id
    end
    
    methods
        function this = FishVideoCapture(filename)
            %VIDEOCAPTURE  Create a new FishVideoCapture object
            %
            %    cap = FishVideoCapture
            %    cap = FishVideoCapture(devid)
            %    cap = FishVideoCapture(filename)
            %
            % FishVideoCapture create a new instance. With no argument, it
            % connects to the default camera device found in the system.
            % You can specify camera devices by devid, an integer value
            % starting from 0. You can also specify a filename to open a
            % video file.
            %
            % See also FishVideoCapture
            %
            if nargin < 1, filename = 0; end
            this.id = FishVideoCapture_(filename);
        end
        
        function delete(this)
            %DELETE  Destructor of VideoCapture object
            FishVideoCapture_(this.id, 'delete');
        end
        
        function varargout = read(this)
            %READ  Grabs, decodes and returns the next video frame
            %
            %    frame = cap.read()
            %
            % The method captures the next video frame and return it. If
            % capturing fails, empty array will be returned instead.
            %
            % See also FishVideoCapture
            %
            if nargout==1
              varargout{1}  = FishVideoCapture_(this.id, 'read');
            else
              [varargout{1},varargout{2}]  = FishVideoCapture_(this.id, 'read');
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
            % See also FishVideoCapture_
            %
            if nargout==1
              varargout{1} = FishVideoCapture_(this.id, 'readScaledS', rgbvec_scale,rgbvec_delta);
            else
              [varargout{1},varargout{2}] = FishVideoCapture_(this.id, 'readScaledS', rgbvec_scale,rgbvec_delta);
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
            % See also FishVideoCapture
            %
            if nargout==1
              varargout{1} = FishVideoCapture_(this.id, 'readScaledU', rgbvec_scale,rgbvec_delta);
            else
              [varargout{1},varargout{2}] = FishVideoCapture_(this.id, 'readScaledU', rgbvec_scale,rgbvec_delta);
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
            % See also FishVideoCapture
            %
            if nargout==1
              varargout{1} = FishVideoCapture_(this.id, 'readSingle');
            else
              [varargout{1},varargout{2}] = FishVideoCapture_(this.id, 'readSingle');
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
            % See also FishVideoCapture
            %
            if nargout==1
               varargout{1} = FishVideoCapture_(this.id, 'readGraySingle');
            else
              [varargout{1},varargout{2}] = FishVideoCapture_(this.id, 'readGraySingle');
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
            % See also FishVideoCapture
            %
            if nargout==1
              varargout{1} = FishVideoCapture_(this.id, 'readInvertedGraySingle');
            else
              [varargout{1},varargout{2}] = FishVideoCapture_(this.id, 'readInvertedGraySingle');
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
            % See also FishVideoCapture
            %
            if nargout==1
              varargout{1} = FishVideoCapture_(this.id, 'readGray');
            else
              [varargout{1},varargout{2}] = FishVideoCapture_(this.id, 'readGray');
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
            % See also FishVideoCapture
            %
            if nargout==1
              varargout{1} = FishVideoCapture_(this.id, 'readInvertedGray');
            else
              [varargout{1},varargout{2}] = FishVideoCapture_(this.id, 'readInvertedGray');
            end

        end

        function value = get(this, key)
            %GET  Returns the specified FishVideoCapture property
            %
            %   value = cap.get(PropertyName)
            %
            % The method returns a Property value of FishVideoCapture.
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
            % See also FishVideoCapture
            %
            value = FishVideoCapture_(this.id, 'get', key);
        end
        
        function set(this, key, value)
            %SET  Sets a property in the FishVideoCapture.
            %
            %    cap.set(PropertyName, value)
            %
            % The method set a Property value of FishVideoCapture.
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
            % See also FishVideoCapture
            %
            FishVideoCapture_(this.id, 'set', key, value);
        end
    end

    methods (Hidden = true)
        function successFlag = open(this, filename)
            %OPEN  Open video file or a capturing device for video capturing
            %
            % ## Output
            % * __retval__ bool, true if successful
            %
            % The methods release() to close the already opened file or camera.
            % Parameters are the same as in the constructor.
            %
            % See also FishVideoCapture.isOpened
            %
            if nargin < 1, filename = 0; end
            successFlag = FishVideoCapture_(this.id, 'open', filename);
        end

        function flag = isOpened(this)
            %ISOPENED  Returns true if video capturing has been initialized already.
            %
            % ## Output
            % * __retval__ bool, return value
            %
            % If the previous call to constructor or open() succeeded, the method returns true.
            %
            % See also FishVideoCapture.open
            %
            flag = FishVideoCapture_(this.id, 'isOpened');
        end

        function release(this)
            %RELEASE  Closes video file or capturing device.
            %
            % The methods are automatically called by subsequent open() and by destructor.
            %
            % See also FishVideoCapture.open
            %
            FishVideoCapture_(this.id, 'release');
        end

        function successFlag = grab(this)
            %GRAB  Grabs the next frame from video file or capturing device.
            %
            % ## Output
            % * __successFlag__ bool, success flag
            %
            % The function grabs the next frame from video file or camera and returns
            % true (non-zero) in the case of success.
            %
            % The primary use of the function is in multi-camera environments, especially
            % when the cameras do not have hardware synchronization. That is, you call
            % grab() for each camera and after that call the slower method retrieve() to
            % decode and get frame from each camera. This way the overhead on demosaicing
            % or motion jpeg decompression etc. is eliminated and the retrieved frames
            % from different cameras will be closer in time.
            %
            % See also FishVideoCapture.read FishVideoCapture.retrieve
            %
            successFlag = FishVideoCapture_(this.id, 'grab');
        end

        function frame = retrieve(this)
            %RETRIEVE  Decodes and returns the grabbed video frame.
            %
            % The function decodes and returns the just grabbed frame. If no frames has
            % been grabbed (camera has been disconnected, or there are no more frames
            % in video file), the function returns an empty matrix.
            %
            % See also FishVideoCapture.read FishVideoCapture.grab
            %
            frame = FishVideoCapture_(this.id, 'retrieve');
        end
    end

end
