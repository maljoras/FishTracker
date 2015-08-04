classdef FishVideoReaderMatlab < FishVideoReader;
% implements the reader using the matlab included video reader
  
  
  properties 
    reader = [];
  end
  
  methods
      
   function dur = a_getDuration(self);
     dur = self.reader.Duration;
   end
   
   function nframes = a_getNFrames(self);
     nframes = floor(self.reader.Duration*self.reader.FrameRate);
   end

   function framesize = a_getFrameSize(self);
     framesize = [self.reader.Height,self.reader.Width];
   end

   function framerate = a_getFrameRate(self);
     framerate = self.reader.FrameRate;
   end
   
   function frame = a_readScaledSFrame(self,scale,delta);      
     frame = self.reader.readFrame();
     frame = im2single(frame);
     frame = bsxfun(@plus,bsxfun(@times,frame,shiftdim(scale(:),-2)),shiftdim(delta(:),-2));
     frame = sum(frame,3);
   end
   
   function frame = a_readScaledUFrame(self,scale,delta);      
     frame = self.reader.readFrame();
     frame = im2single(frame);
     frame = bsxfun(@plus,bsxfun(@times,frame,shiftdim(scale(:),-2)),shiftdim(delta(:),-2));
     frame = uint8(255*sum(frame,3));
   end

   function frame = a_readUFrame(self);      
     % returns uint color frame
     frame = self.reader.readFrame();
   end
   
   function frame = a_readGrayUFrame(self);      
     frame = self.reader.readFrame();
     frame = rgb2gray(frame);
   end
   
   function frame = a_readGraySFrame(self);      
     frame = self.reader.readFrame();
     frame = im2single(rgb2gray(frame));
   end
   
   function frame = a_readSFrame(self);      
     % returns single color frame 
     frame = self.reader.readFrame();
     frame = im2single(frame);
   end
   
   function bool = a_hasFrame(self);
     bool = self.reader.hasFrame()
   end
   
   
   function a_setCurrentTime(self,time);
     self.reader.CurrentTime = time;
   end

   function a_delete(self)
     self.reader = [];
   end
   
    function playMsk(self);
      vp = vision.VideoPlayer('Name',self.videoFile);
      cont = self.hasFrame();
      fgbg = vision.Autothresholder();
      fgbg.InputRange = [0,1];
      while cont
        self.readFrameFormat('GRAYS');
        msk = fgbg.step(self.frame);
        step(vp,msk);
        cont = self.hasFrame() && isOpen(vp);
      end
      release(vp);
      self.reset()
    end
    
 
   function self = FishVideoReaderMatlab(vid,trange,varargin) 
     
     if nargin==1
       trange = [];
       varargin = {};
     end
       
     self@FishVideoReader(vid,varargin{:});
     self.reader = VideoReader(self.videoFile);

     self.init();
     self.timeRange = trange;
     self.reset();
     self.verbose();
     
   end
   
  end
  
  
end
