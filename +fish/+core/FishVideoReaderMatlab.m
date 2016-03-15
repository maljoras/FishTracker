classdef FishVideoReaderMatlab < FishVideoReader;
% implements the reader using the matlab included video reader
  
  
% $$$   properties(Access=protected)
% $$$     reader = [];
% $$$     originalif = false;
% $$$     scale = [1,1,1]/3;
% $$$     delta = 0;  
% $$$   end
  
  methods(Access=protected)
      
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
   
   function [frame oframe] = a_readScaledSFrame(self,scale,delta);      
     oframe = self.reader.readFrame();
     frame = im2single(oframe);
     frame = bsxfun(@plus,bsxfun(@times,frame,shiftdim(scale(:),-2)),shiftdim(delta(:),-2));
     frame = sum(frame,3);
   end
   
   function [frame oframe] = a_readScaledUFrame(self,scale,delta);      
     oframe = self.reader.readFrame();
     frame = im2single(oframe);
     frame = bsxfun(@plus,bsxfun(@times,frame,shiftdim(scale(:),-2)),shiftdim(delta(:),-2));
     frame = uint8(255*sum(frame,3));
   end

   function [frame,oframe] = a_readUFrame(self);      
     % returns uint color frame
     oframe = self.reader.readFrame();
     frame = oframe;
   end
   
   function [frame oframe] = a_readGrayUFrame(self);      
     oframe = self.reader.readFrame();
     frame = rgb2gray(oframe);
   end
   
   function [frame oframe] = a_readGraySFrame(self);      
     oframe = self.reader.readFrame();
     frame = im2single(rgb2gray(oframe));
   end
   
   function [frame oframe]= a_readInvertedGrayUFrame(self);      
     oframe = self.reader.readFrame();
     frame = uint8(255) - rgb2gray(oframe);
   end
   
   function [frame oframe] = a_readInvertedGraySFrame(self);      
     oframe = self.reader.readFrame();
     frame = single(1.0) - im2single(rgb2gray(oframe));
   end

   
   function [frame oframe] = a_readSFrame(self);      
     % returns single color frame 
     oframe = self.reader.readFrame();
     frame = im2single(oframe);
   end
   
   function bool = a_hasFrame(self);
     bool = self.reader.hasFrame();
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

    function a_startReader(self)
      self.reader = VideoReader(self.videoFile);
    end

  end
  
  
  methods
    
 
   function self = FishVideoReaderMatlab(vid,trange,varargin) 
     
     if nargin==1
       trange = [];
       varargin = {};
     end
       
     self@FishVideoReader(vid,varargin{:});


      
   end
   
  end
  
  
end
