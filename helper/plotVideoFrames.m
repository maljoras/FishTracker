function plotVideoFrames(vidfile)


  if ~nargin 
    vidfile = '/home/malte/Videos/5Zebrafish_nocover_22min_trackingVideo.avi';
  end

  vid = cv.VideoCapture(vidfile);
  
  clf;
  skip = 3;
  nframes = 16;
  startFrame = 940;
  set(vid,'PosFrames',startFrame);
  border = 200;
  shift = 0.02;  
  
  [r1,r2] = getsubplotnumber(nframes);
  
  
  for i = 1:nframes;
    
    if skip
      set(vid,'PosFrames',get(vid,'PosFrames')+skip);
    end
    frame = vid.read();
    frame = frame(:,border:end-border,:);
    

    a = subsubplot(1,1,1,r1,r2,i);
    shiftaxes(a,[-shift/2,-shift/2,shift,shift])
    imagesc(frame);
    axis off
  
  end
  