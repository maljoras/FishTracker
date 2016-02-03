function playVideo(self,timerange,writefile)
% PLAYVIDEO(SELF,TIMERANGE,WRITEFILE) plays (and saves) a video
% together with the tracking results
  
  res = self.getTrackingResults();

  if hasOpenCV() && self.useOpenCV
    vr = @FishVideoReader;
  else
    vr = @FishVideoReaderMatlab;
  end

  videoReader = vr(self.videoFile,self.timerange);

  if isinf(videoReader.duration) 
    verbose('Cannot find videofile: %s',self.videoFile);
    self.videoFile = getVideoFile(self.videoFile);
    videoReader = vr(self.videoFile,self.timerange);
  end

  if ~exist('timerange','var')
    t = res.tracks.t(:,1);
    timerange = t([1,end]);
  end

  if ~exist('writefile','var') || isempty(writefile)
    writefile = '';
  end
  if isscalar(writefile) && writefile
    [a,b,c] = fileparts(self.videoFile);
    writefile = [a '/' b '_playVideo' c];
  end
  videoWriter = [];
  if ~isempty(writefile) 
    videoWriter = self.newVideoWriter();
  end

  if isempty(self.videoPlayer)
    self.videoPlayer = self.newVideoPlayer();
  end
  if  ~isOpen(self.videoPlayer)
    self.videoPlayer.show();
  end

  currentTime = videoReader.currentTime;      
  videoReader.timeRange = timerange;
  videoReader.reset();
  videoReader.setCurrentTime(timerange(1)); 

  t_tracks =res.tracks.t(:,1);
  tidx = find(t_tracks>=timerange(1) & t_tracks<timerange(2)); 
  t_tracks = t_tracks(tidx);
  cols = uint8(255*jet(self.nfish));
  s = 0;
  while videoReader.hasFrame() && s<length(tidx) && isOpen(self.videoPlayer)
    uframe = videoReader.readFrameFormat('RGBU');
    s = s+1;
    
    t = videoReader.currentTime;
    %assert(abs(t_tracks(s)-t)<=1/videoReader.frameRate);

    % Get bounding boxes.
    bboxes = shiftdim(res.tracks.bbox(tidx(s),:,:),1);
    
    % Get ids.

    ids = res.tracks.fishId(tidx(s),:);
    foundidx = ~isnan(ids);
    ids = int32(ids(foundidx));
    labels = cellstr(int2str(ids'));
    clabels = cols(ids,:);

    highcost = isinf(res.tracks.assignmentCost(tidx(s),:));
    highcost = highcost(foundidx);
    clabels(highcost,:) = 127; % grey
                               % Draw the objects on the frame.
    uframe = insertObjectAnnotation(uframe, 'rectangle', bboxes(foundidx,:), labels,'Color',clabels);
    

    clprob = shiftdim(res.tracks.classProb(tidx(s),:,:),1);
    clprob(isnan(clprob)) = 0;
    dia = self.fishlength/self.nfish;
    for i_cl = 1:self.nfish
      for j = 1:length(foundidx)
        if ~foundidx(j)
          continue;
        end
        pos = [bboxes(j,1)+dia*(i_cl-1)+dia/2,bboxes(j,2)+ dia/2 ...
               + bboxes(j,4),dia/2];
        uframe = insertShape(uframe, 'FilledCircle',pos, 'color', cols(i_cl,:), 'opacity',clprob(j,i_cl));
      end
    end

    self.videoPlayer.step(uframe);
    if ~isempty(videoWriter)
      videoWriter.write(uframe);
    end
    
    d = seconds(t);
    d.Format = 'hh:mm:ss';
    verbose('CurrentTime %s\r',char(d))
  end
  fprintf('\n')
  clear videoWriter videoReader
end
