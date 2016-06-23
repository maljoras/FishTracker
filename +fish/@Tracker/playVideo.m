function playVideo(self,timerange,writefile)
% PLAYVIDEO(SELF,TIMERANGE,WRITEFILE) plays (and saves) a video together with the tracking
% results. Rectanble annotation is based in the SWITCH-based tracking results, filled
% color is based in the DAG results. Small dots at the bottom show the current fram
% probability of being a particular fish identity (based in the classifier).
  
  res_dag = self.getTrackingResults(1,0,1);
  res_swb = self.getTrackingResults(1,0,0);

  if ~iscell(self.videoFile)
    videoFile = self.videoFile;
  else
    if ~isempty(self.videoFile{2})
      videoFile = self.videoFile{2};
    else
      error(['Stimulus presentation was done without saving the ' ...
             'videoFile. Cannot play video.']);
    end
  end
  
  
  if fish.helper.hasOpenCV() && self.useOpenCV
    vr = @fish.core.FishVideoReader;
  else
    vr = @fish.core.FishVideoReaderMatlab;
  end

  if ~exist(videoFile) 
    fish.helper.verbose('Cannot find videofile: %s',videoFile);
    fish.helper.verbose('Please select.');
    videoFile = fish.helper.getVideoFile(videoFile);
  end

  videoReader = vr(videoFile,self.timerange);
  
  while isinf(videoReader.duration) 
    fish.helper.verbose('Cannot load videofile: %s',videoFile);
    fish.helper.verbose('Please select.');
    videoFile = fish.helper.getVideoFile(videoFile);
    videoReader = vr(videoFile,self.timerange);
  end

  
  if iscell(self.videoFile) || ~strcmp(videoFile,videoFile)
    yes = input('Overwrite the video file of the class with new one (y/n) [Y]?','s');
    if isempty(yes) || lower(yes(1)=='y')
      self.videoFile = videoFile;
      fish.helper.verbose('Set new videoFile');
    end
  end
  
      
     
  if ~exist('timerange','var')
    t = self.res.tabs(:,1);
    timerange = t([1,end]);
  end

  if ~exist('writefile','var') || isempty(writefile)
    writefile = '';
  end
  if isscalar(writefile) && writefile
    [a,b,c] = fileparts(videoFile);
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

  t_tracks =self.res.tabs(:,1);
  tidx = find(t_tracks>=timerange(1) & t_tracks<timerange(2)); 
  t_tracks = t_tracks(tidx);
  cols = uint8(255*jet(self.nfish));
  s = 0;
  
  plotstm = self.stmif && isfield(self.res,'stmInfo') && isprop(self.stimulusPresenter,'IDX_BBOX');

  while videoReader.hasFrame() && s<length(tidx) && isOpen(self.videoPlayer)
    uframe = videoReader.readFrameFormat('RGBU');
    s = s+1;
    
    t = videoReader.currentTime;
    %assert(abs(t_tracks(s)-t)<=1/videoReader.frameRate);
    ind = find(t_tracks<=t,1,'last');
    if isempty(ind)
      continue;
    end


    %% DAG
    % Get bounding boxes.
    bboxes = shiftdim(res_dag.tracks.bbox(tidx(ind),:,:),1);
    %cboxes = [bboxes(:,1:2) + bboxes(:,3:4)/2, min(bboxes(:,3:4)/2,[],2)];
    % Get ids.
    ids = res_dag.tracks.fishId(tidx(ind),:);
    foundidx = ~isnan(ids);
    ids = int32(ids(foundidx));
    labels = cellstr(int2str(ids'));
    clabels = cols(ids,:);

    uframe = insertShape(uframe, 'Filledrectangle', bboxes(foundidx,:),'Color',clabels,'Opacity',0.5);

    %% SW
    % Get bounding boxes.
    bboxes = shiftdim(res_swb.tracks.bbox(tidx(ind),:,:),1);
    
    % Get ids.
    ids = res_swb.tracks.fishId(tidx(ind),:);
    foundidx = ~isnan(ids);
    ids = int32(ids(foundidx));
    labels = cellstr(int2str(ids'));
    clabels = cols(ids,:);

    uframe = insertObjectAnnotation(uframe, 'rectangle', bboxes(foundidx,:), labels,'Color',clabels);

    %%stimulus
    if plotstm
      stmbboxes = shiftdim(self.res.stmInfo(tidx(ind),:,self.stimulusPresenter.IDX_BBOX),1);
      stmidx = ~isnan(stmbboxes(:,1));
      stmbboxes(stmidx,:) = self.stimulusPresenter.fromScreenBbox( stmbboxes(stmidx,:));
      if isprop(self.stimulusPresenter,'IDX_FISHID')
        stmFishIds = squeeze(self.res.stmInfo(tidx(ind),:,self.stimulusPresenter.IDX_FISHID));
      else
        stmFishIds = self.res.stmFishId(tidx(ind),:);
      end
      clabels = cols(stmFishIds(stmidx),:);
      uframe = insertShape(uframe, 'Filledrectangle', stmbboxes(stmidx,:),'Color',clabels,'Opacity',0.2);
    end
    
    
    %%classprob
    clprob = shiftdim(res_swb.tracks.classProb(tidx(ind),:,:),1);
    clprob(isnan(clprob)) = 0;
    dia = self.fishlength/max(self.nfish,2);
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
    fish.helper.verbose('CurrentTime %s\r',char(d))
  end
  fprintf('\n')
  clear videoWriter videoReader
end
