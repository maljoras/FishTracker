function [screenBoundingBox] = calibrateStimulusScreen(self);

  if ~self.stmif || ~iscell(self.videoFile)
    error('Need to be in stimulation grabbing mode');
  end


  verbose('Starting calibration. Make sure no IR filter is installed! Hit Enter!');
  pause;

  % save
  stimPresenter = self.stimulusPresenter;
  nfish = self.nfish;
  nFrames = self.opts.classifier.nFramesForInit;
  inverted = self.videoHandler.inverted;
  displayif = self.displayif;
  sf = self.saveFields;

  self.displayif = 1;
  self.stimulusPresenter = FishStimulusPresenterCalibration(self.opts.stimulus);
  if ~isempty(self.fishlength)
    self.stimulusPresenter.width = self.fishwidth/2;
  end
  self.stimulusPresenter.applyInit(stimPresenter);

  self.stimulusPresenter.tmax = 30;
  self.stimulusPresenter.freq = 0.2;

  % stmiulature is already init. so take its values

  self.nfish = 4;
  self.opts.classifier.nFramesForInit = Inf;
  self.videoHandler.inverted = 1;
  self.addSaveFields('bbox');

  self.track(); % track the markings
  self.stimulusPresenter.flip(); % turn stim off

  pos = self.deleteInvisible('pos');
  bbox = self.deleteInvisible('bbox');

  % delete possible artifacts at beginning
  offset = ceil(size(pos,1)/10); 
  pos = pos(offset+1:end-offset,:,:);
  bbox = bbox(offset+1:end-offset,:,:);

  mpos = squeeze(nanmean(pos,1));

  % NOTE: For simplicity we assume that the screen is rectangular and
  % that the rectangular video is parallel to the edges of the
  % screen. 
  % We could add a more general affine tranformation if
  % necessary... 

  idx = [];
  [~,idx(1)] = nanmin(mpos(1,:)); % left
  xmin = [bbox(:,idx(1),1), bbox(:,idx(1),2)+bbox(:,idx(1),4)/2];
  [~,idx(2)] = nanmax(mpos(1,:)); % right
  xmax = [bbox(:,idx(2),1)+bbox(:,idx(2),3), bbox(:,idx(2),2)+bbox(:,idx(2),4)/2];

  [~,idx(3)] = nanmin(mpos(2,:)); % bottom
  ymin = [bbox(:,idx(3),1)+bbox(:,idx(3),3)/2, bbox(:,idx(3),2)];
  [~,idx(4)] = nanmax(mpos(2,:)); % top
  ymax = [bbox(:,idx(4),1)+bbox(:,idx(4),3)/2, bbox(:,idx(4),2)+bbox(:,idx(4),4)];

  if length(unique(idx))~=length(idx)
    error(['Something went wrong with the calibration. Screen ' ...
           'not aligned to video edges ?'])
  end

  nbins = 20;
  fsize = self.videoHandler.frameSize;
  xedges = linspace(1,fsize(2),nbins);
  yedges = linspace(1,fsize(1),nbins);

  [~,loc] = histc(ymin(:,1),xedges);
  loc(loc==nbins+1) = 0;
  mymin = accumarray(loc(~~loc),ymin(~~loc,2),[nbins,1],@median,NaN);
  
  [~,loc] = histc(ymax(:,1),xedges);
  loc(loc==nbins+1) = 0;
  mymax = accumarray(loc(~~loc),ymax(~~loc,2),[nbins,1],@median,NaN);

  [~,loc] = histc(xmin(:,2),yedges);
  loc(loc==nbins+1) = 0;
  mxmin = accumarray(loc(~~loc),xmin(~~loc,1),[nbins,1],@median,NaN);

  [~,loc] = histc(xmax(:,2),yedges);
  loc(loc==nbins+1) = 0;
  mxmax = accumarray(loc(~~loc),xmax(~~loc,1),[nbins,1],@median,NaN);

  ydelta = diff(yedges(1:2))/2;
  xdelta = diff(xedges(1:2))/2;
  xyframe = cat(3,[xedges'+xdelta,mymin],[mxmax,yedges'+ydelta],[xedges'+xdelta,mymax],[mxmin,yedges'+ydelta]);
  xyframe = permute(xyframe,[1,3,2]); % x/y 3rd

  screenBoundingBox = [nanmean(mxmin),nanmean(mymin),...
                      nanmean(mxmax)-nanmean(mxmin),...
                      nanmean(mymax)-nanmean(mymin)];
  
  if self.opts.display.calibration
    figure;
    imagesc(self.videoHandler.getCurrentFrame);
    hold on;
    plot(xyframe(:,:,1),xyframe(:,:,2),'linewidth',1);
    rectangle('position',screenBoundingBox,'linewidth',1,'edgecolor','r','facecolor','none')
  
    title('Estimated Screen size');
  end

  verbose('Found Bounding Box [%d,%d,%d,%d]',round(screenBoundingBox));

  %restore and set
  self.stimulusPresenter = stimPresenter;
  self.stimulusPresenter.setScreenSize(screenBoundingBox);
  self.nfish = nfish;
  self.opts.classifier.nFramesForInit = nFrames;
  self.videoHandler.inverted = inverted;
  self.displayif = displayif;
  self.saveFields = sf;
end


