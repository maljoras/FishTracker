

function calibrateStimulusScreen(self);

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

  self.displayif = 0;
  self.stimulusPresenter = FishStimulusPresenterCalibration(self.opts.stimulus);
  if ~isempty(self.fishlength)
    self.stimulusPresenter.width = self.fishlength;
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

  % get the left track:        
  idx = [];
  [~,idx(1)] = nanmin(mpos(1,:));
  xmin = [bbox(:,idx(1),1)-bbox(:,idx(1),3)/2, bbox(:,idx(1),2)];
  [~,idx(2)] = nanmax(mpos(1,:));
  xmax = [bbox(:,idx(2),1)+bbox(:,idx(2),3)/2, bbox(:,idx(2),2)];

  [~,idx(3)] = nanmin(mpos(2,:));
  ymin = [bbox(:,idx(3),1), bbox(:,idx(3),2)-bbox(:,idx(3),4)/2];
  [~,idx(4)] = nanmax(mpos(2,:));
  ymax = [bbox(:,idx(4),1), bbox(:,idx(4),2)+bbox(:,idx(4),4)/2];

  if length(unique(idx))~=length(idx)
    error(['Something went wrong with the calibration. Screen ' ...
           'not aligned to video edges ?'])
  end




  keyboard
  %restore
  self.stimulusPresenter = stimPresenter;
  self.nfish = nfish;
  self.opts.classifier.nFramesForInit = nFrames;
  self.videoHandler.inverted = inverted;
  self.displayif = displayif;
  self.saveFields = sf;
end


