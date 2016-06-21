function tracks = step(self,tracks,framesize,t)
% this function will be called from fish.Tracker after each
% round. Calls the stepStimulus method which should be overloaded!
  
  if isempty(tracks)
    return;
  end


  if isempty(self.screenBoundingBox)
    self.setScreenSize([1,1,framesize([2,1])-1]);
  end
  sbbox = self.screenBoundingBox;

  x = nan(length(tracks),1);
  y = nan(length(tracks),1);

  for i =1:length(tracks)
    x(i) = (tracks(i).centroid(1)-sbbox(1))/sbbox(3);
    y(i) = (tracks(i).centroid(2)-sbbox(2))/sbbox(4);
  end

  fishIds = self.getFishIdsFromTracks(tracks);
  stmInfo = self.stepStimulus(x,y,t,fishIds);

  for i =1:length(tracks)
    tracks(i).stmInfo = stmInfo(i,:);
  end

  self.flip();
  
  if self.progressBar
    if ~mod(self.istep,self.progressBarUpdateInt)
      self.updateProgressBar(t,self.tmax);
    end
    self.istep = self.istep + 1;
  end

end
