function res = getRawTrackData(self)
% RES = GETRAWTRACKDATA(SELF) gets the unsorted savedTracks. They
% are NOT corrected for identity!! ONLY FOR THOSE WHO KNOW WHAT
% THEY ARE DOING.  Otherwise use GETTRACKINGRESULTS instead.
  
  
  if isempty(self.savedTracks.id) 
    xy.helper.verbose('WARNING: cannot generate results!')
    xy.helper.verbose('WARNING: not all fish detected. Maybe adjust "nanimals" setting.');
    return;
  end

  nFrames = size(self.savedTracks.id,3)/self.nanimals;
  if self.currentFrame~=nFrames
    xy.helper.verbose(['WARNING: %d frames got lost (premature abort while ' ...
             'tracking?)'],self.currentFrame-nFrames);
  end
  self.currentFrame = nFrames;
 
  res = [];
  for f = fieldnames(self.savedTracks)'
      if isempty(self.savedTracks.(f{1}))
        continue;
      end
      sz = size(self.savedTracks.(f{1}));
      d = length(sz); % at least 3
      trackdat = permute(self.savedTracks.(f{1}),[d,2,1,3:d-1]);
      fishdat = reshape(trackdat(:,:),[self.nanimals,nFrames,sz(2),sz(1),sz(3:d-1)]);
      res.(f{1}) = permute(fishdat,[2,1,3:d+1]);
    end
end
