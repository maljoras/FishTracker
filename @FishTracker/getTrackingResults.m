function res = getTrackingResults(self,delinvif,forceif)
% RES = FT.GETTRACKINGRESULTS() returns the current results.  RES =
% FT.GETTRACKINGRESULTS(DELINVIF) sets times in RES.POS where a track was lost to
% NaN.
  if isempty(self.res) || (exist('forceif','var') && forceif)
    self.generateResults(); % maybe not done yet
  end
  if ~exist('delinvif','var')
    delinvif = 0;
  end
  if isempty(self.res)
    error('No results available. First track()...');
  else
    res = self.res;
    
    % delete beyond border pixels
    posx = squeeze(self.res.pos(:,1,:));
    posy = squeeze(self.res.pos(:,2,:));

    sz = self.videoHandler.frameSize;
    posx(posx>sz(2) | posx<1) = NaN;
    posy(posy>sz(1) | posy<1) = NaN;

    res.pos(:,1,:) = posx;
    res.pos(:,2,:) = posy;


    if delinvif
      p = self.deleteInvisible('pos');
      res.pos(isnan(p)) = NaN;
    end
    
    % could add an interpolation for NaN here
    % could also add some smoothening of the track pos
  end
end

