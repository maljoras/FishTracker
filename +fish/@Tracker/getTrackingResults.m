function res = getTrackingResults(self,delinvif,forceif,dagif)
% RES = GETTRACKINGRESULTS(DELINVIF,FORCEIF,DAGIF) returns the current
% results.  DELINVIF sets times in RES.POS where a track was lost to
% NaN.  FORCIF==1 forces a re-generation of the pos/res
% structure. DAGIF==0 uses the switch-based tracks, otherwise the
% DAG-based tracks
  
  if isempty(self.res) || (exist('forceif','var') && ~isempty(forceif) && forceif)
    generateResults(self); % maybe not done yet
  end
  if ~exist('delinvif','var') || isempty(delinvif)
    delinvif = 0;
  end
  if ~exist('dagif','var') || isempty(dagif)
    dagif = self.opts.tracks.useDagResults;
  end
  
  if isempty(self.res)
    error('No results available. First track()...');
  else
    if dagif 
      fish.helper.verbose('Getting DAGraph tracking results')
      res = self.res.dag;
    else
      fish.helper.verbose('Getting Switch-based tracking results')
      res = self.res.swb;
    end
  
    for f = fieldnames(self.res)'
      if any(strcmp(f{1},{'swb','dag'}))
        continue;
      end
      res.(f{1}) = self.res.(f{1});
    end
    res = subGetCleanPos(res,delinvif);
  
  end


  function res = subGetCleanPos(res,dif);
  
  
    % delete beyond border pixels
    posx = squeeze(res.pos(:,1,:));
    posy = squeeze(res.pos(:,2,:));

    sz = self.videoHandler.frameSize;
    posx(posx>sz(2) | posx<1) = NaN;
    posy(posy>sz(1) | posy<1) = NaN;

    res.pos(:,1,:) = posx;
    res.pos(:,2,:) = posy;
    

    if dif
      p = self.deleteInvisible('pos',[],res);
      res.pos(isnan(p)) = NaN;
    end
    % could add an interpolation for NaN here
    % could also add some smoothening of the track pos
  
  end
  
  


end

