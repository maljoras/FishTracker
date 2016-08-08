function setDefaultResultType(self,dagif)
% SETDEFAULTRESULTTYPE(DAGIF) sets the type of results one get's from
% the gatherring results methods per default if DAGIF = 1 set's to
% DaGraph tracking results, otherwise switchbased.
% SETDEFAULTRESULTTYPE('DAG') or SETDEFAULTRESULTTYPE('SWITCH') can
% also be used.
%
% see also GETTRACKINGRESULTS,GETINVISIBLEMASK)

  if ~ischar(dagif)
    if dagif
      dagif = 'dag';
    else
      dagif = 'switch';
    end
  end
    
  switch lower(dagif)
    case {'dag','dagbased','dag-based','dagraph'}
      xy.helper.verbose('Set default track results to DAGRAPH-based.');
      self.opts.tracks.useDagResults = 1;
    case {'sw','switchbased','switch-based','switch'}
      xy.helper.verbose('Set default track results to SWITCH-based.');
      self.opts.tracks.useDagResults = 0;
    otherwise
      error('Unknown default type (either ''dag'' or ''switch'')');
  end

  
    
  
        