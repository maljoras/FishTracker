function out = deleteInvisible(self,field,timeRange,res)
% X = DELETEINVISIBLE(SELF,field) returns the requested field
% in res.tracks with deleted (NaN) the consecutive invisible
% frames 
  
  if isempty(self.res)
    return;
  end
  if ~exist('timeRange','var') || isempty(timeRange)
    timeRange = self.timerange;
  end
  
  if nargin<4
    res = self.getTrackingResults();
  end
  
  t = res.tracks.t(:,1);
  idx = t>=timeRange(1) & t<timeRange(2);

  msk = self.getInvisibleMsk(res);
  msk = msk(idx,:);
  
  if strcmp(field,'pos')
    out = res.pos(idx,:,:);
    out(permute(cat(3,msk,msk),[1,3,2])) = NaN;
  elseif isfield(res.tracks,field)
    % always nFrames x nFish x nOther
    out = double(res.tracks.(field)(idx,:,:,:,:,:));
    sz = size(out);
    out = reshape(out,numel(msk),[]);
    out(msk,:) = NaN;
    out = reshape(out,sz);
  else
    error(sprintf('Field %s does not exist in res.tracks.',field));
  end
  
end
