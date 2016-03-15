function plotClassProb(self,plotTimeRange,fishIds);
% PLOTCLASSPROB(SELF,PLOTTIMERANGE,FISHIDS) plots the difference of
% actualy classprob and the maximal classprob if fish are
% reordered. If there is a consistent difference for some time
% range it might indicate a wrong classficiation.

  if ~exist('fishIds','var') || isempty(fishIds)
    fishIds = 1:self.nfish;
  end

  if ~exist('plotTimeRange','var') || isempty(plotTimeRange)
    plotTimeRange = self.timerange;
  end

  clf;
  res = self.getTrackingResults();
  t = res.tracks.t(:,1);
  plotidx = t>=plotTimeRange(1) & t<plotTimeRange(2);

  if ~isfield(res.tracks,'classProb') || isempty(res.tracks.classProb)
    error('No classprob data');
  end

  prob = res.tracks.classProb(plotidx,fishIds,fishIds);
  probid = nanmean(prob(:,1:length(fishIds)+1:end),2);
  if self.nfish<=5
    p = perms(1:self.nfish);
    probmax = probid;
    for i = 1:size(p)
      idx = sub2ind([self.nfish,self.nfish],(1:self.nfish)',p(i,:)');
      probmax = nanmax(probmax,nanmean(prob(:,idx),2));
    end
  else
    probmax = nanmax(prob(:,:,:),[],3);
  end

  nconv = 50;
  probid = conv(probid,ones(nconv,1)/nconv,'same');
  probmax = conv(probmax,ones(nconv,1)/nconv,'same');

  tt = t(plotidx);


  a(1) = subplot(3,1,1);
  plot([probid,probmax]);

  a(2) = subplot(3,1,2);
  seq = res.tracks.consecutiveInvisibleCount(plotidx,fishIds);
  plot(seq);

  a(3) = subplot(3,1,3);
  x = res.pos(plotidx,:,fishIds);
  v2 = squeeze(sum(abs(diff(x,1,1)),2)); 
  plot(v2);
  linkaxes(a,'x');


end
