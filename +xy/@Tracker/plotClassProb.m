function plotClassProb(self,plotTimeRange,identityIds);
% PLOTCLASSPROB(SELF,PLOTTIMERANGE,IDENTITYIDS) plots the difference of
% actualy classprob and the maximal classprob if identity are
% reordered. If there is a consistent difference for some time
% range it might indicate a wrong classficiation.

  if ~exist('identityIds','var') || isempty(identityIds)
    identityIds = 1:self.nindiv;
  end

  if ~exist('plotTimeRange','var') || isempty(plotTimeRange)
    plotTimeRange = self.timerange;
  end

  clf;
  res = self.getTrackingResults(plotTimeRange);
  t = res.t(:,1);

  if ~isfield(res.tracks,'classProb') || isempty(res.tracks.classProb)
    error('No classprob data');
  end

  prob = res.tracks.classProb(:,identityIds,identityIds);
  probid = nanmean(prob(:,1:length(identityIds)+1:end),2);
  if self.nindiv<=5
    p = perms(1:self.nindiv);
    probmax = probid;
    for i = 1:size(p)
      idx = sub2ind([self.nindiv,self.nindiv],(1:self.nindiv)',p(i,:)');
      probmax = nanmax(probmax,nanmean(prob(:,idx),2));
    end
  else
    idx = find(eye(self.nindiv));
    probmax = prob;
    probmax(:,idx) = 0;
    probmax = nanmean(nanmax(probmax(:,:,:),[],3),2);
  end

  nconv = 50;
  probid = conv(probid,ones(nconv,1)/nconv,'same');
  probmax = conv(probmax,ones(nconv,1)/nconv,'same');

  a(1) = subplot(3,1,1);
  plot(t,[probid,probmax]);
  legend('chosen','max other');
  title('Average identity class probabilty')  
  
  a(2) = subplot(3,1,2);
  seq = res.tracks.consecutiveInvisibleCount(:,identityIds);
  plot(t,seq);
  title('Consecutive invisible counts per identity')  ;
  
  a(3) = subplot(3,1,3);
  x = res.pos(:,:,identityIds);
  v2 = squeeze(sqrt(sum(abs(diff(x,1,1)).^2,2))); 
  plot(t(1:end-1),v2);
  linkaxes(a,'x');
  title('Postion distance adjacent frames')  ;
  ylabel('Pixels');
end
