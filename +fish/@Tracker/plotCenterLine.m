function plotCenterLine(self,fishIds,plotTimeRange)
%   PLOTCENTERLINE(SELF,FISHIDS,PLOTTIMERANGE) plots the traces
%   with center line information. 
  
  if ~exist('fishIds','var') || isempty(fishIds)
    fishIds = 1:self.nfish;
  end

  if ~exist('plotTimeRange','var') || isempty(plotTimeRange)
    plotTimeRange = self.timerange;
  end

  if diff(plotTimeRange) > 100
    error('Time range longer than 100 seconds!');
  end


  res = self.getTrackingResults();
  t = res.tracks.t(:,1);
  plotidx = t>=plotTimeRange(1) & t<plotTimeRange(2);

  if  ~isfield(res.tracks,'centerLine') || isempty(res.tracks.centerLine)
    error('No centerLine data');
  end

  clx = permute(res.tracks.centerLine(plotidx,fishIds,1,:),[4,1,2,3]);
  cly = permute(res.tracks.centerLine(plotidx,fishIds,2,:),[4,1,2,3]);
  clx = convn(clx,ones(2,1)/2,'valid');
  cly = convn(cly,ones(2,1)/2,'valid');
  lap = abs(nanmean(diff(clx,2,1))) + abs(nanmean(diff(cly,2,1),1));
  mx = nanmean(clx(:,:,:),1);
  my = nanmean(cly(:,:,:),1);

% $$$       % find sudden changes
% $$$       clf;
% $$$       %tt = t(plotidx); 
% $$$       tt = find(plotidx);% t is wrong in the txt file...
% $$$       v2 = squeeze(abs(diff(mx,1,2)) + abs(diff(my,1,2))); 
% $$$       a(1) = subplot(3,1,1);
% $$$       plot(tt(1:end-1),v2);
% $$$       title('Difference position')
% $$$       a(2) = subplot(3,1,2);
% $$$       if isfield(res.tracks,'consecutiveInvisibleCount')
% $$$         seq = res.tracks.consecutiveInvisibleCount(plotidx,fishIds);
% $$$         plot(tt,seq);
% $$$         title('Consec. invisible counts');
% $$$       end
% $$$       
% $$$       a(3) = subplot(3,1,3);
% $$$       prob = res.tracks.classProb(plotidx,fishIds,fishIds);
% $$$       plot(tt,prob(:,1:length(fishIds)+1:end));
% $$$       title('Class prob');
% $$$       linkaxes(a,'x');


  figure;
  cmap = jet(self.nfish);
  szFrame = self.videoHandler.frameSize;

  cla;

  for i = 1:length(fishIds)
    col = cmap(fishIds(i),:);
    
    plot(clx(:,:,i),cly(:,:,i),'color',col,'linewidth',2);
    hold on;        
    plot(clx(1,:,i),cly(1,:,i),'o','color',col,'linewidth',1,'markersize',4);

    idx = lap(1,:,i) > 1.5;
    %scatter(mx(1,idx,i),my(1,idx,i),50,'r','o','filled'); 
    plot(clx(:,idx,i),cly(:,idx,i),'color','r','linewidth',1);
  end
  xlabel('X-position [px]')
  ylabel('Y-position [px]');
  xlim([1,self.videoHandler.frameSize(2)])
  ylim([1,self.videoHandler.frameSize(1)])
end
