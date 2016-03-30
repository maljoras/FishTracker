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
  b = getBendingLaplace(self,fishIds,plotTimeRange);


  mx = nanmean(b.clx(:,:,:),1);
  my = nanmean(b.cly(:,:,:),1);


  
  clf;
  cmap = jet(self.nfish);
  szFrame = self.videoHandler.frameSize;

  lapthres = 0.75;

  for i = 1:length(fishIds)
    col = cmap(fishIds(i),:);
    
    plot(b.clx(:,:,i),b.cly(:,:,i),'color',col,'linewidth',2);
    hold on;        
    plot(b.clx(1,:,i),b.cly(1,:,i),'o','color',col,'linewidth',1,'markersize',4);

    col1 = 'rk';
    for j = 1:2
      if j==1
        idx = b.lap(:,i) > lapthres;
      else
        idx = b.lap(:,i) < -lapthres;
      end
      
      plot(b.clx(:,idx,i),b.cly(:,idx,i),'color',col1(j),'linewidth',1);
      plot(b.clx(1,idx,i),b.cly(1,idx,i),'o','color',col1(j),'linewidth',1,'markersize',6);

      x =  [b.clx(1,idx,i)', b.clx(1,idx,i)' + 10*cos(b.ori(idx,i))];
      y =  [b.cly(1,idx,i)', b.cly(1,idx,i)' - 10*sin(b.ori(idx,i))];
      plot(x',y','color',col1(j),'linewidth',1);
    end
  end

  
  xlabel('X-position [px]')
  ylabel('Y-position [px]');
  xlim([1,self.videoHandler.frameSize(2)])
  ylim([1,self.videoHandler.frameSize(1)])
end
