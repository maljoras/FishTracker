function plotByType(self,plottype,plotTimeRange,fishIds)
% PLOTBYTYPE(SELF,PLOTTYPE,PLOTTIMERANGE,FISHIDS) plots the tracking
% results for the given FISHIDS. PLOTTYPE is one of TRACE,
% VELOCITY, PROBMAP, and DOMAINS. These plots are equivalent
% to the coresponding methods PLOTTRACE, PLOTVELOCITY,
% PLOTPROBMAP, ans PLOTDOMAINS. 

  if ~exist('plottype','var') || isempty(plottype)
    plottype = 'TRACE';
  end

  if isempty(self.res)
    warning('No results available. First track()...');
    return;
  end

  if ~exist('fishIds','var') || isempty(fishIds)
    fishIds = 1:self.nfish;
  end

  if ~exist('plotTimeRange','var') || isempty(plotTimeRange)
    plotTimeRange = self.timerange;
  end

  cla;
  res = self.getTrackingResults();
  t = res.tracks.t(:,1);
  plotidx = t>=plotTimeRange(1) & t<plotTimeRange(2);

  centroid = res.tracks.centroid;
  centroidx = centroid(plotidx,fishIds,1);
  centroidy = centroid(plotidx,fishIds,2);

  posx = squeeze(res.pos(plotidx,1,fishIds));
  posy = squeeze(res.pos(plotidx,2,fishIds));
  posx = conv2(posx,ones(5,1)/5,'same');
  posy = conv2(posy,ones(5,1)/5,'same');

  dt = 1/self.videoHandler.frameRate;

  switch upper(plottype)

    case 'TRACE'
      plot(posx,posy);
      title('Traces');
      xlabel('x-Position [px]')
      ylabel('y-Position [px]')
      
      msk = isnan(posx);
      hold on;
      plot(centroidx(msk),centroidy(msk),'.r','linewidth',1);
      
    case 'VELOCITY'
      

      plot(diff(posx)/dt,diff(posy)/dt,'.');
      xlabel('x-Velocity [px/sec]')
      ylabel('y-Velocity [px/sec]')
      
      title('Velocity');
      
    case {'PROBMAP','DOMAINS'}


      dxy = 10;
      szFrame = self.videoHandler.frameSize;
      sz = ceil(szFrame/dxy);
      sz = sz(1:2);
      
      dposx = min(max(floor(posx/dxy)+1,1),sz(2));
      dposy = min(max(floor(posy/dxy)+1,1),sz(1));
      P = zeros([sz,length(fishIds)]);
      for i = 1:size(posx,2)
        msk = ~(isnan(dposy(:,i)) | isnan(dposx(:,i)));
        P(:,:,i) = accumarray([dposy(msk,i),dposx(msk,i)],1,sz);
        P(:,:,i) = P(:,:,i)/sum(sum(P(:,:,i)));
      end
      P(1,1,:) = 0;
      
      if strcmp(upper(plottype),'PROBMAP')
        if size(P,3)==3
          Z = P;
        else
          Z = sum(P,3);
        end
        Z = imfilter(Z,fspecial('disk',3),'same');
        if size(P,3)==3
          image(1:szFrame(2),1:szFrame(1),Z/quantile(P(:),0.99));
        else
          imagesc(1:szFrame(2),1:szFrame(1),Z,[0,quantile(P(:),0.99)]);
        end
        
        title('Overall probability');
        xlabel('x-Position [px]')
        ylabel('y-Position [px]')
        axis xy;
      else

        [~,cl] = max(P,[],3);
        cl = cl + (sum(P,3)~=0);
        imagesc(1:szFrame(2),1:szFrame(1),cl);
        colormap(jet(size(P,3)+1));
        xlabel('x-Position [px]')
        ylabel('y-Position [px]')
        title('Domain of fish')
        axis xy;
      end
      
    otherwise
      error('do not know plot type');
  end
end