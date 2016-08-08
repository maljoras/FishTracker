function plotByType(self,plottype,plotTimeRange,identityIds)
% PLOTBYTYPE(SELF,PLOTTYPE,PLOTTIMERANGE,IDENTITYIDS) plots the tracking
% results for the given IDENTITYIDS. PLOTTYPE is one of TRACE,
% VELOCITY, PROBMAP, and DOMAINS. These plots are equivalent
% to the coresponding methods PLOTTRACE, PLOTVELOCITY,
% PLOTPROBMAP, ans PLOTDOMAINS. 

  if ~exist('plottype','var') || isempty(plottype)
    plottype = 'TRACE';
  end

  if isempty(self.res)
    xy.helper.verbose('WARNING: No results available. First track()...');
    return;
  end

  if ~exist('identityIds','var') || isempty(identityIds)
    identityIds = 1:self.nanimals;
  end

  if ~exist('plotTimeRange','var') || isempty(plotTimeRange)
    plotTimeRange = self.timerange;
  end

  cla;
  res = self.getTrackingResults(plotTimeRange);
  t = res.t(:,1);


  pos = self.getResField(res,'pos',0);
  msk = self.getInvisibleMsk(res);
  invmsk = squeeze(isnan(pos(:,1,identityIds))) | msk(:,identityIds);
  centroid = self.getResField(res,'centroid',0);
  centroidx = centroid(:,identityIds,1);
  centroidy = centroid(:,identityIds,2);

  pos = self.interpolateInvisible(res,'pos',3);
  posx = squeeze(pos(:,1,identityIds));
  posy = squeeze(pos(:,2,identityIds));

  dt = 1/self.videoHandler.frameRate;

  switch upper(plottype)

    case 'TRACE'
      plot(posx,posy);
      title('Traces');
      xlabel('x-Position [px]')
      ylabel('y-Position [px]')
      
      msk = invmsk;
      hold on;
      plot(centroidx(msk),centroidy(msk),'.r','linewidth',1);
      
    case 'VELOCITY'
      

      plot(diff(posx)/dt,diff(posy)/dt,'.');
      xlabel('x-Velocity [px/sec]')
      ylabel('y-Velocity [px/sec]')
      
      title('Velocity');
      
    case {'PROBMAP','DOMAINS'}


      dxy = 15;
      szFrame = self.videoHandler.frameSize;
      sz = ceil(szFrame/dxy);
      sz = sz(1:2);
      
      dposx = min(max(floor(posx/dxy)+1,1),sz(2));
      dposy = min(max(floor(posy/dxy)+1,1),sz(1));
      P = zeros([sz,length(identityIds)]);
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
        p = sum(P,3);
        P(P<max(p(:))*1e-2) = 0;
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
