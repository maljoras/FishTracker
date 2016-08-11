function varargout = plotVelocityMap(self,plotTimeRange,identityIds,minmaxvel)
% PLOTVELOCITYMAP(SELF,PLOTTIMERANGE,IDENTITYIDS,MAXVEL) plots a spatial
% map with averaged velocity [in bodylengths/seconds]. MAXVEL is the
% cutoff velocity.  If MAXVEL<1 it is assumed that MAXVEL is given in
% quantile of velocities. if numel(MAXVEL)==2 MAXVEL(1) is min val and
% MAXVEL(2) is the max vel. In this case a probability map is plotted
% of having a velocity in this range.

  if isempty(self.res)
    xy.helper.verbose('WARNING: No results available. First track()...');
    return;
  end

  if ~exist('identityIds','var') || isempty(identityIds)
    identityIds = 1:self.nindiv;
  end

  if ~exist('plotTimeRange','var') || isempty(plotTimeRange)
    plotTimeRange = self.timerange;
  end


  res = self.getTrackingResults(plotTimeRange);
  t = res.t(:,1);

  pos = self.interpolateInvisible(res,'pos',3);
  velocity = [zeros(1,2,self.nindiv);bsxfun(@rdivide,diff(pos,1),diff(res.tabs))];
  avelocity = squeeze(sqrt(sum(velocity.^2,2)))/self.bodylength;
  avelocity = avelocity(:,identityIds);

  if ~exist('minmaxvel','var')
    maxvel = min(quantile(avelocity(:),0.99),self.maxVelocity);
    minvel = 0;
  else
    if numel(minmaxvel)==1
      maxvel = minmaxvel;
      minvel = 0;
    else
      maxvel = minmaxvel(2);
      minvel = minmaxvel(1);
    end
  end

  if maxvel<=1 % assume quantile
    maxvel = min(quantile(avelocity(:),maxvel),self.maxVelocity);
    minvel = min(quantile(avelocity(:),minvel),self.maxVelocity);
  end
  
  dlmsk = avelocity>maxvel | avelocity<minvel;
  avelocity(dlmsk) = NaN;
  posx = squeeze(pos(:,1,identityIds));
  posy = squeeze(pos(:,2,identityIds));

  
  dxy = 10;
  szFrame = self.videoHandler.frameSize;
  sz = ceil(szFrame/dxy);
  sz = sz(1:2);
  
  dposx = min(max(floor(posx/dxy)+1,1),sz(2));
  dposy = min(max(floor(posy/dxy)+1,1),sz(1));
  P = zeros([sz,length(identityIds)]);
  for i = 1:size(posx,2)
      if nargin<4 || numel(minmaxvel)<2
        msk = ~(isnan(dposy(:,i)) | isnan(dposx(:,i)) | ...
                isnan(avelocity(:,i)));
 
        P(:,:,i) = accumarray([dposy(msk,i),dposx(msk,i)],avelocity(msk,i),sz,@nanmean);
      else
        %probability
        msk = ~(isnan(dposy(:,i)) | isnan(dposx(:,i)));
        P(:,:,i) = accumarray([dposy(msk,i),dposx(msk,i)],~isnan(avelocity(msk,i)),sz,@mean);
      end
  end
  
  P(1,1,:) = 0;
  
  if nargout
    varargout{1} = P;
  else
    
    %plot
    Z = P;
    Z = imfilter(Z,fspecial('disk',3),'same');
    
    [r1,r2] = getsubplotnumber(length(identityIds));

    if length(identityIds)>1
      figure;
      clf;
    else
      cla;
    end

    pmax = quantile(P(:),0.99);
    for i = 1:length(identityIds)
      if length(identityIds)>1
        subplot(r1,r2,i,'align');
      end
      
      imagesc(1:szFrame(2),1:szFrame(1),Z(:,:,i),[0,pmax]);
      
      title(sprintf('#%d',identityIds(i)));
      if ~mod(i-1,r2)
        ylabel('y-Position [px]');
      end
      if i>=r1*r2-r2
        xlabel('x-Position [px]')
      end
      axis xy;
      
    end
    
  end
  