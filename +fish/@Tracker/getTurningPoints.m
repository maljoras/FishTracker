function out = getTurningPoints(self,fishIds, plotTimeRange,res)

  if nargin<3 || isempty(plotTimeRange)
    plotTimeRange = [-inf,inf];
  end
  
  if nargin<2 || isempty(fishIds)
    fishIds = 1:self.nfish;
  end
  
  if nargin<4 || isempty(res)
    res = self.getTrackingResults();
  end
  
  t = res.t;
  plotidx = t>=plotTimeRange(1) & t<plotTimeRange(2);

  if  ~isfield(res.tracks,'centerLine') || isempty(res.tracks.centerLine)
    error('No centerLine data');
  end

  PLOTIF = 0; % FOR DEBUG
  
  centerLine = self.deleteInvisible('centerLine',[],res);

  clx = permute(centerLine(plotidx,fishIds,1,:),[4,1,2,3]);
  cly = permute(centerLine(plotidx,fishIds,2,:),[4,1,2,3]);
  
  mclx = self.interpolateInvisible(shiftdim(mean(clx,1),1),[],res);
  mcly = self.interpolateInvisible(shiftdim(mean(cly,1),1),[],res);

  % how many tpoints for fit
  n_pk = max(floor(self.avgTimeScale),1);
  %n_thr = max(floor(self.avgTimeScale)*2,1);

  % threshols for getting the peaks (in std)
  thres_pks = 0.1;
  %thres_thr = 0.1;

  if PLOTIF
    figure;
  end
  
  
  out =[];
  for i = 1:length(fishIds)
    y = cly(:,:,i);
    x = clx(:,:,i);

    [idx_pk,ori] = subGetResPeaks(x,y,1,n_pk,thres_pks);
    %idx_thr = subGetResPeaks(x,y,-1,n_thr,thres_thr);

    mx = mclx(:,i);
    my = mcly(:,i);

    out(i).x = mx(idx_pk);
    out(i).y = my(idx_pk);
    out(i).tidx = idx_pk;
    out(i).ori = ori;
    
    
    if PLOTIF
      subplot(2,3,i)
      
      plot(mx,my)
      hold on;
      plot(mx(idx_pk),my(idx_pk),'or','linewidth',2)
      %plot(mx(idx_thr),my(idx_thr),'xb','linewidth',2)
      
      L = 50;
      xx = [mx(idx_pk), mx(idx_pk) + L*cos(ori)]';
      yy = [my(idx_pk), my(idx_pk) + L*sin(ori)]';
      plot(xx,yy,'m','linewidth',2)

      title(i)
    end
    
    
  end
  



  function [idx,ori] = subGetResPeaks(x,y,pm,nreg,pkthres);

    ncp = size(x,1);
    nconv = 2*ncp; % smothining over body
    [w,resx,u,resy] = fish.helper.contLinearRegression(y,x,nreg*ncp);

    [mres,loc] = min([resx,resy],[],2);
    
    residual = conv(mres,ones(nconv,1)/nconv,'same');
    residual = residual-nanmean(residual);

    nidx = isnan(residual);
    residual(nidx) = 0;
    
    % peaks are the turns
    pks=fish.helper.getPeaks(pm*residual,pkthres,0);
    pks(nidx) = 0;
    pk = find(pks);
    idx = unique(pk - mod(pk,ncp))/ncp;

    
    if nargout
      % connect peaks with line and get the ori
      px  = mean(x(:,idx),1); % mean over ncp
      py  = mean(y(:,idx),1);
      ori = atan2(diff(py),diff(px))';
      ori(end+1) = NaN;
    
    end

    
  end
  
  
end




