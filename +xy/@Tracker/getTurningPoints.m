function out = getTurningPoints(self, res, plotif)
% POINTSTRUC = GETTURNINGPOINTS(SELF,RES) caculated the turning
% points.  GETTRACKINGRESULTS (which uses the default tracking mode
% set by SETDEFAULTRESULTTYPE).  GETTURNINGPOINTS(...,PLOTIF) make a
% dubug plot (default off)
  
  if nargin<3 || isempty(plotif)
    plotif = 0;
  end
  
  if  ~isfield(res.tracks,'centerLine') || isempty(res.tracks.centerLine)
    error('No centerLine data');
  end

  
  centerLine = self.getResField(res,'centerLine',0);
  clx = permute(centerLine(:,:,1,:),[4,1,2,3]);
  cly = permute(centerLine(:,:,2,:),[4,1,2,3]);
  
  cic = res.tracks.consecutiveInvisibleCount;
  mclx = self.interpolateInvisible(shiftdim(mean(clx,1),1),cic);
  mcly =  self.interpolateInvisible(shiftdim(mean(cly,1),1),cic);

  % how many tpoints for fit
  fa = floor(self.avgTimeScale);
  n_pk = max(fa - mod(fa+1,2)); % make ODD
  %n_thr = max(floor(self.avgTimeScale)*2,1);

  % threshols for getting the peaks (in std)
  thres_pks = 0.1;
  %thres_thr = 0.1;

  if plotif
    figure;
  end
  
  
  out =[];
  for i = 1:size(clx,3)
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
    
    
    if plotif
      subplot(2,3,i)
      
      plot(clx(:,:,i),cly(:,:,i),'color',[0.5,0.5,0.5]);
      hold on;
      plot(mx,my,'b')
      
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
    [w,resx,u,resy] = xy.helper.contLinearRegression(y,x,nreg*ncp);

    [mres,loc] = min([resx,resy],[],2);
    nconv = ncp; % smothining over body
    residual = conv(mres,ones(nconv,1)/nconv,'same');
    residual = residual-nanmean(residual);

    nidx = isnan(residual);
    residual(nidx) = 0;
    
    % peaks are the turns
    pks=xy.helper.getPeaks(pm*residual,pkthres,0);
    pks(nidx) = 0;
    pk = find(pks);
    idx = unique(pk - mod(pk,ncp))/ncp + 1; % plus 1 seems to be
                                            % correct !?!

    
    if nargout
      % connect peaks with line and get the ori
      px  = mean(x(:,idx),1); % mean over ncp
      py  = mean(y(:,idx),1);
      ori = atan2(diff(py),diff(px))';
      ori(end+1) = NaN;
    
    end

    
  end
  
  
end




