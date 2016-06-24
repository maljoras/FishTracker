function out = getBendingLaplace(self,fishIds, plotTimeRange,res)

  if nargin<3 || isempty(plotTimeRange)
    plotTimeRange = self.timerange;
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

  centerLine = self.deleteInvisible('centerLine',[],res);
  clx = permute(centerLine(plotidx,fishIds,1,:),[4,1,2,3]);
  cly = permute(centerLine(plotidx,fishIds,2,:),[4,1,2,3]);

  mclx = bsxfun(@minus,clx,nanmean(clx,1));
  mcly = bsxfun(@minus,cly,nanmean(cly,1));
  mclx = convn(mclx,ones(2,1)/2,'valid');
  mcly = convn(mcly,ones(2,1)/2,'valid');

  
  % estimate body orientation from the center line because it is
  % already corrected into the direction of movement (in contrast to
  % segment_orientation field)
  ori = atan2(mclx(1,:,:) - mclx(end,:,:),mcly(1,:,:) - mcly(end,:,:));
  
  
  % direction of movement 
  vel =  permute(self.deleteInvisible('velocity',[],res),[3,1,2]);
  vel = vel(:,plotidx,fishIds);
  vori = atan2(vel(1,:,:),vel(2,:,:));

  %project onto normal of body axis
  ocly = bsxfun(@times,cos(ori),mclx) - bsxfun(@times,sin(ori),mcly);
  %oclx = bsxfun(@times,cos(ori),mclx) + bsxfun(@times,sin(ori),mcly);
  
  % laplace on normal
  olap = nanmean(diff(ocly,2,1));
  
  
  out.t = t(plotidx);
  out.lap = shiftdim(olap,1);
  out.clx = clx;
  out.cly = cly;
  out.ori = shiftdim(ori,1);
  out.vori = shiftdim(vori,1);
  out.vel = vel;
  