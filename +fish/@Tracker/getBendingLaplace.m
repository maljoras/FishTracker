function out = getBendingLaplace(self,res)

  
  if  ~isfield(res.tracks,'centerLine') || isempty(res.tracks.centerLine)
    error('No centerLine data');
  end

  centerLine = self.deleteInvisible(res,'centerLine');
  clx = permute(centerLine(:,:,1,:),[4,1,2,3]);
  cly = permute(centerLine(:,:,2,:),[4,1,2,3]);

  mclx = bsxfun(@minus,clx,nanmean(clx,1));
  mcly = bsxfun(@minus,cly,nanmean(cly,1));
  mclx = convn(mclx,ones(2,1)/2,'valid');
  mcly = convn(mcly,ones(2,1)/2,'valid');

  
  % estimate body orientation from the center line because it is
  % already corrected into the direction of movement (in contrast to
  % segment_orientation field)
  ori = atan2(mclx(1,:,:) - mclx(end,:,:),mcly(1,:,:) - mcly(end,:,:));
  
  
  % direction of movement 
  vel =  permute(self.deleteInvisible(res,'velocity'),[3,1,2]);
  vori = atan2(vel(1,:,:),vel(2,:,:));

  %project onto normal of body axis
  ocly = bsxfun(@times,cos(ori),mclx) - bsxfun(@times,sin(ori),mcly);
  %oclx = bsxfun(@times,cos(ori),mclx) + bsxfun(@times,sin(ori),mcly);
  
  % laplace on normal
  olap = nanmean(diff(ocly,2,1));
  
  
  out.t = res.t;
  out.lap = shiftdim(olap,1);
  out.clx = clx;
  out.cly = cly;
  out.ori = shiftdim(ori,1);
  out.vori = shiftdim(vori,1);
  out.vel = vel;
  