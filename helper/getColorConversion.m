function [scale,delta] = getColorConversion(bwmsks,cframes)
% [SCALE,DELTA] = GETCOLORCONVERSION({BWMSKS},{CFRAMES}). Gets a good rgb to
% color conversion
  
  if ~iscell(bwmsks) || ~iscell(cframes)
    error('Expect cell array');
  end
  
  assert(length(bwmsks)==length(cframes));
  
  cfpooled = reshape(permute(cat(4,cframes{:}),[1,2,4,3]),[],3);
  spread = diff(quantile(cfpooled,[0.05,0.95]));
  mu = mean(cfpooled,1);

  col1 = [];
  col2 = [];
  for i = 1:length(cframes)
    cframe = cframes{i};
    cframe1 = reshape(cframe,[],3);
    cframe1 = bsxfun(@rdivide,cframe1,spread);    
    cframe1 = bsxfun(@minus,cframe1,mu);

    bwmsk = bwmsks{i};    
    oidx = find(bwmsk);
    ridx = find(~bwmsk);
    ridx = ridx(floor(rand(length(oidx),1)*length(ridx))+1);
    
    col1 = [col1;double(cframe1(oidx,:))];
    col2 = [col2;double(cframe1(ridx,:))];
  end


  [fv,~,proj,r] = lfd([col1;col2],[ones(size(col1,1),1);zeros(size(col2,1),1)],1,0);
  fv = fv(:)/norm(fv);

  % test directions (objects should be black, background white)
  frame1 = sum(bsxfun(@times,cframes{end},shiftdim(fv./spread(:),-2)),3) + sum(-mu.*fv');
  if mean(frame1(oidx))>mean(frame1(ridx))
    fv = -fv;
  end

  scale = double(fv'./spread);
  delta = double(-mu.*fv' + 0.5/3); % !! needs to be double
  
  if sum(abs(fv))< 1e-5 % makes no sense
    scale = [];
    delta = [];
  end
  
  