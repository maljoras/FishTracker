function [scale,delta] = getColorConversion(bwmsk,cframe)
 % [SCALE,DELTA] = GETCOLORCONVERSION(BWMSK,CFRAME). Gets a good rgb to
 % color conversion
 
  
  oidx = find(bwmsk);
  ridx = find(~bwmsk);
  ridx = ridx(floor(rand(length(oidx),1)*length(ridx))+1);
  

      
  cframe1 = reshape(cframe,[],3);
  spread = diff(quantile(cframe1,[0.05,0.95]));
  cframe1 = bsxfun(@rdivide,cframe1,spread);
  mu = mean(cframe1,1);
  cframe1 = bsxfun(@minus,cframe1,mu);


  
  col1 = double(cframe1(oidx,:));
  col2 = double(cframe1(ridx,:));


  [fv,~,proj,r] = lfd([col1;col2],[ones(size(col1,1),1);zeros(size(col2,1),1)],1,0);
  fv = fv(:)/norm(fv);
  
  frame1 = sum(bsxfun(@times,cframe,shiftdim(fv./spread(:),-2)),3) + sum(-mu.*fv');

  if mean(frame1(oidx))>mean(frame1(ridx))
    fv = -fv;
  end

  scale = double(fv'./spread);
  delta = double(-mu.*fv' + 0.5/3); % !! needs to be double
  
  if sum(abs(fv))< 1e-5 % makes no sense
    scale = [];
    delta = [];
  end
  
    