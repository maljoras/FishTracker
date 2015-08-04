function [scale,delta] = getColorConversion(bwmsk,cframe)
 % [SCALE,DELTA] = GETCOLORCONVERSION(BWMSK,CFRAME). Gets a good rgb to
 % color conversion
 
  
  oidx = find(bwmsk);
  ridx = find(~bwmsk);
  ridx = ridx(floor(rand(length(oidx),1)*length(ridx))+1);
  

      
  cframe1 = reshape(cframe,[],3);
  col1 = double(cframe1(oidx,:));
  col2 = double(cframe1(ridx,:));


  [fv,~,proj,r] = lfd([col1;col2],[ones(size(col1,1),1);zeros(size(col2,1),1)],1,0);
  
  frame1 = sum(bsxfun(@times,bsxfun(@minus,cframe,shiftdim(r.mu,-1)),shiftdim(fv,-2)),3);

  if mean(frame1(oidx))>mean(frame1(ridx))
    fv = -fv;
  end

  scale = fv(:)'/norm(fv);
  delta = -r.mu.*scale;

  if sum(abs(scale))< 1e-5 % makes no sense
    scale = [];
    delta = [];
  end
  
    