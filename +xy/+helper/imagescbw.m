function varargout = imagescbw(x,y,X,scale)
% IMAGESCBW(X,SCALE) as IMAGESC but plots in black&white without using
% the colormap
  
  if nargin<=2
    X = x;

    if exist('y','var') 
      scale = y;
    end
  end
  
  if ~exist('scale','var') 
    scale = xy.helper.minmax(double(X(:)'));
  end
  
  
  
  scX = max(min((X-scale(1))/diff(scale),1),0);

  if nargin>2
    h = image(x,y,repmat(scX,[1,1,3]));
  else
    h = image(repmat(scX,[1,1,3]));    
  end
  
  