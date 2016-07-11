function hc = smallcolorbar(ax,varargin);
% HC = SMALLCOLORBAR(AX) as COLORBAR but makes a somewhat smaller
% vertical colorbar on the right inside.
% MJR 2/2008

  if ~exist('ax','var')
    ax = gca;
  end
  axes(ax);
  hc = colorbar;
  p2 = get(ax,'position');
  p1 = get(hc,'position');
  p1(1) = p1(1)-0.02;
  p1(3) = p1(3)*0.5;
  set(hc,'position',p1);
  set(ax,'position',p2);

  
  if nargin>1
    delete(hc)
    set(ax,'position',p2);
  end
    
    