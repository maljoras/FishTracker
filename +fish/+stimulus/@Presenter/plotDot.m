function timestamp = plotDot(self,x,y,inSize,inColor)
% plots a dot in normalized coordinates. Width in pixel
  if isempty(x) || isempty(y)
    timestamp = NaN;
    return
  end

  if ~exist('inColor','var')
    inColor = self.defaultColor;
  end
  if ~exist('inSize','var')
    inSize = 20;
  end
  assert(length(x)==length(y))

  xx = self.toScreenX(x);
  yy = self.toScreenY(y);

  for i = 1:length(xx)
    s2 = inSize(min(i,end))/2;
    rect = [xx(i)-s2,yy(i)-s2,xx(i)+s2,yy(i)+s2];
    Screen('FillOval', self.window, inColor(min(i,end),:)*255, rect);
  end

  if nargout
    timestamp = self.flip();
  end

end
