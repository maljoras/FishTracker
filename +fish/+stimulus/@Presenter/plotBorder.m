function timestamp = plotBorder(self,width,outColor,inColor)
% [TIMESTAMP] = PLOTBORDER(SELF,WIDTH,OUTCOLOR,INCOLOR) plots a border
% of color OUTCOLOR with width WIDTH. Inside color is INCOLOR.

  if ~exist('outColor','var')
    outColor = self.defaultColor;
  end
  if ~exist('inColor','var')
    inColor = BlackIndex(self.screen);
  end

  width = min(width,0.5);
  self.patch(0,0,1,1,outColor);
  self.patch(width,width,1-2*width,1-2*width,inColor);

  if nargout
    timestamp = self.flip();
  end
end
