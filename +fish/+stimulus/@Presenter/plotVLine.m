function timestamp = plotVLine(self,x,inWidth,inColor)
% plots a dot in normalized coordinates. Width in pixel
  if ~exist('inColor','var')
    inColor = self.defaultColor;
  end
  if ~exist('inWidth','var')
    inWidth = 10;
  end

  xw = inWidth;
  xw = max(min(xw,10),0.5); % max supported

  xx = self.toScreenX(x);

  Screen('DrawLine', self.window, inColor*255,xx,self.toScreenY(0),...
         xx,self.toScreenY(1),inWidth);

  if nargout 
    timestamp =self.flip();
  end

end
