function timestamp = plotVPlane(self,x,inColor1,inColor2)
% plots a vertical half plane at x with lexyT size color1 right color2
  if ~exist('inColor1','var')
    inColor1 = self.defaultColor;
  end
  if ~exist('inColor2','var')
    inColor2 = BlackIndex(self.screen);
  end

  self.patch(0,0,x,1,inColor1);
  self.patch(x,0,1,1,inColor2);

  if nargout
    timestamp = self.flip();
  end
end
