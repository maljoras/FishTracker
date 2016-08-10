function timestamp = plotHPlane(self,y,inColor1,inColor2)
% plots a horizontal half plane at y with left size color1 right color2
  if ~exist('inColor1','var')
    inColor1 = self.defaultColor;
  end
  if ~exist('inColor2','var')
    inColor2 = BlackIndex(self.screen);
  end

  self.patch(0,0,1,y,inColor1);
  self.patch(0,y,1,1,inColor2);

  if nargout
    timestamp = self.flip();
  end
end
