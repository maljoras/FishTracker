function patch(self,x,y,wx,wy,inColor)
% plots a patch WITHOUT flipping!

  if ~exist('inColor','var')
    inColor = self.defaultColor;
  end
  xx = self.toScreenX(x);
  yy = self.toScreenY(y);
  wxx = self.toScreenWidth(wx);
  wyy = self.toScreenHeight(wy);

  if self.xreversed
    xx = xx-wxx;
  end
  if self.yreversed
    yy = yy-wyy;
  end

  Screen('FillPoly', self.window, inColor*255,...
         [xx,yy;xx+wxx-1,yy;xx+wxx-1,yy+wyy-1;xx,yy+wyy-1;xx,yy]);

end

