function  varargout = drawTexture(self,relidx,rect)
  if relidx>length(self.textureIdx) || relidx<1
    error('Cannot find texture. Forgot to init ?');
  end
  if nargin<3
    rect = [0,0,1,1];
  end

  wrect = self.toScreenRect(rect);


  Screen('drawTexture',self.window,self.textureIdx(relidx),[],wrect);

  if nargout
    varargout{1} = self.flip();
  end
end

