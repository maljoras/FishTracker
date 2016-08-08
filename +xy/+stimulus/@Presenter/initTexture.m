function relidx = initTexture(self,mat)
  self.textureIdx(end+1) = Screen('makeTexture',self.window,mat*255);
  relidx = length(self.textureIdx);
end
