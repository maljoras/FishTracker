function stmInfo = stepStimulus(self,x,y,t,identityIds);
% plots the stimulus (does NOT need a flip). 

  if self.borderWidth
    self.plotBorder(self.borderWidth);
  end
  
  col = zeros(length(identityIds),3);
  for i = 1:length(identityIds)
    col(i,:) = [1-(identityIds(i)==1),(identityIds(i)==2),1 - (identityIds(i)-1)/length(identityIds)];
  end
  self.plotDot(x,y,50,col);     

  stmInfo = [x,y,identityIds(:)];
end