function stmInfo = stepStimulus(self,x,y,t,fishIds);
% plots the stimulus (does NOT need a flip). 

  col = zeros(length(fishIds),3);
  for i = 1:length(fishIds)
    col(i,:) = [1-(fishIds(i)==1),(fishIds(i)==2),1 - (fishId(i)-1)/length(fishIds)];
  end
  self.plotDot(x,y,50,col);     

  stmInfo = [x,y];
end