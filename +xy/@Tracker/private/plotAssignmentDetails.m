function plotAssignmentDetails(self)          


  figure(5);
  clf;
  [r1,r2] = xy.helper.getsubplotnumber(length(self.segments));
  for i = 1:length(self.segments)
    subplot(r1,r2,i);
    imagesc(self.segments(i).FilledImage);
    if any(i==self.unassignedDetections)
      title(sprintf('unassigned [%1.3f]',outcost(i)));
    else
      title(sprintf('trackID %d [%1.3f]',self.assignments(self.assignments(:,2)==i,1),outcost(i)));
    end
    
  end
  drawnow;

  if self.displayif && self.opts.display.assignmentCost
    figure(6);
    subplot(self.nbody+2,1,1);
    hold on;   
    plot(self.currentFrame,self.meanAssignmentCost','xk');
    title('Assignment costs');
    
    subplot(self.nbody+2,1,2);
    hold on;
    set(gca,'colororderindex',1);
    plot(self.currentFrame,self.meanCost','x');
    title('mean cost');
    
    
    plotcost = nan(self.nbody,size(fcost,3));
    ssfcost = bsxfun(@times,fcost,shiftdim(scale(:),-2));
    for i = 1:size(self.assignments,1)
      plotcost(self.assignments(i,1),:) = ssfcost(self.assignments(i,1),self.assignments(i,2),:);
    end
    

    
    col = 'rgbykmc';
    for i = 1:size(plotcost,1)
      subplot(self.nbody+2,1,2+i);
      hold on;
      for j = 1:size(plotcost,2)
        h(j) = plot(self.currentFrame,plotcost(i,j),['x',col(j)]);
      end
      idx = find(i==self.assignments(:,1));
      if ~isempty(idx)
        plot(self.currentFrame,sfcost(self.assignments(idx,1),self.assignments(idx,2)),'ok');
      end
      ylim([0,self.meanAssignmentCost*5])
    end
% $$$           if self.currentFrame==5
% $$$             legend(h,self.costinfo);
% $$$           end
    

  end
  drawnow;
end
