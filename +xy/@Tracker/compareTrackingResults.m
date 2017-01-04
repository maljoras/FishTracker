function compareTrackingResults(self)
% COMPARETRACKINGRESULTS compares the results with the switch-based
% method an DAG method. Basically plots the overlaps


  eqmsks = self.daGraph.checkOverlap();
  
  r = self.getSwitchBasedTrackingResults();
  rdag = self.getDagTrackingResults();
  n = size(eqmsks,2);
  
  [r1,r2] = xy.helper.getsubplotnumber(n);
  for i = 1:n
    eqmsk = eqmsks(:,i);
    
    mt = size(eqmsk,1);
    nt = self.currentFrame;
    
    d = diff([0;eqmsk;0]);
    stop = find(d==-1)-1 + nt-mt;
    start = find(d==1)-1 + nt-mt ;
    L = stop - start; 
    sL = sort(L,'descend');
    
    if ~isempty(sL)
      mL = sL(1);
      idx = find(L==mL);
      idx = idx(1);
    else
      start = max(self.currentFrame-200,1);
      stop = max(self.currentFrame-100,1);
      idx = 1;
    end
    
    subplot(r1,r2,i)
    subPlotTrack(idx);
  end

  
  
  function subPlotTrack(idx1)
    
    cla;
    offs = 200 ;
    st1 = max(start(idx1) - offs,1);
    st2 = min(stop(idx1) + offs,self.currentFrame);
    dim = 1;
    
    inds = st1:st2;
    %plot(squeeze(r.pos(inds,1,:)),squeeze(r.pos(inds,2,:)),'-');
    hold on;
    x =self.videoHandler.frameSize;
    x = x([2,1]);
    z = [start(idx1),start(idx1),stop(idx1),stop(idx1)];
    patch(z,[0,x(dim),x(dim),0],'r','edgecolor','none','facealpha',0.5);
    
    plot(inds,squeeze(r.pos(inds,dim,:)),'-');
    hold on;
    set(gca,'colororderindex',1);
    %plot(squeeze(rdag.pos(inds,1,:)),squeeze(rdag.pos(inds,2,:)),'--');
    plot(inds,squeeze(rdag.pos(inds,dim,:)),'--');
    xlim(inds([1,end]));
    ylim([0,x(dim)])
    
    
  end
  
end




