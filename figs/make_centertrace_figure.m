
COMPUTE = 0;
PLOT = 1;

if ~exist('xyT','var') || COMPUTE
  [~,~,~,~,xyT] = xy.Tracker.runTest(50);
end



if PLOT

  figure;
  trange = [20,25];

  
   a = xyT.plotCenterLine([],trange,0.75,1);
  
  
  
end

  
  