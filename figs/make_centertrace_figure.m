
COMPUTE = 0;
PLOT = 1;

if ~exist('ft','var') || COMPUTE
  [~,~,~,~,ft] = fish.Tracker.runTest(50);
end



if PLOT

  figure;
  trange = [20,25];

  
   a = ft.plotCenterLine([],trange,0.75,1);
  
  
  
end

  
  