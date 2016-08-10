
COMPUTE = 0;
PLOT = 1;

if ~exist('T','var') || COMPUTE
  [~,~,~,~,T] = xy.Tracker.runTest(50);
end



if PLOT

  figure;
  trange = [20,25];

  
   a = T.plotCenterLine([],trange,0.75,1);
  
  
  
end

  
  