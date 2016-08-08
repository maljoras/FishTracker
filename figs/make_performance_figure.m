


LOAD = 0;
PLOT = 1;
SAVEIF = 1;
COMPUTE = 1;

if LOAD || ~exist('xyT2','var')
  

end

if COMPUTE 
  t = [];
  tic;xyT3.track([0,100]);t(3) = toc;  
  tic;xyT2.track([0,100]);t(2) = toc;  
  tic;xyT1.track([0,100]);t(1) = toc;
  
  t
end

if PLOT

end
