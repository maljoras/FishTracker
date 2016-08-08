


LOAD = 0;
PLOT = 1;
SAVEIF = 1;
COMPUTE = 1;

if LOAD || ~exist('ft2','var')
  

end

if COMPUTE 
  t = [];
  tic;ft3.track([0,100]);t(3) = toc;  
  tic;ft2.track([0,100]);t(2) = toc;  
  tic;ft1.track([0,100]);t(1) = toc;
  
  t
end

if PLOT

end
