


LOAD = 0;
PLOT = 1;
SAVEIF = 1;
COMPUTE = 1;

if LOAD || ~exist('T2','var')
  

end

if COMPUTE 
  t = [];
  tic;T3.track([0,100]);t(3) = toc;  
  tic;T2.track([0,100]);t(2) = toc;  
  tic;T1.track([0,100]);t(1) = toc;
  
  t
end

if PLOT

end
