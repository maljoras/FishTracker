


LOAD = 0;
PLOT = 1;
SAVEIF = 1;
COMPUTE = 1;

if LOAD || ~exist('ft2','var')
  
  % NOTE: CANNOT PRELOAD THE OPENCV LIBRARIES FOR THE MATLASB READER TO WORK...
  v =  '/home/malte/Videos/5Zebrafish_nocover_22min.avi';
  ft2 = fish.Tracker(v,'nfish',5,'displayif',0,'detector.adjustThresScale',1,'useMex',0,'useOpenCV',0); 
  ft3 = fish.Tracker(v,'nfish',5,'displayif',0,'detector.adjustThresScale',1,'useMex',0,'useOpenCV',1); 
  
  
  ft1 = fish.Tracker(v,'nfish',5,'displayif',0,'detector.adjustThresScale',1,'useMex',1); 



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
