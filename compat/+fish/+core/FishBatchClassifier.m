classdef FishBatchClassifier
% dummy for loading the old objects
   
  
  
  methods(Static)
    
    function obj = loadobj(S);
      S.nbody = S.nfish;
      obj = xy.core.BatchClassifier.loadobj(S);
    end
  end
end
