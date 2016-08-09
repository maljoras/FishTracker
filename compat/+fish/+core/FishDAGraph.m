classdef FishDAGraph
% dummy for loading the old objects
   
  
  
  methods(Static)
    
    function obj = loadobj(S);
      S.nbody = S.nfish;
      obj = xy.core.DAGraph.loadobj(S);
    end
  end
end
