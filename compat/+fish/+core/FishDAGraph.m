classdef FishDAGraph
% dummy for loading the old objects
   
  
  
  methods(Static)
    
    function obj = loadobj(S);
      S.nindiv = S.nfish;
      obj = xy.core.DAGraph.loadobj(S);
    end
  end
end
