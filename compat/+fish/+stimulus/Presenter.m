classdef Presenter
% dummy for loading the old objects
   
  methods(Static)
    
    function obj = loadobj(S);
      obj = xy.stimulus.Presenter();
      for f = fieldnames(S)'
        if isprop(obj,f{1})
          obj.(f{1}) = S.(f{1});
        end
      end
    end
  end
end
