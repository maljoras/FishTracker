classdef PresenterTrackTextureRegions
% dummy for loading the old objects
   
  methods(Static)
    
    function obj = loadobj(S);
      obj = xy.stimulus.PresenterTrackTextureRegions();
      for f = fieldnames(S)'
        if isprop(obj,f{1})
          try
            obj.(f{1}) = S.(f{1});
          end
        end
      end

    end
  end
end
