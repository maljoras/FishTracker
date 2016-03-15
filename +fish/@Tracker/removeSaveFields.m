function removeSaveFields(self,varargin);
%  FT.REMOVESAVEFIELD('FIELDNAME1','FIELDNAME2',...) removed a tracks field from the saved structure
%  FT.REMOVESAVEFIELD() removes all but the necessary ones.
  minimalSaveFields = {'centroid','classProb','bbox','assignmentCost', 'velocity', 'consecutiveInvisibleCount', ...
                      'segment.Orientation','segment.Size', 'centerLine', 'thickness'};

  if nargin==1
    fish.helper.verbose('Reset saveFields to minimal fields.');
    self.saveFields = minimalSaveFields;
  else
    for i = 1:length(varargin)
      field = varargin{i};
      assert(ischar(field));
      
      if any(strcmp(field,minimalSaveFields))
        warning(sprintf('Cannot remove field %s since it is required.',field));
      else
        
        idx = find(strcmp(field,self.saveFields));
        if isempty(idx)
          warning(sprintf('Field %s not found in saveFields',field));
        else
          self.saveFields(idx)=[];
        end
      end
    end
  end
end
