function addSaveFields(self,varargin)
%  XYT.ADDSAVEFIELD('FIELDNAME1','FIELDNAME2',...) adds a tracks field to the saved structure

  f = fieldnames(self.initializeTracks())'; 
  forbiddenFields = {'id','segment','classProbHistory','predictor',...
                     'crossedTrackIds','batchFeatures','features'};
  for i = 1:length(forbiddenFields)
    f(strcmp(f,forbiddenFields{i})) = [];
  end

  for i = 1:length(varargin)
    field = varargin{i};
    assert(ischar(field));
    
    if ~any(strcmp(field,self.saveFields))
      if ~any(strcmp(field,f)) && isempty(strfind(field,'segment.')) 
        xy.helper.verbose('Available track fields:\n');
        disp(f')
        error(sprintf('Field "%s" no available. Chose one of the above or a "segment." field.',field))
      else
        % add
        self.saveFields{end+1} = field;
      end
    end
    
  end
end
