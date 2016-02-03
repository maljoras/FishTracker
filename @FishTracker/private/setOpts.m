
function setOpts(self,opts)

  if nargin<2
    opts = self.opts;
  end

  for f = fieldnames(opts)'
    if any(strcmp(f{1},{'tracks','classifier','blob','detector','reader','stimulus'}))
      self.opts(1).(f{1}) = opts.(f{1});
    end
  end
  self.videoHandler.setOpts(self.opts);
  if ~isempty(self.stimulusPresenter) && isfield(opts,'stimulus')
    for f = fieldnames(opts.stimulus)'
      if isprop(self.stimulusPresenter,f{1})
        self.stimulusPresenter.(f{1}) = opts.stimulus.(f{1});
      end
    end
  end

end
