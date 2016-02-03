
function setOpts(self,varargin)

  if nargin<2
    opts = self.opts;
  elseif nargin==2
    assert(isstruct(varargin{1}))
    opts = varargin{1};
  else
    opts = [];
    assert(~mod(length(varargin),2))
    for i = 1:2:length(varargin)
      assert(ischar(varargin{i}));
      eval(sprintf('opts(1).%s = varargin{i+1};',varargin{i}));
    end
  end

  for f = fieldnames(opts)'
    if any(strcmp(f{1},{'tracks','classifier','blob','detector','reader','stimulus'}))
      for ff = fieldnames(opts.(f{1}))'
        self.opts(1).(f{1}).(ff{1}) = opts.(f{1}).(ff{1});
      end
    end
  end
  self.videoHandler.setOpts(opts);
  if ~isempty(self.stimulusPresenter) && isfield(opts,'stimulus')
    for f = fieldnames(opts.stimulus)'
      if isprop(self.stimulusPresenter,f{1})
        self.stimulusPresenter.(f{1}) = opts.stimulus.(f{1});
      end
    end
  end

end
