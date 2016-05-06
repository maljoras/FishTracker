
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
    if any(strcmp(f{1},{'tracks','classifier','blob','detector','reader','stimulus','dag'}))
      for ff = fieldnames(opts.(f{1}))'
        self.opts(1).(f{1}).(ff{1}) = opts.(f{1}).(ff{1});
      end
    end
  end
  self.videoHandler.setOpts(self.opts);
  if ~isempty(self.stimulusPresenter) 
    for f = fieldnames(self.opts.stimulus)'
      if isprop(self.stimulusPresenter,f{1})
        self.stimulusPresenter.(f{1}) = self.opts.stimulus.(f{1});
      end
    end
  end

  if ~isempty(self.daGraph)
    for f = fieldnames(self.opts.dag)'
      if isprop(self.daGraph,f{1})
        self.daGraph.(f{1}) = self.opts.dag.(f{1});
      end
    end
  end
  
end
