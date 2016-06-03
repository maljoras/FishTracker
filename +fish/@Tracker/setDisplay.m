function setDisplay(self,varargin)
% FT.SETDISPLAY(VALUE) sets the amount of plotting during the tracking of the fish. Use 0
% for disabling all plotting. Additional, particular plots can be turned on (see
% help of fish.Tracker).
% Example:
% >> ft.setDisplay(0); % turnoff plotting
% >> ft.setDisplay(3); % turn on and set track plotting level to 3
% >> ft.setDisplay('tracks',true); % turns on tracks display


  if nargin==2
    number = varargin{1};
    assert(isnumeric(number));
    
    self.displayif = ~~number;
    if number >1
      self.opts.display.level = number;
    end
    if number>1
      self.opts.display.tracks = true;
    end
  elseif ~mod(length(varargin),2) && length(varargin)>0
    for i = 1:2:length(varargin)
      assert(ischar(varargin{i}))
      if ~isfield(self.opts.display,varargin{i})
        fprintf('\n');
        disp(fieldnames(self.opts.display));
        error(sprintf('Provide valid display types (see above). "%s" unknown.',varargin{i}));
      end
      assert((isnumeric(varargin{i+1}) ||islogical(varargin{i+1})) ...
             && length(varargin{i+1})==1)
      self.opts.display.(varargin{i}) = varargin{i+1};
    end
    
  else
    fprintf('\n');
    disp(fieldnames(self.opts.display));
    error('Choose one of the above field names to set, E.g. ft.setDisplay(''tracks'',1)');
  end

  if ~self.displayif && nargin>2
    fish.helper.verbose('WARNING: Displaying is turned off. Turn on with ft.setDisplay(1)')
  end

end
