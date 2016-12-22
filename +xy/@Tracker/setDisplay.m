function setDisplay(self,varargin)
% XYT.SETDISPLAY(VALUE) sets the amount of plotting during the tracking of the bodies. Use 0
% for disabling all plotting. Additional, particular plots can be turned on (see
% help of xy.Tracker).
% Example:
% >> T.setDisplay(0); % turnoff plotting
% >> T.setDisplay(3); % turn on and set track plotting level to 3
% >> T.setDisplay('tracks',true); % turns on tracks display


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
    if number>9
      self.opts.display.displayEveryNFrame = 1;
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
    error('Choose one of the above field names to set, E.g. T.setDisplay(''tracks'',1)');
  end
  
  if ~isempty(self.videoHandler)
    self.videoHandler.plotting(~~self.opts.display.videoHandler);
  end
  
  if ~self.displayif && nargin>2
    xy.helper.verbose('WARNING: Displaying is turned off. Turn on with T.setDisplay(1)')
  end

end
