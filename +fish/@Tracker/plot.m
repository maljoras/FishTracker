function plot(self,varargin)
% FISHTRACKER/PLOT([PLOTTIMERANGE],[FISHIDS]) plots an overview of
% the tracking results

  clf;
  subplot(2,2,1)
  self.plotTrace(varargin{:});

  subplot(2,2,2)
  self.plotVelocity(varargin{:});

  subplot(2,2,3)
  self.plotProbMap(varargin{:});

  subplot(2,2,4)
  self.plotDomains(varargin{:});

end
