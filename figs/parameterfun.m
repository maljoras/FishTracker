function [d,ddag] = parameterfun(opts1,parname,parvalue,tmax,dfname);

  if nargin>4 && ~isempty(dfname)
    global VERBOSEDIARY;
    VERBOSEDIARY = dfname;
  end
  
  
  maxNumCompThreads(4);
  opts = opts1;

  opts = xy.helper.setfield(opts,parname,parvalue);
  xy.helper.verbose('Start parameter variation: %s = %f',parname,parvalue);

  ft = [];
  [~,t_elapsed,ftposdag,idpos, ft] = xy.Tracker.runTest(tmax,opts,[],[],0);

  
  ddag = squeeze(sqrt(sum((ftposdag - idpos).^2,2)));
  r = ft.getSwitchBasedTrackingResults();
  ftpos = r.pos;
  d = squeeze(sqrt(nansum((ftpos - idpos).^2,2)));
  
  if nargin>4 && ~isempty(dfname)
    VERBOSEDIARY = 0;
  end

  
end

