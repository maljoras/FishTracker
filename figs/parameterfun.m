function [d,ddag] = parameterfun(opts1,parname,parvalue,tmax,dfname);

  if nargin>4 && ~isempty(dfname)
    global VERBOSEDIARY;
    VERBOSEDIARY = dfname;
  end
  
  
  maxNumCompThreads(4);
  opts = opts1;

  opts = xy.helper.setfield(opts,parname,parvalue);
  xy.helper.verbose('Start parameter variation: %s = %f',parname,parvalue);

  T = [];
  [~,t_elapsed,xyposdag,idpos, T] = xy.Tracker.runTest(tmax,opts,[],[],0);

  
  ddag = squeeze(sqrt(sum((xyposdag - idpos).^2,2)));
  r = T.getSwitchBasedTrackingResults();
  xypos = r.pos;
  d = squeeze(sqrt(nansum((xypos - idpos).^2,2)));
  
  if nargin>4 && ~isempty(dfname)
    VERBOSEDIARY = 0;
  end

  
end

