function [d,ddag] = parameterfun(opts1,parname,parvalue,tmax,dfname);

  if nargin>4 && ~isempty(dfname)
    global VERBOSEDIARY;
    VERBOSEDIARY = dfname;
  end
  
  
  maxNumCompThreads(4);
  opts = opts1;

  opts = xy.helper.setfield(opts,parname,parvalue);
  xy.helper.verbose('Start parameter variation: %s = %f',parname,parvalue);

  xyT = [];
  [~,t_elapsed,xyTposdag,idpos, xyT] = xy.Tracker.runTest(tmax,opts,[],[],0);

  
  ddag = squeeze(sqrt(sum((xyTposdag - idpos).^2,2)));
  r = xyT.getSwitchBasedTrackingResults();
  xyTpos = r.pos;
  d = squeeze(sqrt(nansum((xyTpos - idpos).^2,2)));
  
  if nargin>4 && ~isempty(dfname)
    VERBOSEDIARY = 0;
  end

  
end

