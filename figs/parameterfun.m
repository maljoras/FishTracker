function [d,ddag] = parameterfun(opts1,parname,parvalue,tmax,dfname);

  if nargin>4 && ~isempty(dfname)
    global VERBOSEDIARY;
    VERBOSEDIARY = dfname;
  end
  
  
  maxNumCompThreads(4);
  opts = opts1;

  opts = fish.helper.setfield(opts,parname,parvalue);
  fish.helper.verbose('Start parameter variation: %s = %f',parname,parvalue);

  ft = [];
  [~,t_elapsed,ftposdag,idpos, ft] = fish.Tracker.runTest(tmax,opts,[],[],0);

  
  ddag = squeeze(sqrt(sum((ftposdag - idpos).^2,2)));
  r = ft.getTrackingResults(0,[],0);
  ftpos = r.pos;
  d = squeeze(sqrt(sum((ftpos - idpos).^2,2)));
  
  if nargin>4 && ~isempty(dfname)
    VERBOSEDIARY = 0;
  end

  
end

