function [d,ddag] = parameterfun(opts1,parname,parvalue,tmax);

    maxNumCompThreads(4);
    opts = opts1;

    opts = fish.helper.setfield(opts,parname,parvalue);

    ft = [];
    [~,t_elapsed,ftposdag,idpos, ft] = fish.Tracker.runTest(tmax,opts,[],[],0);
    ddag = squeeze(sqrt(sum((ftposdag - idpos).^2,2)));
    r = ft.getTrackingResults([],[],0);
    ftpos = r.pos;
    d = squeeze(sqrt(sum((ftpos - idpos).^2,2)));

  end
  
