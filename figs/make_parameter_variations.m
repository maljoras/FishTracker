

COMPUTE = 1;
PLOT = 1;

if COMPUTE
  
  nbins = 50;
  tmax = [];
  
  pars= {'tracks.crossingCostScale','classifier.nFramesAfterCrossing','classifier.crossCostThres','tracks.costOfNonAssignment','tracks.invisibleCostScale','tracks.crossingCostScale','dag.probScale','classifier.reassignProbThres'};
  pararr = {linspace(0.1,10,nbins),[3:nbins],linspace(0.1,10,nbins),linspace(0.1,10,nbins),linspace(0.1,10,nbins),linspace(0.1,10,nbins),linspace(0.01,2,nbins),linspace(0.01,1,nbins)};


  
  ftpos = [];
  ftposdag = [];
  meandist = {};
  meandistdag = {};
  t_elapsed = {};
  parfor i = 1:length(pars)

    lastn = maxNumCompThreads(4);
    for j = 1:length(pararr{i})
      opts = [];
      opts.verbosity = 1;
      opts.nfish = 5;
      opts.fishlength = 100;
      opts.fishwidth = 20;
      opts = fish.helper.setfield(opts,pars{i},pararr{i}(j));

      [~,t_elapsed{i}(j),ftposdag,idpos, ft] = fish.Tracker.runTest(tmax,opts,[],[],0);

      ddag = nanmax(sqrt(sum((ftposdag - idpos).^2,2)),[],3);
      meandistdag{i}(j) = nanmean(ddag);

      r = ft.getTrackingResults([],[],0);
      ftpos = r.pos;
      d = nanmax(sqrt(sum((ftpos - idpos).^2,2)),[],3);
      meandist{i}(j) = nanmean(d);
    end
    
  end
end


if PLOT
  figure;
  [r1,r2] = fish.helper.getsubplotnumber(length(pars));
  for i = 1:length(pars)
    subplot(r1,r2,i)
    plot(pararr{i},meandist{i});
    hold on;
    plot(pararr{i},meandistdag{i});
    xlabel(pars{i});
    ylabel('Mean-max dist');
  end
end
