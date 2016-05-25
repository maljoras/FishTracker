

COMPUTE = 1;
PLOT = 0;



if COMPUTE
  
  nbins = 50;
  tmax = [];
  
  pars= {'tracks.crossingCostScale','classifier.nFramesAfterCrossing','classifier.crossCostThres','tracks.costOfNonAssignment','tracks.invisibleCostScale','tracks.crossingCostScale','dag.probScale','classifier.reassignProbThres','classifier.handledProbThres','classifier.nlfd','detector.adjustThresScale'};
  pararr = {linspace(0.1,10,nbins),[3:nbins],linspace(0.1,10,nbins),linspace(0.1,10,nbins),linspace(0.1,10,nbins),linspace(0.1,10,nbins),linspace(0.01,2,nbins),linspace(0.01,1,nbins),linspace(0.01,1,nbins),1:nbins,linspace(0.5,1,nbins)};


  

  opts = [];
  opts.verbosity = 1;
  opts.nfish = 5;
  opts.fishlength = 100;
  opts.fishwidth = 20;
  opts.classifier.npca = max(40,nbins);

  clear f
  s = 0; 
  translate = {};
  for i = 1:length(pars)
    for j = 1:length(pararr{i})
      s= s+1;
      f(s) = parfeval(@parameterfun,2,opts,pars{i},pararr{i}(j),tmax);
      translate{s} = [i,j];
    end
  end
  fish.helper.verbose('Submitted %d tasks.',s);
  res = [];

  for k = 1:length(f)
    try
      [completedIdx, d, ddag] = fetchNext(f);
    end
    idx = translate{completedIdx};
    i = idx(1);
    j = idx(2);
    res(i,j).d = d;
    res(i,j).ddag = ddag;
    res(i,j).parval = pararr{i}(j);
    res(i,j).parname = pars{i};
  

    fish.helper.verbose('Fetched %d/%d result.\r',k,length(f));
  end
  save('~/ftparvars.mat','res')
  
end

  
  
if PLOT
  figure;
  [r1,r2] = fish.helper.getsubplotnumber(size(res,1));
  for i = 1:size(res,1)
    parr = cat(2,res(i,:).parval);
    md = shiftdim(nanmean(nanmax(cat(3,res(i,:).d),[],2),1));
    sd = shiftdim(nanstd(nanmax(cat(3,res(i,:).d),[],2),1));
    mdag = shiftdim(nanmean(nanmax(cat(3,res(i,:).ddag),[],2),1));
    sdag = shiftdim(nanstd(nanmax(cat(3,res(i,:).ddag),[],2),1));

    subplot(r1,r2,i)
    fish.helper.errorbarpatch(parr,[md;mdag],[sd;sdag]);
    xlabel(res(i,1).parname);

    ylabel('Mean-max dist');
    
  end
end
