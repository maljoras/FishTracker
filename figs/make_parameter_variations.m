

COMPUTE = 0;
PLOT = 0;
COLLECT = 1;

CONTINUE = 1;
DIARY = 1;


clear global VERBOSELEVEL
clear global VERBOSEDIARY

if COMPUTE
  
  if ~exist('res','var')
    CONTINUE = 0;
  end
  
  
  if ~CONTINUE
    nbins = 100;
    tmax = [];

    pv =[];
    
    % define parameters
    pv.classifier.nFramesAfterCrossing = [3:nbins];
    pv.classifier.crossCostThres = linspace(0.1,10,nbins);
    pv.tracks.costOfNonAssignment = linspace(0.1,10,nbins);
    pv.tracks.invisibleCostScale = linspace(0.1,10,nbins);
    pv.tracks.crossingCostScale = linspace(0.001,2,nbins);
    pv.dag.probScale = linspace(0.01,1,nbins); 
    pv.classifier.reassignProbThres = linspace(0.01,1,nbins);
    pv.classifier.handledProbThres = linspace(0.01,1,nbins);
    pv.classifier.nlfd = 1:39;
    pv.detector.adjustThresScale = linspace(0.5,1,nbins);
    pv.tracks.costtau = round(linspace(1,800,nbins));
    pv.tracks.crossBoxLengthScale = linspace(0.1,3,nbins);
    pv.classifier.nFramesForUniqueUpdate = round(linspace(10,300,nbins));
    pv.classifier.tau = round(linspace(100,10000,nbins));
    pv.detector.history = round(linspace(10,5000,nbins));
    pv.tracks.probThresForFish = linspace(0.01,1,nbins);
    pv.classifier.clpMovAvgTau = linspace(1,50,nbins);

    
    % re-format
    pars= fish.helper.allfieldnames(pv)';
    pararr = {};
    for i = 1:length(pars)
      pararr{i} = fish.helper.getfield(pv,pars{i});
    end
  

    opts = [];
    opts.verbosity = 3;
    opts.nfish = 5;
    opts.fishlength = 100;
    opts.fishwidth = 20;
    opts.classifier.npca = 40;
    
  end
  
    
  clear f
  s = 0; 

  translate = {};

    
  if DIARY
    basepath = '/home/malte/diary/';
  end

  
  for i = 1:length(pars)
    for j = 1:length(pararr{i})
      
      if CONTINUE
        if i<=size(res,1) &&  j<=size(res,2) && ~isempty(res(i,j).d)
          continue;
        end
      end
      
      s= s+1;
      fn = sprintf('%sdiary.%d.log',basepath,s);
      f(s) = parfeval(@parameterfun,2,opts,pars{i},pararr{i}(j),tmax,fn);
      translate{s} = [i,j];
    end
  end
  fish.helper.verbose('Submitted %d tasks.',s);

  if ~CONTINUE
    res = [];
  end
  
end

if COLLECT
  
  for k = 1:length(f)

    ONLINE = 1;
    if ONLINE
      [completedIdx, d, ddag] = fetchNext(f);
    else
      completedIdx = k;
      try
        [d, ddag] = fetchOutputs(f(k));
      catch
        continue;
      end
    end
  
    idx = translate{completedIdx};
    i = idx(1);
    j = idx(2);
    res(i,j).d = d;
    res(i,j).ddag = ddag;
    res(i,j).parval = pararr{i}(j);
    res(i,j).parname = pars{i};
    res(i,j).log = sprintf('%sdiary.%d.log',basepath,completedIdx);
    %read log with: textread(res(i,j).log,'%s','delimiter','\n');
    fish.helper.verbose('Fetched %d/%d results.\r',k,length(f));
  end
  save('~/ftparvars1.mat','res')
end

  
if PLOT
  figure;
  [r1,r2] = fish.helper.getsubplotnumber(size(res,1));
  %fun = @(x)nanmax(x,[],2);
  fun = @(x)nanmean(x,2);

  for i = 1:size(res,1)
    parr = cat(2,res(i,:).parval);
    md = shiftdim(nanmean(fun(cat(3,res(i,:).d)),1));
    sd = shiftdim(stderr(fun(cat(3,res(i,:).d)),1));
    mdag = shiftdim(nanmean(fun(cat(3,res(i,:).ddag)),1));
    sdag = shiftdim(stderr(fun(cat(3,res(i,:).ddag)),1));

    
    
    subplot(r1,r2,i);
    nc = 5;
    fish.helper.errorbarpatch(parr,imfilter([md,mdag],ones(nc,1)/nc,'same','replicate'),[sd,sdag]);
    xlabel(unique(cat(1,res(i,:).parname),'rows'));

    ylabel('Mean-max dist');
    
  end
end
