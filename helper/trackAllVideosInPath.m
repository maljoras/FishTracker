function trackAllVideosInPath(path,varargin)

  def.path = '/home/malte/Videos/';
  def.opts.excludeDone = 1;
  def.opts.excludeNames = {'toolboxes'};
  def.opts.extension = 'avi';
  def.opts.nameadd = 'a_';
  def.opts.args = {}; % for FishTracker
  def.opts.parallelif = 0;
  
  parseInputs;
  if HELP; return;end
  
  
  p = genpath(path);
  p = strsplit(p,':');
  
  if ~iscell(opts.excludeNames)
    opts.excludeNames = {opts.excludeNames};
  end

  idx = 0;
  for i = 1:length(opts.excludeNames)
    idx = idx | ~cellfun('isempty',strfind(p,opts.excludeNames{i}));
  end
  p(idx) = [];
  
  verbose('Following paths are searched:');
  disp(strjoin(p,'\n'))
  
  
  for i = 1:length(p)

    d = dir([p{i} filesep '*.' opts.extension]);
    
    if isempty(d)
      continue;
    end
    
    for j = 1:length(d)
      fname = [p{i} filesep d(j).name];
      [a,b,c] = fileparts(fname);
      matname = [opts.nameadd, b,'.mat'];
      
      if opts.excludeDone && exist(matname,'file');
        continue;
      end
      keyboard
      if opts.parallelif
        parfeval(@subTrack,0,fname,matname,opts.args{:});
      else
        subTrack(fname,matname,opts.args{:});
      end
    end
    
    
    
  end
  
  
    
  
function subTrack(fname,matname,varargin);

  ft = FishTracker(fname,'displayif',0,varargin{:});
  ft.track();
  ft.save(matname)
