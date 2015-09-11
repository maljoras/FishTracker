function trackAllVideosInPath(path,varargin)

  def.path = '/home/malte/Videos/';
  def.opts.excludeDone = 1;
  def.opts.checkDone = 1;

  def.opts.selectManual = 1;

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
      matname = [p{i} filesep opts.nameadd, b,'.mat'];
      nfish = [];
      if opts.excludeDone && exist(matname,'file');
        
        if opts.checkDone
          % check whether done correctly
          v = load(matname);
          nfish = v.ft.nfish;
          nfish_selected = subShowFrame(fname,nfish);
          if nfish==nfish_selected
            continue;
          else
            nfish = nfish_selected;
          end
          
        else
          continue;
        end
      end
      
      if isempty(nfish) && opts.selectManual
        nfish = subShowFrame(fname,3);
        if isempty(nfish)
          continue; % do not track this video
        end
        
      end
      
      if ~isempty(nfish)
        add = {'nfish',nfish};
      else
        add = {};
      end
      
      
      if opts.parallelif
        parfeval(@subTrack,0,fname,matname,add{:},opts.args{:});
      else
        subTrack(fname,matname,add{:},opts.args{:});
      end
    end
    
    
    
  end
  

function nfish = subShowFrame(fname,nfish)
  if hasOpenCV
    vid = FishVideoReaderCV(fname);
  else
    vid = FishVideoReaderMATLAB(fname);
  end
  
  frame = readFrame(vid);
  
  f = figure;
  p1 = get(f,'position');
  %p1(4) = p1(3)/2;
  set(f,'name','Choose the number of fish:','WindowStyle','modal');

  a = subplot(1,2,1);
  image(frame);
  axis off
  handles = [];
  handles.cancel= 0;
  daspect(a,[1,1,1]);
  p1 = get(a,'position');
  p1(3) = p1(3)*2;
  p1(1) = 0.05;
  set(a,'pos',p1);

  
  title('Choose the number of fish:');
  
  gap = 0.02;
  bw = 0.07;
  p2 = [p1(1)+p1(3)+gap, p1(2), 1-p1(3)-p1(1)-2*gap, p1(4)];
  handles.sld = uicontrol('Style', 'listbox', 'Min',1,'Max',1,'Value',nfish,...
                          'units','normalized','Position', p2,...
                          'String',num2cell([1:50]),'Fontsize',18,...
                  'Callback', @(s,c) eval(['f=gcbf;f.UserData.txt.String = num2str(floor(s.Value));'])); 

  handles.btn = uicontrol('Style', 'pushbutton', 'String', 'OK',...
        'units','normalized','Position',[p2(1),p2(2)-gap-bw,p2(3),bw],...
        'Callback', 'uiresume(gcbf)');       
  handles.btn2 = uicontrol('Style', 'pushbutton', 'String', 'Do Not Track',...
        'units','normalized','Position',[p2(1)-p2(3)-gap,p2(2)-gap-bw,p2(3),bw],...
        'Callback', @(s,c) eval('f=gcbf;f.UserData.cancel = 1;uiresume(gcbf)'));       
  set(f,'UserData',handles);

  uiwait(f); 
  
  if ~ishandle(f)
    error('Window closed');
  end
  
  handles = get(f,'UserData');
  if handles.cancel
    nfish= [];
  else
    nfish = handles.sld.Value;
  end
  close(f);
  
  
function subTrack(fname,matname,varargin);

  ft = FishTracker(fname,'displayif',0,varargin{:});
  ft.track();
  ft.save(matname)
