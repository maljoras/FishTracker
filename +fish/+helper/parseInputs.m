% PARSEINPUTS sets the default values in a M-file.
%
% USAGE in M-file:
%
%   function myfun(arg1,arg2,varargin) %%%note the varargin for the options
%   %help text
%   
%   %set default values. %%%Note that ordering DOES matter: set from first
%                        %%%argument to last 
%   def.arg1 = 2;                           
%   doc.arg1 = 'documentation of first argument';  
%
%   def.arg2 = [];       
%   doc.arg2 = {'the second argument' 'has a very long documentation'...
%               'which can be split in multiple lines'}; 
%
%   def.opts.option1 = 1;
%   doc.option1 = 'the first option'; %%% NOTE the "doc.option1" instead of "doc.opts.option1"
%   
%   def.opts.option2 = 2;
%   doc.option2 = 'the second option'; 
%
%   %%%... arbitrary more options
%
%   parseInputs;           %%% CAUTION: define no other variables before that
%                          %%% call, than 'def' 'doc' and 'MFILENAME' 
%                          %%% (and additionally you can set VERBOSELEVEL
%                          %%% to reduce the verbosity, which can also be
%                          %%% given as option) 
%
%   if HELP;return;end     %%% to avoid calculation when no argument is
%                          %%% given. Instead helptext and default values
%                          %%% together with the documentation will be displayed  
%
%   % function body
%   %%% - now we dont need to check if an argument exist: it always exists
%   %%%   (with default value if not set explicitly in function call)
%   %%% - all options are available in the 'opts' struc.
%   %%% - additionally PARSEINPUTS provides an INPUTS structure which includes
%   %%%   all given inputs.
%
%   [....]
%
%
%
% CALL of M-file
% Now there are multiple possibilities to call MYFUN from above.  
%
%  1)  MYFUN  (without parameters) will display help, default parameters
%      and associated documentation
%
%  2)  MYFUN(value1)  (missing arguments) will use value1 for first argument, but
%      defaults for options and second argument  
%
%  3)  MYFUN(value1,value2,'option2',ovalue2) (options
%      as par-value pairs) sets option2 to ovalue2 and arguments to
%      value1 and value2,  but uses defaults for 'option1'
%
%  4)  MYFUN(value1,[],'option2',ovalue2,'option1',ovalue1) (omitting
%      args) will set arg2 to default value. 
%
%  5)  MYFUN(value1,value2,opts) (options given as structure) sets the
%      options according to a structure opts, which has (a subset of)
%      option names as field names,e.g:
%      >> opts.option1 = 1;opts.option2 = 2;
%      
%  6)  MYFUN(value1,value2,opts,'option2',4)  (options as struc and
%      par-value pairs) will set option2 to 4 instead of 2 (with the
%      opts structure from 5) ).

%
% MJR  14.2.2006,  V1.01,  NOW: empty opts values get default values
%      21.3.2006           supports doc strings
%      15.5.2006           gets help for object methods with MYFUN(obj,'HELP')
%      25.12.2007          gets help for object methods with HELP=1
%                          defined in the preamble. One should use
%                          MFILENAME = [class(self) mfilename] to get it working
%      15.12.2010          omit MFILENAME. now with dbstack. however have
%                          not tested for classes yet
%      26.1.2011           reduced verbosity when called within functions
%      24.8.2012           Fixed the var order in verbose when called

if ~exist('STATIC','var')
  if exist('WW','var') | exist('INPUTS','var') | exist('j','var') | ...
      exist('i','var') | exist('SDEF','var') | exist('WWW','var') | ...
      exist('INFOSTR','var') |       exist('EMPTYSTR','var')
    error('please use different input var names!')
  end
end


if ~exist('MFILENAME','var')
  MFILENAME = dbstack;
  MFILENAME = MFILENAME(2).name;
end

  
%set empty values to default
WW = whos;
WW = {WW.name};
if ~exist('HELP','var') % this can be used to get help directly
  HELP = 0;
end
if exist('def','var');
  INPUTS = def; % get the order right
end
if ~exist('ALLFIELDNAMES','var')
  ALLFIELDNAMES = false;
end

for i = 1:length(WW)
  if strcmp(WW{i},'def'), continue; end;
  if strcmp(WW{i},'doc'), continue; end;
  if strcmp(WW{i},'varargin'), continue; end;
  if strcmp(WW{i},'MFILENAME'), continue; end;
  if strcmp(WW{i},'VERBOSELEVEL'), continue; end;
  if strcmp(WW{i},'ALLFIELDNAMES'), continue; end;
  if strcmp(WW{i},'STRICT'), continue; end; 
  if strcmp(WW{i},'HELP'), continue; end; 
  
  INPUTS.(WW{i}) = eval(WW{i});

  if isempty(eval(WW{i})) && isfield(def,WW{i}) && ~isobject(eval(WW{i}))
    INPUTS.(WW{i}) = def.(WW{i});
  end
  
  if ~isfield(def,WW{i}) && ~exist('STATIC','var')
    fprintf('WARNING: no default value for input var "%s" set!',WW{i});
    pause(1);
  end
  
end

%set non-existent inputs to defaults
if ALLFIELDNAMES
  WW = fish.helper.allfieldnames(def);
else
  WW = fieldnames(def);
end

for i = 1:length(WW)
  if strcmp(WW{i},'opts') || strcmp(WW{i},'varargin') , continue; end;
  
  if ~exist(WW{i})
    if ALLFIELDNAMES
      eval(sprintf('INPUTS.%s = def.%s;',WW{i},WW{i}));
    else
      INPUTS.(WW{i}) = def.(WW{i}); % same as above, but will result in error
                                    % for higher depth. So keep it. 
    end
  end
end



if ~exist('opts','var')
  opts = [];
end

if exist('varargin','var');
  if mod(length(varargin),2)
    %odd: look for a struc not preceded by a string. this will be the
    %option struc. additional ("par" value) field will override options
    %in opts-struc

    if strcmp(lower(varargin{1}),'help')
      HELP = 1; %for object methods
    else 
      %ordinary input
      i = 1;
      while i<length(varargin)+2
        if i<length(varargin) && ischar(varargin{i})
          i = i+2;
        elseif isstruct(varargin{i})
          opts = varargin{i};
          varargin(i) = [];
          break
        else
          error(['ERROR: wrong input options structure: (''option'',value)-pairs and/or  ',...
                 ' _one_ option structue opts (otps has lower precedence).']);
        end
      end
    end
  end
    
  for i = 1:2:length(varargin)
    if ischar(varargin{i})
      if ALLFIELDNAMES
        eval(sprintf('opts.%s = varargin{i+1};',varargin{i}));
      else
        opts.(varargin{i}) = varargin{i+1};
      end
    else
      error('ERROR: wrong input options structure: (''option'',value)-pairs and/or _one_ option structue opts (opts has lower precedence).');
    end
  end
  clear varargin
end

      

%set default options
if isfield(def,'opts')

  if ~isstruct(def.opts)
    error('expect structure as options');
  end
  
  if ALLFIELDNAMES
    WW = fish.helper.allfieldnames(def.opts);
    if ~isempty(opts)
      WWW = fish.helper.allfieldnames(opts);
    else
      WWW = {};
    end
  else
    WW = fieldnames(def.opts);
    if ~isempty(opts)
      WWW = fieldnames(opts);
    else
      WWW = {};
    end
  end
    
  
  for i = 1:length(WW)
    if ~any(strcmp(WWW,WW{i})) % same as isfield(opts,WW{i})
      if ALLFIELDNAMES
        eval(sprintf('opts.%s = def.opts.%s;',WW{i},WW{i}));
      else
        opts.(WW{i}) = def.opts.(WW{i});
      end
    end
  end
  if (isfield(opts,'STRICT') && opts.STRICT) || ...
        (exist('STRICT','var') && STRICT)
    
    %check for options with no definitions in def.opts (maybe spelling
    %error!) 
    warnfields = {};
    for i = 1:length(WWW)
      if any(strcmp(WWW{i},{'STRICT','VERBOSELEVEL','HELP','ALLFIELDNAMES'})),continue;end;
      if ~any(strcmp(WWW{i},WW))
        warnfields{end+1} = WWW{i};
      end
    end
    
    if ~isempty(warnfields)
      %issue a warning
      if exist('MFILENAME','var')
        fprintf('WARNING %s.parseInput(): These options are invalid: ',MFILENAME);
      else
        fprintf('WARNING parseInput(): These options are invalid: ');
      end
      fprintf(' "%s"',warnfields{:});
      fprintf('\n');
      
      if (isfield(opts,'STRICT') && opts.STRICT>1) || ...
            (exist('STRICT','var') && STRICT>1)
        error('Wrong option types given and STRICT>1!')
      end
    end
    
    clear WWW warnfields STRICT
    if isfield(opts,'STRICT')
      opts = rmfield(opts,'STRICT');
    end
  end
end

if ~isempty(opts)
  INPUTS.opts = opts;
end

clear i WW 

if ~exist('doc','var');
  doc = [];
end

  
if ~evalin('caller','nargin') || HELP


  %display help and defaults 
  if exist('MFILENAME','var')
    fprintf('\n')
    %disp('******************************************************************************')
    help(MFILENAME)
    disp('  ************')
  end
  fprintf('  INPUTS and default values:\n')
  
  %displaydefdoc(def,doc,1)
  SDEF = evalc('disp(def);');
  SDEF = strread(SDEF,'%s','delimiter','\n');
  SDEF = strvcat(SDEF{:});
  EMPTYSTR = char(' '*ones(1,size(SDEF,2)));

  WW = fieldnames(def);
  fprintf('\n')
  for i = 1:length(WW)
    if strcmp(WW{i},'opts'), continue;end
    if isfield(doc,WW{i})
      INFOSTR = doc.(WW{i});
      if ~iscell(INFOSTR)
        INFOSTR = {INFOSTR};
      end
      
      fprintf('\t%d) %s\t%% %s\n',i,SDEF(i,:),INFOSTR{1});
      for j = 2:length(INFOSTR)
        fprintf('\t   %s\t%% .. %s\n',EMPTYSTR,INFOSTR{j});
      end
    else
      fprintf('\t%d) %s\n',i,SDEF(i,:));
    end
  end
  fprintf('\n')
  clear WW SDEF EMPTYSTR INFOSTR
  
  
  
  if ~isempty(opts)
    fprintf('  OPTIONS and default values:\n')

    if ALLFIELDNAMES
      WWW = fish.helper.allfieldnames(opts);
    else
      WWW = fieldnames(opts);
    end

    %displaydefdoc(opts,doc)
    WW = [];
    for i = 1:length(WWW)
      INFOSTR=regexp(evalc(sprintf('WW.info = opts.%s',WWW{i})),'(:\W.*)\n\n','tokens');
      SDEF{i} = [sprintf('%s ',WWW{i}),INFOSTR{1}{1}];
    end
    
    SDEF = strvcat(SDEF{:});
    EMPTYSTR = char(' '*ones(1,size(SDEF,2)));
    fprintf('\n');
    if ALLFIELDNAMES
      WW = fish.helper.allfieldnames(doc);
    else
      if isempty(doc)
        WW = {};
      else
        WW = fieldnames(doc);
      end
    end
    
    for i = 1:length(WWW)
      if any(strcmp(WW,WWW{i}))
        eval(sprintf('INFOSTR = doc.%s;',WWW{i}));
        if ~iscell(INFOSTR)
          INFOSTR = {INFOSTR};
        end
        fprintf('   %s %% %s\n',SDEF(i,:),INFOSTR{1});
        for j = 2:length(INFOSTR)
          %multiple rows
          fprintf('   %s %% .. %s\n',EMPTYSTR,INFOSTR{j});
        end
      else
        fprintf('   %s\n',SDEF(i,:));
      end
    end
    fprintf('\n')
    clear SDEF EMPTYSTR INFOSTR 
  end


  HELP = 1; %for the program to stop

else
  
% $$$   WW = whos;
% $$$   WW = {WW.name};
% $$$   %make nice struct for displaying
% $$$   for i = 1:length(WW)
% $$$     if strcmp(WW{i},'def'),continue;end;
% $$$     INPUTS.(WW{i}) = eval(WW{i});
% $$$   end

  %reduce verbosity if function is called within a function if not
  %explicitly set
  if ~(isfield(opts,'VERBOSELEVEL') && opts.VERBOSELEVEL==0) && ...
        ~(exist('VERBOSELEVEL','var') && VERBOSELEVEL==0)
    if size(dbstack,1)>2 % call within a function
      VERBOSELEVEL = 0;
    end
  end
  
  
  %verbose
  if ~(isfield(opts,'VERBOSELEVEL') && opts.VERBOSELEVEL==0) && ...
        ~(exist('VERBOSELEVEL','var') && VERBOSELEVEL==0)
    if exist('MFILENAME','var')
      fprintf('\nINPUTS to %s:\n\n',MFILENAME);
    else
      fprintf('\nINPUTS:\n\n')
    end

    %displaydefdoc(INPUTS,doc,1);     %should not depend on external files
    SDEF = evalc('disp(INPUTS);');
    SDEF = strread(SDEF,'%s','delimiter','\n');
    SDEF = strvcat(SDEF{:});
    EMPTYSTR = char(' '*ones(1,size(SDEF,2)));
    WW = fieldnames(INPUTS);
    %fprintf('\n')
    for i = 1:length(WW)
      if strcmp(WW{i},'opts');continue;end;

      if isfield(doc,WW{i})
        INFOSTR = doc.(WW{i});
        if ~iscell(INFOSTR)
          INFOSTR = {INFOSTR};
        end
        fprintf('\t%d) %s %% %s\n',i,SDEF(i,:),INFOSTR{1});
        for j = 2:length(INFOSTR)
          fprintf('\t   %s %% .. %s\n',EMPTYSTR,INFOSTR{j});
        end
      else
        fprintf('\t%d) %s\n',i,SDEF(i,:));
      end
    end
    fprintf('\n')
    clear WW SDEF EMPTYSTR INFOSTR
    
    if ~isempty(opts)

      if ALLFIELDNAMES
        WWW = allfieldnames(opts);
      else
        WWW = fieldnames(opts);
      end

      fprintf('  OPTIONS:\n\n')
      %displaydefdoc(INPUTS.opts,doc)
      WW = [];
      for i = 1:length(WWW)
        INFOSTR=regexp(evalc(sprintf('WW.info = opts.%s',WWW{i})),'(:\W.*)\n\n','tokens');
        SDEF{i} = [sprintf('%s ',WWW{i}),INFOSTR{1}{1}];
      end
      SDEF = strvcat(SDEF{:});
      EMPTYSTR = char(' '*ones(1,size(SDEF,2)));
      %fprintf('\n')
      if ALLFIELDNAMES
        WW = fish.helper.allfieldnames(doc);
      else
        if isempty(doc)
          WW = {};
        else
          WW = fieldnames(doc);
        end
      end
      for i = 1:length(WWW)
        if any(strcmp(WW,WWW{i}))
          eval(sprintf('INFOSTR = doc.%s;',WWW{i}));
          if ~iscell(INFOSTR)
            INFOSTR = {INFOSTR};
          end
          fprintf('   %s %% %s\n',SDEF(i,:),INFOSTR{1});
          for j = 2:length(INFOSTR)
            fprintf('   %s %%  .. %s\n',EMPTYSTR,INFOSTR{j});
          end
        else
          fprintf('   %s\n',SDEF(i,:));
        end
      end
      fprintf('\n')
      clear SDEF EMPTYSTR INFOSTR
    end
  end

  if exist('MFILENAME')
    INPUTS.MFILENAME = MFILENAME;
    
    %order fields
    WW = length(fieldnames(INPUTS));
    INPUTS = orderfields(INPUTS,[WW,1:WW-1]);
    clear WW
  end

  %append doc to INPUTS var
  INPUTS.doc = doc;

  
  %INPUTS.opts = opts;
  clear i j WW MFILENAME VERBOSELEVEL def doc WWW ALLFIELDNAMES
  
  %put everything in workspace
  if ~exist('STATIC','var')
    for WW = fieldnames(INPUTS)';
      if ~any(strcmp(WW{1},{'MFILENAME','doc'}))
        eval([WW{1} '= INPUTS.(WW{1});']);
      end
    end
  end
  clear WW 
  
  
end

