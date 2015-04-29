function verbose(str,varargin)

  global VERBOSELEVEL
  
  if ~VERBOSELEVEL
    return
  end

  nostack = 0;
  
  if str(1)=='0'
   str = str(2:end);
   nostack = 1;
  end

  if isempty(findstr(str,'\r'))
    str = [str  '\n']; % return;
  end
  
  d = dbstack;
  if length(d)==1 || nostack
    fprintf(['[%s]:  ' str ],datestr(now),varargin{:});
  else
    fprintf(['[%s,%s]:  ' str],datestr(now),d(2).name,varargin{:});
  end
  