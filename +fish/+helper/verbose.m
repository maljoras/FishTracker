function verbose(str,varargin)

  global VERBOSELEVEL
  global VERBOSEDIARY  

  
  if ~VERBOSELEVEL
    return
  end

  if any(VERBOSEDIARY)
    if ischar(VERBOSEDIARY)
      fid = fopen(VERBOSEDIARY,'a+');
    end
  else
    fid = 1;
  end
  
  
  nostack = 0;
  
  if str(1)=='0'
   str = str(2:end);
   nostack = 1;
  end

  if isempty(findstr(str,'\r'))
    str = [str  '\n']; % return;
  end

  if fid~=1
    str(findstr(str,'\r')+1) = 'n';
  end
  
  if ispc()
    persistent LASTDISPLENGTH
    if LASTDISPLENGTH
      for i = 1:length(LASTDISPLENGTH)
        fprintf('\b');
      end
    end
    
    idx = findstr(str,'\r');
    if ~isempty(idx)
      LASTDISPLENGTH = idx(end);
      str([idx(end),idx(end)+1]) = [];
    else
      LASTDISPLENGTH = 0;
    end
  end
  
  d = dbstack;
  if length(d)==1 || nostack
    fprintf(fid,['[%s]:  ' str ],datestr(now),varargin{:});
  else
    fprintf(fid,['[%s,%s]:  ' str],datestr(now),d(2).name,varargin{:});
  end
  
  if any(VERBOSEDIARY)
    fclose(fid);
  end
  