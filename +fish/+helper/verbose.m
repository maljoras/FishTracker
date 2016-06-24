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
      for i = 1:LASTDISPLENGTH
        fprintf('\b');
      end
    end
    
    idx = findstr(str,'\r');
    if ~isempty(idx)
      LASTDISPLENGTH = 1;
      str([idx(end),idx(end)+1]) = [];
    else
      LASTDISPLENGTH = 0;
    end
  end
  
  d = dbstack;
  if length(d)==1 || nostack
    strout = sprintf(['[%s]:  ' str ],datestr(now),varargin{:})
    fprintf(fid,strout);
  else
    strout = sprintf(['[%s,%s]:  ' str],datestr(now),d(2).name,varargin{:});
    fprintf(fid,strout);
  end
  
  if ispc() && LASTDISPLENGTH
    LASTDISPLENGTH = length(strout);
  end
  
  if any(VERBOSEDIARY)
    fclose(fid);
  end
  