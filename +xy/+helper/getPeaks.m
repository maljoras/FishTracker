function peaks = getPeaks(x,varargin);
% gets all peaks along first dimension
% outputs equals input in size and has 1 at peak position and 0 elsewhere
% (each side is zero padded once)
% varargin{1} :  if present, normalize and cutoff below varargin before
%                peak search
% varargin{2} :  if 0: shuts down the ra preprocessing

  if size(x,1)==1;
    x=x';
    trans = 1;
  else
    trans=0;
  end
  
  
  if nargin>1
    x = bsxfun(@minus,x,mean(x));
    x = bsxfun(@rdivide,x,std(x));
    if nargin <= 2 || varargin{2}
      %running average (shifts the peaks quit a bit, one probably should
      %adjust parameters)
      x = xy.helper.movavg(x,3);
    end
    x(x<varargin{1}) = 0;
  end
  
  
  x(end+1,:) = 0; %zero padded 
  x(end+1,:) = 0;
  peaks = (x(2:end,:) - x(1:end-1,:)) > 0;
  peaks = peaks(1:end-1,:) & ~peaks(2:end,:);

  peaks = circshift(peaks,[1,0]);
  
  if trans
    peaks = peaks';
  end
  
