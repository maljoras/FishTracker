function [varargout]= errorbarpatch(x,y,u,varargin)
%  H = ERRORBARLINE(X,Y,U,...) plots like errorbar a line and an errorbar
   
  

  %just use plot
  h1 = plot(x,y,varargin{:});
  hold on

  
  if (min(size(y))>1) && (min(size(x))==1)
    if size(x,1)==1
      x = x';
    end
    x = repmat(x,[1,size(y,2)]);
  end

  
  if min(size(y))==1
    if size(y,1) == 1
      y = y';
    end
    
    if size(x,1) == 1
      x = x';
    end

    if size(u,1) == 1
      u = permute(u,[2,1,3]);
    end
  end
  
  
  for i = 1:size(y,2)
    %for multiple data
    %just make loop
    
    if size(u,3)>1
      % directly given confidence
      tmpy = [u(:,i,1), u(:,i,2)];
    else
      tmpy = [y(:,i)-u(:,i), y(:,i)+u(:,i)];
    end
    
    tmpx = [x(:,i) , x(:,i)];

    c = get(h1(i),'Color');
    c = c + 0.8*(1-c);
    h2(i) = patch([tmpx(:,1)' tmpx(end:-1:1,2)'],[tmpy(:,2)',tmpy(end:-1:1,1)'],c,'EdgeColor',c,'FaceColor',c);
  end
  
  
  %put patches in the back.
  ch = get(gca,'children');
  p = findobj(gca,'type','patch');
  idx = find(ismember(ch,p));
  others = setdiff(1:length(ch),idx);
  ch = ch([others,idx']);
  set(gca,'children',ch);
  
  
  %delete(h1);
  %h1 = plot(x,y,varargin{:});

  
  
  
  varargout{1} = h1;
  varargout{2} = h2;
  
  

  
 