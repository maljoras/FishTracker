function a2 = submargin(a1,varargin);
%  A1 = SUBMARGIN(A,X,Y,...) plots a second plot below the actual plot
%  and reduces the size of the first axis. 

  location = 1; %or 0 for lexyT;
  SPACE = 0.20; % in percent of the original plot
  MARGIN = 0.05; %margin between plots
  
  PERCENT = 1;
  SWAP = 0;  % whether to swap lexyT right, or bottom, top
  
  args = varargin;
  for i = length(args)-1:-2:1
    if strcmp(lower(args{i}),'location')
      location = args{i+1};
      args{i} = [];
      args{i} = [];
    end
    
    if strcmp(lower(args{i}),'space')
      SPACE = args{i+1};
      args{i} = [];
      args{i} = [];
    end

    if strcmp(lower(args{i}),'percent')
      PERCENT = args{i+1};
      args{i} = [];
      args{i} = [];
    end
    
    if strcmp(lower(args{i}),'margin')
      MARGIN = args{i+1};
      args{i} = [];
      args{i} = [];
    end

    if strcmp(lower(args{i}),'swap')
      SWAP = args{i+1};
      args{i} = [];
      args{i} = [];
    end

    
  end
  
  
  if ischar('location')
    switch lower(location)
     case 'north'
      location = 1;
      SWAP = 1;
     case 'south'
      location = 1;
      SWAP = 0;
     case 'east'
      location = 0;
      SWAP = 1;
     case 'west'
      location = 0;
      SWAP = 0;
     otherwise
      %default
      location = 1;
      SWAP = 0;
      
    end
  end
  
    
      
      
  
  
  
  p1 = get(a1,'position');
  p2 = p1;  

  ind1 = 1+location;
  ind2 = 3+location;

  set(a1,'xticklabel',[]);
  
    
  %plot below
  if PERCENT
    wy = SPACE*p1(ind2) - MARGIN/2;
    p1(ind1) = wy+ p1(ind1) + MARGIN;
    p1(ind2) = (1-SPACE)*p1(ind2) - MARGIN/2;
  else
    %absolut values
    wy = SPACE; 
    p1(ind1) = SPACE + p1(ind1) + MARGIN;
    p1(ind2) = p1(ind2) - SPACE - MARGIN;
  end
  
    
  set(a1,'position',p1);
  p2(ind2) = wy;
  a2 = axes;
  set(a2,'position',p2);

  
  
  if location 
    xl = get(a1,'xlim');
    set(a2,'xlim',xl);
  else
    %on side
    yl = get(a1,'ylim');
    set(a2,'ylim',yl);
  end
  

  if SWAP
    p1 = get(a1,'position');
    p2 = get(a2,'position');

    p3 = p1;
    if location % bottom top
      p1(2) = p2(2);
      p2(2) = p3(4)+p3(2) - p2(4);
      set(a2,'xaxisloc','top')
    else
      p1(1) = p2(1);
      p2(1) = p3(3)+p3(1) - p2(3);
      set(a2,'yaxisloc','right')
    end
    
    set(a1,'position',p1);
    set(a2,'position',p2);
  
  end
  
  

  
