function a =subsubplot(varargin)
% AX = SUBSUBPLOT(ORGR1,ORGR2,ORGI,R1,R2,I) plots subplots of size [r1,r2]
% into the subplot araea given by orgz
  
  
  phandle = gcf;
  if nargin==7 && ishandle(varargin{1}) 
    phandle = varargin{1};
    args = varargin(2:7);
  else
    args = varargin(1:6);
  end
  
  opts.margin = 0.1;
       
  [orgr1,orgr2,orgi,r1,r2,i] = deal(args{:});


  % get the position as in subplot
  nRows = orgr1;
  nCols = orgr2;
  plotId = orgi;
  def_pos = get(phandle, 'DefaultAxesPosition');
  % This is the percent offset from the subplot grid of the plotbox.
  inset = [.2, .18, .04, .1]; % [left bottom right top]
  
  row = (nRows - 1) - fix((plotId - 1) / nCols);
  col = rem(plotId - 1, nCols);

  % compute outerposition and insets relative to figure bounds
  rw = max(row) - min(row) + 1;
  cw = max(col) - min(col) + 1;
  width = def_pos(3) / (nCols - inset(1) - inset(3));
  height = def_pos(4) / (nRows - inset(2) - inset(4));
  inset = inset .* [width, height, width, height];
  outerpos = [def_pos(1) + min(col) * width - inset(1), ...
              def_pos(2) + min(row) * height - inset(2), ...
              width * cw, height * rw];
  
  % compute inner position
  position = [outerpos(1 : 2) + inset(1 : 2), ...
              outerpos(3 : 4) - inset(1 : 2) - inset(3 : 4)];

  
  p1 = position;

  %% compute new position
  
  x = p1(3);
  y = p1(4);
  
  

  width = (x-(r2>1)*opts.margin*x)/r2;
  height = (y-(r1>1)*opts.margin*y)/r1;
  

  wmargin = opts.margin*x/(max(r2,2)-1);
  hmargin = opts.margin*y/(max(r1,2)-1);
  
  [k,j] = ind2sub([r2,r1],i);
  
  xpos = (width+wmargin)*(k-1) + p1(1);
  ypos = (height+hmargin)*(r1-j) + p1(2);
  
  p = [xpos,ypos,width,height];
  
  a =axes('parent',phandle);
  set(a,'position',p);
  
  