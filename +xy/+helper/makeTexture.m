function texture = makeTexture(typestr,dim,varargin)
% [TEXTURE] = MAKETEXTURE(TYPESTR,DIM,VARARGIN) generates texture for
% stimulus for Zoltan's spatial info dynamics experiment.
%
% 21.5.2006, MJR,  malte@igi.tu-graz.ac.at
% 7.1.2011, freqs changed to visual degrees for illu bars 
  
  VERBOSELEVEL = 0;
  
  def.typestr = 'gratings';
  
  def.dim = [300,300];
  
  def.opts.freq = 10;    
  
  doc.freq = 'in cyc/visual degree';
  
  def.opts.ori = pi/4;
  def.opts.phase = 0;

  def.opts.carrierori = 0;
  def.opts.carrierfreq = 10;
  def.opts.carrierwidth = 1;
  doc.carrierwidth = 'in degrees';
  
  def.opts.contourfreq = 2;
  def.opts.contourori = pi/2;
  def.opts.contourphase = 0;

  def.opts.visspace = [1,1];
  doc.visspace = 'in visual degrees';

  def.opts.barsthres = 0.815;
  doc.barsthres = 'threshold to get b&w bars from gratings (for illubars)';
  
  
  def.opts.gamma =  3; 
  doc.gamma = 'elongation in freqspace in case of orientation';
  
  def.opts.cutoff = 0.8;
  doc.cutoff = 'cutoff for gaussian in freqspace';
  
  
  def.opts.contrast = 1;

  def.opts.filter = 'butter(4,0.20,''low'')';
  def.opts.filter2d = {0.2,'low',4};
  
  
  %mirroraxis
  def.opts.ndots = 10;
  def.opts.dfreq = 10; % size of dots
  def.opts.gap = 0.1;
  def.opts.border = 0.1;
  def.opts.alignedif = 0;

  
  MFILENAME = mfilename;
  xy.helper.parseInputs;
  if HELP;return;end
  
  
  %convert freq into cyc/unit
  % given in cyc per visual degrees [approximate (exact if square)]
  carrierfreq = opts.carrierfreq*mean(opts.visspace);
  contourfreq = opts.contourfreq*mean(opts.visspace);

  
  carrierori = opts.carrierori;
  
  switch lower(typestr)
   case 'gratings'

    [X,Y] = meshgrid((0:dim(1)-1)/dim(1)*opts.visspace(1),(0:dim(2)-1)/dim(2)*opts.visspace(2));
    
    texture = sin(2*pi*opts.freq* (cos(opts.ori).*X + sin(opts.ori).*Y) + opts.phase);

   case 'gratingscircle'
    % as gratings but in circle with border
    
    [X,Y] = meshgrid((0:dim(1)-1)/dim(1),(0:dim(2)-1)/dim(2));
    texture = sin(2*pi*opts.freq* (cos(opts.ori).*X + sin(opts.ori).*Y) + opts.phase);

    msk =(X-0.5).^2 + (Y-0.5).^2 <= (0.5-opts.border).^2;
    texture = texture.*msk;
    
    case 'gaborcircle'
    % as gratings but in circle with border
    
    [X,Y] = meshgrid((0:dim(1)-1)/dim(1),(0:dim(2)-1)/dim(2));
    texture = sin(2*pi*opts.freq* (cos(opts.ori).*X + sin(opts.ori).*Y) + opts.phase);

    msk =exp(-((X-0.5).^2 + (Y-0.5).^2)/2/(0.5-opts.border*2).^2);
    texture = texture.*msk;
    
    
   case 'sparsedotnoise'
    texture = rand(dim);
    texture = texture.*(texture>0.95) - ones(size(texture)).*(texture<=0.95);

   case 'noise2d'
    texture = randn(dim);
    if ~isempty(opts.filter2d)
      if iscell(opts.filter2d)
        texture = filtbutter2(texture,opts.filter2d{:});
      else
        texture = feval(opts.filter2d,texture);
      end
    end

    
   case 'doubledots'


    gap = opts.gap/mean(opts.visspace);
    
    dots = rand(opts.ndots,2)*sqrt(2);

    dots = [dots;bsxfun(@plus,dots,[gap,0])];

    %orientation
    dots = (dots-sqrt(2)/2)*R2d(opts.ori)' + sqrt(2)/2;
    dots(any(dots<0,2) | any(dots>=1,2),:) = [];
    
    
    idots = xy.helper.s2i(dim,floor(bsxfun(@times,dots,dim)) + 1);
    
    texture = zeros(dim);
    texture(idots) = 1;

    
    %filter
    r = 1./(opts.freq*opts.visspace)/2/2;
    
    n = ceil(r.*dim);
    [X,Y] = meshgrid(0:2*n(1),0:2*n(2));

    H = double((((X-n(1))/n(1)).^2 + ((Y-n(2))/n(2)).^2) <=1);
    
    
    texture = conv2(texture,H,'same');
    
    texture = min(texture,1);
    
    
    
   case 'mirroraxis'
    % some mirrored points (ndots) with size of freq, gap between and oriented
    % somewhat (ori).
    
    % border on all sites is border/2
    % ndots will be mirrored on other half 

    if opts.border>0.5 || opts.border<0
      error('border between 0..0.5')
    end
    
    if (0.5-opts.gap/2-opts.border) <=0
      error('gap has to be smaller than 0.5-border')
    end

    dots = rand(opts.ndots,2);
    
    if opts.alignedif
      dots(:,1) = dots(1,1);
    end
    
    
    mrdots = dots;
    mrdots(:,1) = 1-mrdots(:,1);
    dots(:,2) = dots(:,2)*(1-2*opts.border)+opts.border;
    dots(:,1) = dots(:,1)*(0.5-opts.gap/2-opts.border) + opts.border;

    %mirror
    mrdots(:,2) = mrdots(:,2)*(1-2*opts.border)+opts.border;
    mrdots(:,1) = mrdots(:,1)*(0.5-opts.gap/2-opts.border)+0.5+opts.gap/2;
    

    dots = [dots;mrdots];
    
    %orientation
    dots = (dots-0.5)*R2d(opts.ori)' + 0.5;
    dots(any(dots<0,2) | any(dots>=1,2),:) = [];
    

    
    idots = xy.helper.s2i(dim,floor(bsxfun(@times,dots,dim)) + 1);
    
    texture = zeros(dim);
    texture(idots) = 1;

    
    %filter
    r = 1./(opts.dfreq*opts.visspace)/2/2;
    
    n = ceil(r.*dim);
    [X,Y] = meshgrid(0:2*n(1),0:2*n(2));

    H = double((((X-n(1))/n(1)).^2 + ((Y-n(2))/n(2)).^2) <=1);
    
    
    texture = conv2(texture,H,'same');
    
    texture = min(texture,1);
   
    
    
    
   case 'noise'
    texture = randn(dim);

    p0 = 1:length(dim);
    bg = texture;
    if ~isempty(opts.filter)
      for i = 1:length(dim)
        if iscell(opts.filter)
          [b,a] = eval(opts.filter{length(dim)-i+1});
        else
          [b,a] = eval(opts.filter);
        end
        
        p = circshift(p0,[0,i]);
        bg = permute(bg,p);
        bg = filter(b,a,bg(:));
        bg = reshape(bg,dim(p));
        bg = ipermute(bg,p);
      end
    end
    
    texture = bg;
    
   case 'orinoise'
    
    texture = subGetOriNoise(dim,opts);
   
   case {'illusory_gratings','illusory','illu'}
    %illusory contour with gratings as cue

    phasediff = pi;
    
    [X,Y] = ndgrid((0:dim(1)-1)/dim(1),(0:dim(2)-1)/dim(2));
    
    grat1 = sin(2*pi*carrierfreq* (cos(carrierori).*X + sin(carrierori).*Y) + opts.phase);
    grat2 = sin(2*pi*carrierfreq* (cos(carrierori).*X + sin(carrierori).*Y) + opts.phase+phasediff);
    
    contour = sin(2*pi*contourfreq* (cos(opts.contourori).*X + sin(opts.contourori).*Y) + opts.contourphase);
    
    texture = grat1.*(contour>=0) + grat2.*(contour<0);
    
    
   case {'carrierbar','carrierbars','bars','bar'}
    [X,Y] = meshgrid((0:dim(1)-1)/dim(1),(0:dim(2)-1)/dim(2));
    texture = sin(2*pi*carrierfreq* (cos(opts.carrierori).*X + sin(opts.carrierori).*Y) + opts.phase);

    w = opts.carrierwidth/mean(opts.visspace);

    if w>1/carrierfreq
      xy.helper.verbose('WARNING: carrierwidth is specified too large')
    end
    
    x1 = cos(pi*carrierfreq*w);
    texture = double(texture>x1);


    

   case {'illusory_bars','illubar','illubars'}
    %illusory contour with contrasty bars as cue
   
    phasediff = pi;
    
    ntimes = 1;
    dimb = dim*ntimes;
    [X,Y] = ndgrid((0:dimb(1)-1)/dimb(1),(0:dimb(2)-1)/dimb(2));

    grat1 = sin(2*pi*carrierfreq* (cos(carrierori).*X + sin(carrierori).*Y) + opts.phase);
    grat2 = sin(2*pi*carrierfreq* (cos(carrierori).*X + sin(carrierori).*Y) + opts.phase+phasediff);
    
    contour = sin(2*pi*contourfreq* (cos(opts.contourori).*X + sin(opts.contourori).*Y) + opts.contourphase);
    
    texture = grat1.*(contour>=0) + grat2.*(contour<0);
    
    texture = (texture.*(texture>0)).^4; % to get less "oebertoene" ?!
    texture = double(texture>=opts.barsthres);

    
   case 'f2ori'
    % phase as 'orinoise' but with 1/f^2 spectrum in each direction

    %texture = subGetOriNoise(dim,opts);
    p = 2; %power
    theta = opts.ori;
    gamma = opts.gamma;
    
    [f1,f2] = freqspace(dim,'meshgrid');
    z = [find(f1==0 & f2==0)];
    f1(z) = 1; %mean power
    f2(z) = 1; %irrelevant

    xx =   f1 * cos(theta) + f2*sin(theta);
    yy = - f1 * sin(theta) + f2*cos(theta);

    pow = 1./sqrt((abs(xx).^2 + gamma^2*abs(yy).^2)).^p;
    pow((f1.^2+f2.^2)<opts.cutoff) = 0;
    pow = pow/max(pow(:));
    
    
    for i = 1:2
      texture = randn(dim);
      %imagesc(texture);
      
      F = fxyT2(texture);

      F = F.*pow;
    
      %F(z) = 0;
      tex{i} = abs(ifxyT2(F)).^2;
    end
    
    texture = tex{1}-tex{2}; % difference texture
    texture = texture/std(texture(:));

   otherwise 
    error('dont know texture type');
  end
  

  %
  %texture = reshape(mapstd(texture(:)'),dim);
  
  %  m = texture-min(texture(:));
  %  texture = 2*m/max(m(:))-1;
  

  if opts.contrast>0
    %contrast
    texture = tanh(opts.contrast*texture);
    texture = (texture+1)/2;
  end
  
  
  
  
  
  
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function texture = subGetOriNoise(dim,opts)  

   texture = randn(dim);
    
   fn =100;
    
   [f1,f2] = freqspace(fn,'meshgrid');
   
   
   cutoff = opts.cutoff;
   theta = opts.ori;
   
   sigma = 0.3; %3sigma
   gamma = opts.gamma; %elongation
   
   xx =   f1 * cos(theta) + f2*sin(theta);
   yy = - f1 * sin(theta) + f2*cos(theta);
   
   resp = exp(-(xx.^2+gamma^2*yy.^2)/(2*sigma^2))>cutoff;
   
   %imagesc(resp)
  % pause
   h = fwind1(resp,hamming(fn));
   
   
   texture = real(imfilter(texture,h));
    
