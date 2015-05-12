classdef MyForegroundDetector < handle;
  
 
    
  properties 
    
    thres = 'auto';
  
    fgmsk = 1;
    bgimg = 0;

    mtau = 100;
    dtau = 0;
    dratio = 0.25;
    inverse = 0;
    plotif = 1;

    %   manuallyAdjustThres = true;
    
    rgbchannel = [];
    map = [];
  end

  properties (SetAccess = private);

    mframe = [];
    frame = [];
    nauto = 30;
    framecounter = 0;
    bwimg = [];
    msk = [];

    dframe = [];
    fv = [];
    mu = [];
  end
  
  
  methods (Access=private)
    
    
    function frame = convertFrame(self,oframe)

      frame = oframe;

      if size(frame,3)==3 
        if ~isempty(self.fv)
          %frame = imadjust(oframe,[.4 .4 0.2; .6 .5 .3],[],1);
          frame = sum(bsxfun(@times,bsxfun(@minus,oframe,self.mu),self.fv),3);
        else
          if ~isempty(self.map)
            frame = rgb2ind(oframe,self.map,'nodither');
            
          else
            if isempty(self.rgbchannel)
              frame = rgb2gray(frame);   
            else
              frame = frame(:,:,self.rgbchannel);   
            end
          end
        end
        %frame = mean(frame,3);      
      end
      if ~isfloat(frame)
        frame = double(frame)/255;      
      end
    end
    
    
    function [thres,objlength,objwidth, fv, mu,bw] = getAutoThreshold(self,frame,oframe,nObjects)

      % minsize = 5;
      mf = mean(frame(:));

      if self.inverse
        thresarr = shiftdim(linspace(max(frame(:)),mf,self.nauto),-1);
        thresarr =thresarr(end:-1:1);
        bwimg = bsxfun(@ge,frame,thresarr(2:end-1));
      else
        thresarr = shiftdim(linspace(min(frame(:)),mf,self.nauto),-1);
        bwimg = bsxfun(@le,frame,thresarr(2:end-1));
      end
      
      bwimg = bsxfun(@and,bwimg,self.fgmsk);
      
      largest = nan(1,size(bwimg,3));
      nextlargest = nan(1,size(bwimg,3));
      %msize = nan(1,size(bwimg,3));
      %nspots = nan(1,size(bwimg,3));
      
      for i = 1:size(bwimg,3)
        tmp = bwconncomp(bwimg(:,:,i),8);;
        rp = regionprops(tmp,'MajorAxisLength','MinorAxisLength');
        
        
        %sz = cellfun('length',tmp.PixelIdxList);
        sz = cat(1,rp.MajorAxisLength);
        %sz(sz<minsize) = [];
        %nspots(i) = length(sz);
        
        % get the largest objects
        ssz = sort(sz,'descend');

        if length(ssz)
          largest(i) = mean(ssz(1:min(end,1)));
          nextlargest(i) = mean(ssz(min(end,2):min(end,nObjects)));
        end
        
      end
      ratio = nextlargest./largest;
      conds = largest<0.2*max(size(bwimg)) & nextlargest>10;%minimal 10 pixels

      if ~sum(conds)
        imagesc(frame)
        error('Cannot detect objects. Change parameter settings');
      end

      mr = min(max(ratio(conds)),0.9);
      % largest axis should be less than a 5th of the movie size.
      thresidx = find(nextlargest./largest>=mr & conds,1,'last');

      % look for the size of the objects
      bw = bwconncomp(bwimg(:,:,thresidx),4);;
      sz = cellfun('length',bw.PixelIdxList);
      [ssz, idx]= sort(sz,'descend');
      idx = idx(1:min(nObjects,end));
      bw.NumObjects = length(idx);
      bw.PixelIdxList = bw.PixelIdxList(idx);
      rp = regionprops(bw,'MajorAxisLength','MinorAxisLength');
      
      objlength = round(mean(cat(1,rp.MajorAxisLength)));
      objwidth = round(mean(cat(1,rp.MinorAxisLength)));
      
      % get a good rgb to color conversion
      objidx = cat(1,bw.PixelIdxList{:});
      %randidx = ceil(rand(min(length(objidx)*10,prod(bw.ImageSize)),1)*prod(bw.ImageSize));
      objij = i2s(bw.ImageSize,objidx);
      randij = ceil(bsxfun(@plus,objlength*randn([size(objij),10]),objij));
      randidx = s2i(bw.ImageSize,reshape(permute(randij,[1,3,2]),[],2));
      randidx = max(min(randidx,prod(bw.ImageSize)),1);
      randidx = setdiff(randidx,objidx);

      oframe1 = double(reshape(oframe,[],3));
      col1 = oframe1(objidx,:);
      col2 = oframe1(randidx,:);
      [fv,~,proj,r] = lfd([col1;col2],[ones(size(col1,1),1);zeros(size(col2,1),1)],1,0);
      
      frame1 = sum(bsxfun(@times,bsxfun(@minus,oframe,shiftdim(r.mu,-1)),shiftdim(fv,-2)),3);

      if mean(frame1(objidx))>mean(frame1(randidx))
        fv = -fv;
        frame1 = sum(bsxfun(@times,bsxfun(@minus,oframe,shiftdim(r.mu,-1)),shiftdim(fv,-2)),3);
      end
      
      fv = shiftdim(fv,-2);
      mu = shiftdim(r.mu,-1);
      
      
% $$$       if self.manuallyAdjustThres
% $$$         keyboard
% $$$         
% $$$         
% $$$       end


      thres1 = mean(frame1(randidx)) - 2.5*std(frame1(randidx));
      thres2 = mean(frame1(objidx)) + 3*std(frame1(objidx));

      if thres2<thres1
        thres = thres1 + (thres2-thres1)/2;
      else
        thres = max(thres1,mean(frame1(objidx)));
      end
      
      %d = abs([diff(nspots,2),0,0]);
      %idx = nspots>0 & msize<self.maxsize;
      %thres1 = thresarr(1+find(idx & d==min(d(idx)),1,'last'));

      %thres1 = thresarr(thresidx); 
      % set the property
      %thres = thres1-mf;

      
        
      if self.dtau
        thres = (1-self.dratio)*thres; % approx;
      end
      
    end

    
    function bwimg = getBWImage(self,frame)

      frame = frame - self.mframe;
      
      % update and subtract mean frame
      mtau = min(self.framecounter,self.mtau);
      self.mframe = (mtau-1)/mtau * self.mframe + 1/mtau*frame; % running avaerge



      
      if self.dtau>1
        f1 = frame;
        frame = (1-self.dratio)*frame + (self.dratio)*self.dframe;
        % difference frame
        self.dframe = (self.dtau-1)/self.dtau * self.dframe + 1/self.dtau*f1; % running avaerge
        %figure(10);
        %imagesc(frame);
        %drawnow;
      end
      
      
      % bw image
      if self.inverse
        bwimg = (frame >= self.thres) & self.msk;
      else
        bwimg = (frame <= self.thres) & self.msk;
      end
      
      %save info
      self.bwimg = bwimg;
      self.frame = frame;
    end
    
  end
  
  % public methods
  methods    
    
    function self = MyForegroundDetector(varargin) % constructor
      
      self = self@handle();

      nargs = length(varargin);
      if nargs>0 && mod(nargs,2)
        error('expected arguments of the type ("PropName",pvalue)');
      end
      for i = 1:2:nargs
        if ~ischar(varargin{i})
          error('expected arguments of the type ("PropName",pvalue)');
        else
          self.(varargin{i}) = varargin{i+1};
        end
      end
      
    end
    
    
    function  varargout = init(self,oframe,nObjects)

      nTimes = size(oframe,4);
      for i = 1:nTimes
        frame(:,:,i) = convertFrame(self,oframe(:,:,:,i));
      end
% $$$       % subtract background
% $$$       frame1 = zeros(size(frame));
% $$$       if nTimes>1
% $$$         for i = 1:nTimes
% $$$           others = setdiff(1:nTimes,i);
% $$$           frame1(:,:,i) = frame(:,:,i) - mean(frame(:,:,others),3);
% $$$         end
% $$$         frame = frame1;
% $$$         clear frame1;
% $$$       end
      
      % get adaptive threshold 
      verbose('Get auto threshold and approx. object sizes...')
      for i = 1:nTimes
        verbose('Computing adaptive threshold %d/%d..\r',i,nTimes)
        [thres(i),objh(i),objw(i), fv{i}, mu{i},bw{i}] = getAutoThreshold(self,frame(:,:,i),oframe(:,:,:,i),nObjects);
      end

      thres= mean(thres);
      verbose('Found (mean) threshold: %1.2f',thres);
      objh = round(mean(objh));
      objw = round(mean(objw));
      verbose('Mean object sizes: %d x %d pixels',objh,objw);

      self.fv = mean(cat(5,fv{:}),5);
      self.fv = self.fv/norm(self.fv(:));
      self.mu = mean(cat(5,mu{:}),5);

      if ischar(self.thres)  
        self.thres = thres;
      else
        verbose('However: use the given threshold: %1.2f',self.thres);
      end

      if  numel(self.bgimg)>1
        % pixels likely to be above thres;
        notbgmsk = self.bgimg-mean(self.bgimg(:)) <= self.thres;
      else
        notbgmsk = 1;
      end
      self.msk = self.fgmsk & notbgmsk;

      for i = 1:nTimes
        frame(:,:,i) = convertFrame(self,oframe(:,:,:,i));
      end
      
      self.mframe = 0;
      if nTimes>1
        self.mframe = self.mframe + mean(frame(:,:,2:end),3);
      end

      assert(self.mtau>=1);
      self.dframe = zeros(size(frame));      
      self.frame = frame(:,:,1);
      self.framecounter = 0;
      
      if nargout
        varargout = {objh,objw,self.fv,self.mu};
      end

      if self.plotif
        % plot         
        clf;
        [r1,r2] = getsubplotnumber(nTimes);
        for i = 1:nTimes
          subplot(r1,r2,i);
          L = labelmatrix(bw{i});
          others = setdiff(1:nTimes,i);
          framei = frame(:,:,i)  - mean(frame(:,:,others),3) - self.thres;
          framei = (framei-min(framei(:)))/(max(framei(:))-min(framei(:)));

          rgbLabels = label2rgb(L);
          rgbFrame = repmat(255*max(min(framei,1),0),[1,1,3]);
          msk = repmat(mean(rgbLabels,3)==255,[1,1,3]);
          image(uint8(~msk.*single(rgbLabels) + msk.*rgbFrame));

        end
      end

      
      
    end
    
    
    
    
    function [mask,varargout] = step(self,oframe)
    %  [MASK,FRAME] = STEP(SELF,OFRAME) detects a foreground and returns a BW mask

      if isempty(self.mframe)
        init(self,oframe);
      end

      cframe = convertFrame(self,oframe);
      
      self.framecounter = self.framecounter +1;

      mask = getBWImage(self,cframe);

      %output
      if nargout>1
        varargout{1} = cframe; % RGB converted to intensity
      end
      if nargout>2
        varargout{2} = self.frame; % mean / dframe subtracted
      end
      
    end

    
  end
  
      

end
