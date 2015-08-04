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
    minAxisLength = 20;
    minMinorAxisLength = 8;
    excludeBorderPercentForAutoThres = 0.02;
    rgbchannel = [];
    map = [];
  end

  properties (SetAccess = private);

    mframe = [];
    frame = [];
    nauto = 30;
    framecounter = 0;
    msk = [];

    dframe = [];
    fv = [];
    mu = [];
  end
  
  
  methods (Access=private)
    
    
    function [frame cframe]= convertFrame(self,oframe,updateif)

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
    

      cframe = frame; % converted without mframe subtraction

      
      if ~isempty(self.mframe) 
        
        frame = frame - self.mframe;
        
        if updateif && self.mtau>0
          % update and subtract mean frame
          mtau = min(self.framecounter,self.mtau);
          
          self.mframe = (mtau-1)/mtau * self.mframe + 1/mtau*cframe; % running avaerge
        end
      end
      
      if self.dtau>1 
        dframe = (1-self.dratio)*frame + (self.dratio)*self.dframe;
        if  updateif
          % difference frame
          self.dframe = (self.dtau-1)/self.dtau * self.dframe + 1/self.dtau*frame; 
        end
        frame = dframe;
      end

      frame = frame - mean(frame(:)); % to avoid light effects no
                                      % present in the mframe
      %self.frame  = frame;

    end
    
    
    
    function [largest,nextlargest,rp_last,bw_last] = getObjectSizes(self,bwimg,mmsk,nObjects)
      largest = nan(size(bwimg,3),2);
      nextlargest = nan(size(bwimg,3),2);
      se = strel('disk',self.minAxisLength);

      for i = 1:size(bwimg,3)
        bwimg(:,:,i) = bwareaopen(bwimg(:,:,i),self.minMinorAxisLength*self.minAxisLength,8);
        bw = bwconncomp(bwimg(:,:,i),8);
        
        rp = regionprops(bw,'MajorAxisLength','MinorAxisLength','PixelIdxList');
        
        sz = cat(1,rp.MajorAxisLength);
        sz_minor = cat(1,rp.MinorAxisLength);

        if length(rp)
          perinmsk = zeros(length(rp),1);
          for ii = 1:length(rp)
            perinmsk(ii) = mean(mmsk(rp(ii).PixelIdxList));
          end
          
          delidx = sz_minor<self.minMinorAxisLength | perinmsk<0.9;
          sz(delidx ) = [];
          sz_minor(delidx) = [];
        end
        
        if length(sz)
          [ssz,sidx] = sort(sz,'descend');
          ssz_minor = sz_minor(sidx);
          % get the largest objects
          if length(ssz)
            largest(i,:) = [ssz(1:min(end,1)),ssz_minor(1:min(end,1))];
            nextlargest(i,:) = [median(ssz(min(end,2):min(end,nObjects))),median(ssz_minor(min(end,2):min(end,nObjects)))];
          end
        end
      end
      rp(delidx) = [];
      rp_last = rp;
      bw_last = bw;
    end
    
      
    
    function [thres,objlength,objwidth, fv, mu,bw] = getAutoThreshold(self,frame,oframe,nObjects)


      % sometimes there is a black/white mask. Ignore these pixels
      mmsk = ~((frame==max(frame(:))) | (frame==min(frame(:))));

      % minsize = 5;
      

      if self.inverse
        thresarr = shiftdim(linspace(median(frame(mmsk)),max(frame(mmsk)),self.nauto),-1);      
        thresarr =thresarr(end:-1:1);
        bwimg = gather(bsxfun(@ge,frame,thresarr(2:end-1)));
      else
        thresarr = shiftdim(linspace(min(frame(mmsk)),median(frame(mmsk)),self.nauto),-1);      
        bwimg = gather(bsxfun(@le,frame,thresarr(2:end-1)));
      end
      
      borderPixels = round(size(frame)*min(self.excludeBorderPercentForAutoThres,0.4));

      bwimg = bsxfun(@and,bsxfun(@and,bwimg,self.fgmsk),mmsk);
      if borderPixels>0  
        bwimg(1:borderPixels,:,:) = 0;
        bwimg(:,1:borderPixels,:) = 0;
        bwimg(end-borderPixels+1:end,:,:) = 0;
        bwimg(:,end-borderPixels+1:end,:) = 0;
      end
      
      [largest,nextlargest] = self.getObjectSizes(bwimg,mmsk,nObjects);
      ratio = nextlargest(:,1)./largest(:,1);
      conds = largest(:,1)<0.2*max(size(bwimg)) & nextlargest(:,1)>self.minAxisLength;%minimal 10 pixels


      if ~sum(conds)
        warning('Cannot detect objects automatically. Better switch to manual...');
        objlength = NaN;
        objwidth = NaN;
        thres = NaN;
        fv = [];
        mu = [];
        bw = [];
        return
      end

      mr = min(max(ratio(conds)),0.9);
      % largest axis should be less than a 5th of the movie size.
      thresidx = find(nextlargest(:,1)./largest(:,1)>=mr & conds,1,'last');

      % look for the size of the objects
      [largest,nextlargest, rp, bw] = self.getObjectSizes(bwimg(:,:,thresidx),mmsk,nObjects);
      if nObjects==1
        objlength = round(largest(:,1));
        objwidth = round(largest(:,2));
      else
        objlength = round(nextlargest(:,1));
        objwidth = round(nextlargest(:,2));
      end
      
      % get a good rgb to color conversion
      lst = {rp.PixelIdxList};
      objidx = cat(1,lst{:});
      %randidx = ceil(rand(min(length(objidx)*10,prod(bw.ImageSize)),1)*prod(bw.ImageSize));
      objij = i2s(size(frame),objidx);
      randij = ceil(bsxfun(@plus,objlength*randn([size(objij),10]),objij));
      randidx = s2i(size(frame),reshape(permute(randij,[1,3,2]),[],2));
      randidx = max(min(randidx,prod(size(frame))),1);
      randidx = setdiff(randidx,objidx);

      oframe1 = gather(double(reshape(oframe,[],3)));
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

    
    function [thres,objh,objw,bw] = getAdaptiveAuto(self,oframe,nObjects)
      
      %iteratively
      verbose('Get auto threshold and approx. object sizes...')      
      nTimes = size(oframe,4);
      frame = zeros(size(oframe,1),size(oframe,2),nTimes);

      for i_auto = 1:2
      
        for i = 1:nTimes
          frame(:,:,i) = convertFrame(self,oframe(:,:,:,i),0);
        end
      
        % get adaptive threshold 
        for i = 1:nTimes
          verbose('Computing adaptive threshold %d/%d..\r',i,nTimes)
          [thres(i),objh(i),objw(i), fv{i}, mu{i},bw{i}] = getAutoThreshold(self,frame(:,:,i),oframe(:,:,:,i),nObjects);
        end
        if all(isnan(thres))
          
          return
        end
        %delete nans
        delidx = find(isnan(thres));
        if ~isempty(delidx)
          thres(delidx) = [];
          objh(delidx) = [];
          objw(delidx) = [];
          fv(delidx) = [];
          mu(delidx) = [];
        end
        
        
        self.thres= mean(thres);
        
        self.fv = mean(cat(5,fv{:}),5);
        self.fv = self.fv/norm(self.fv(:));
        self.mu = mean(cat(5,mu{:}),5);
        
        
        if  numel(self.bgimg)>1
          % pixels likely to be above thres;
          notbgmsk = self.bgimg-mean(self.bgimg(:)) <= self.thres;
        else
          notbgmsk = 1;
        end
        self.msk = self.fgmsk & notbgmsk;

        % convert again according to the new LFD
        for i = 1:nTimes
          frame(:,:,i) = convertFrame(self,oframe(:,:,:,i),0);
        end

        self.mframe = 0;
        if nTimes>1
          self.mframe = self.mframe + mean(frame(:,:,2:end),3);
        end

      end
      verbose('Found (mean) threshold: %1.2f',self.thres);
      objh = round(mean(objh));
      objw = round(mean(objw));
      verbose('Mean object sizes: %d x %d pixels',objh,objw);


    end
    
    
    
    function bwimg = getBWImage(self,gpuFrame)

      frame = gather(gpuFrame);
      
      % bw image
      if self.inverse
        bwimg = (frame >= self.thres) & self.msk;
      else
        bwimg = (frame <= self.thres) & self.msk;
      end
      
      %save info
      %self.bwimg = bwimg;
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

      oframe = gather(oframe); % do on CPU;
      nTimes = size(oframe,4);
      savethres = self.thres;


      if ischar(savethres)  
         [thres,objh,objw,bw] = getAdaptiveAuto(self,oframe,nObjects);
      else
        self.thres = savethres;
        self.mu = [];
        self.fv = [];
        objh =  [];
        objw = [];
        bw = [];
      end
     
      
      if  numel(self.bgimg)>1
        % pixels likely to be above thres;
        notbgmsk = self.bgimg-mean(self.bgimg(:)) <= self.thres;
      else
        notbgmsk = 1;
      end
      self.msk = self.fgmsk & notbgmsk;
      
      % convert again according to the new LFD
      sz = size(oframe);
      frame = zeros([sz(1:2),nTimes]);
      for i = 1:nTimes
        frame(:,:,i) = convertFrame(self,oframe(:,:,:,i),0);
      end
      
      self.mframe = 0;
      if nTimes>1
        self.mframe = self.mframe + mean(frame(:,:,2:end),3);
      end

      if ~self.mtau
        verbose('Using auto frames as background!');
      end
      
      self.dframe = zeros(size(frame));      
      %self.frame = frame(:,:,1);
      self.framecounter = 0;
      
      if nargout
        varargout = {objh,objw};
      end

      if self.plotif
        % plot         
        clf;
        [r1,r2] = getsubplotnumber(nTimes);
        for i = 1:nTimes
          subplot(r1,r2,i);

          others = setdiff(1:nTimes,i);
          framei = frame(:,:,i)  - mean(frame(:,:,others),3) - self.thres;
          framei = (framei-min(framei(:)))/(max(framei(:))-min(framei(:)));

          if ~isempty(bw)
            L = labelmatrix(bw{i});
            rgbLabels = label2rgb(L);
            rgbFrame = repmat(255*max(min(framei,1),0),[1,1,3]);
            msk = repmat(mean(rgbLabels,3)==255,[1,1,3]);
            image(uint8(~msk.*single(rgbLabels) + msk.*rgbFrame));
          else
            imagesc(framei.*(framei>self.thres));
          end
          
          
        end
      end

      
    end
    
    
    
    
    function [mask,varargout] = step(self,oframe)
    %  [MASK,FRAME] = STEP(SELF,OFRAME) detects a foreground and returns a BW mask

      if isempty(self.mframe)
        init(self,oframe);
      end

      self.framecounter = self.framecounter +1;
      [frame, cframe] = convertFrame(self,oframe,1);
      mask = getBWImage(self,frame);

      %output
      if nargout>1
        varargout{1} = cframe; % RGB converted to intensity
      end
      if nargout>2
        varargout{2} = frame; % mean / dframe subtracted
      end
      
    end

    
  end
  
      

end