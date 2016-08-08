classdef FishBlobAnalysis < handle;
  
  properties 
    
    rprops = {'Centroid','BoundingBox','Orientation','Area'};

    computeMSER = 0;
    minMSERDistance = 3; % in pixels
    computeMSERthres = 3; % when to compute MSER (extent larger than computerMSERthres*fishwidth*fishlength)
    overlapthres = 0.8;

    fishwidth = 20;% approximate value in pixels
    fishlength = 100; 
    headprop = 0.6; % proportion of fishlength to use as feature;
    
    interpif = 1; % much better effect it seems
    readjustposition = false;
    plotif = 0;    
    se = [];
    segm = [];
    colorfeature = false;


  end


  properties (Access=protected)

    mser = [];
    minArea  = [];
    maxArea = [];
    minextent = [];
    maxextent = [];
    minWidth  = [];
    featurewidth = [];
    featureheight = [];
    
  end
  
  
  methods (Access=protected) %%% THIS FUNCTIONS CAN BE OVERLOADED
    
    function a_init(self);
    % initialize mser
      if isempty(self.mser)
        self.mser  =  cv.MSER('MinArea',self.minArea,'MaxArea',self.maxArea,'MaxVariation',0.6);
        if isa(self.se,'strel')
          self.se = getnhood(self.se);
        end
      end
      
    end
    

    
    function  [oimg,msk,oimg_col] = a_imrotate(self,ori,oimg,mback,msk,oimg_col,mback_col)
    % rotate fish to axes
      
      sz = size(oimg);
      p = ori/180*pi;
      a = cos(p);
      b = sin(p);

      szin = sz([2,1]);
      szout = [szin(1)*abs(a)+szin(2)*abs(b),szin(1)*abs(b)+szin(2)*abs(a)];

      %T = cv.getRotationMatrix2D(sz/2,ori,1.);
      x = szin(1)/2;
      y = szin(2)/2;
      offs = (szout-szin)/2;
      T = [ a, b, (1-a)*x-b*y + offs(1); -b, a, b*x+(1-a)*y + offs(2) ];

      oimg = cv.warpAffine(oimg,T,'DSize',szout,'Interpolation','Linear','BorderValue',mback);
      msk = cv.warpAffine(msk,T,'DSize',szout,'Interpolation','Nearest','BorderValue',0);
      
      if ~isempty(oimg_col)
        oimg_col = cv.warpAffine(oimg_col,T,'DSize',szout,'Interpolation','Linear','BorderValue',mback_col);
      end
    end
  

    function [region] = a_computeMSERregions(self,region,bb2)

      featim = region.FilledImage2x;

      if ~isa(featim,'uint8')
        featim = uint8(featim*255);
      end
      featim(~region.Image2x) = uint8(mean(featim(~region.Image2x)));

      
      [p,bboxes] = self.mser.detectRegions(featim);
      loc = zeros(length(p),2);
      ori = zeros(length(p),1);
      pixelList = cell(length(p),1);

      for i = 1:length(p)
        info = cv.fitEllipse(p{i});
        ori(i) = info.angle;
        pixelList{i} = cat(1,p{i}{:});
        loc(i,:) = info.center;
      end
      idxmsk = [];
      if size(loc,1)>1
        % select some of the regions
        idxmsk =  all(bsxfun(@ge,loc,bb2(3:4)/4) & bsxfun(@le,loc,3*bb2(3:4)/4),2);
      
        if any(idxmsk)
          % delete very near points
          d = tril(squareform(pdist(loc)));
          x = 1;
          while ~isempty(x)
            [x,y] = find(d<self.minMSERDistance & d>0,1);
            if ~isempty(x)
              if length(pixelList{x})>length(pixelList{y})
                z = y;
              else
                z = x;
              end
              idxmsk(z) = 0;
              d(:,y) = 0; % y still deleted
            end
          end
        end
        
        
        %region.MSERregions = points;
        s = 0;
        for ii = find(idxmsk)'
          s = s+1;
          region.MSERregions(s).Location = loc(ii,:);
          region.MSERregions(s).Orientation = ori(ii);
          region.MSERregions(s).PixelList = pixelList{ii};
        end
        
        region.MSERregionsOffset = bb2(1:2);
      
      
            
      
        DEBUG = 0;
        if DEBUG
          
          p1= detectMSERFeatures(featim,'MaxAreaVariation',0.6);
          
          subplot(1,2,1);
          cla;
          imshow(featim);hold on;
          plot(p1, 'showPixelList', true, 'showEllipses', false);      
          title('matlab')
          subplot(1,2,2);
          cla;
          imshow(featim);hold on;
          idx =find(idxmsk);
          if ~isempty(idx)
            pixelList =cellfun(@int32,pixelList,'uni',0);
            pcv = MSERRegions(pixelList(idx));
            plot(pcv, 'showPixelList', true, 'showEllipses', true);
          end
          title('openCV')
          drawnow;
          keyboard
        end

      
      end
      
      
      
    end

    
    
    function [bb,center] = a_getMaxAreaRegion(self,bwimg);
    % only area,centroid,boundarybos from the max
      
      contours = cv.findContours(bwimg,'Mode','External','Method','Simple');
      
      if length(contours)>1
        area = zeros(length(contours),1);
        % only from largest needed
        for i = 1:length(contours)
          area(i) =  cv.contourArea(contours{i});;
        end
        [~,imxarea] = max(area);
        c =contours{imxarea};
      else
        c =contours{1};
      end
      
      bb =  cv.boundingRect(c);
      if length(c)>4
        info = cv.fitEllipse(c);      
        center = info.center([1,2]);
      else
        center = bb([1,2])+ bb([3,4])/2;
      end
      
    end
    

    function doimg = a_interp(self,foimg,xshift,mback,type)

      [h,w] = size(foimg);
      indw = single(0:w-1);
      indh = single(0:h-1)';
      indxI = max(min(bsxfun(@plus,indw,xshift),w),1);
      indy = repmat(indh,[1,w]);

      % color should work, too.
      doimg = cv.remap(foimg,indxI,indy,'BorderType','Constant','BorderValue',mback,'Interpolation',type);

    end

    function msk = a_closeHoles(self,msk);
    % not needed. Contours are filled.
    %msk = cv.dilate(msk,'Element',self.se,'BorderValue',0,'Iterations',1);
    %  msk = cv.erode(msk,'BorderValue',0,'Iterations',2);
    %  msk = cv.dilate(msk,'Element',self.se,'BorderValue',0,'Iterations',1);
    %  msk = logical(msk);
    end
    
    
    function rp = a_getRegions(self,bwimg,Iframe,rprops);
      
      contours = cv.findContours(bwimg,'Mode','External','Method','Simple');
% $$$       isuint = isa(Iframe,'uint8');

      rp = [];
      sizei = size(bwimg,1);
      s = 0;
      for i = 1:length(contours)
        
        area = cv.contourArea(contours{i});;
        if area< self.minArea
          continue;
        end
        s = s+1;
        % BoundingBox
        bb = cv.boundingRect(contours{i});
        if length(contours{i})>4
          info = cv.fitEllipse(contours{i});
        else
          info.center = bb([1,2]) + bb([3,4])/2;
          info.size = bb([3,4]);
          info.angle = 0;
        end
        
        rp(s,1).Area = area;

        rp(s,1).Centroid = info.center;
        rp(s,1).BoundingBox = bb;
        rp(s,1).Size = info.size;
        rp(s,1).MajorAxisLength = max(info.size);
        rp(s,1).MinorAxisLength = min(info.size);
        rp(s,1).Orientation = -info.angle + 90;
        %smallmask = bwimg(bb(2)+1:bb(2)+bb(4),bb(1)+1:(bb(1)+bb(3)));
        msk = zeros(bb([4,3]));
        rp(s,1).Image =  cv.drawContours(msk,contours(i),'Offset',-bb([1,2]),'Thickness',-1,'MaxLevel',2);

      end
    end
    
  end

  
  
  
  methods % GENERAL ACCESS METHODS

    function regionext = getMoreFeatures(self,regions,Iframe,Cframe)
      
      if isempty(regions)
        regionext = regions;
        return
      end
      
      regionext = [];

      if ~isfield(regions,'Size');
        [regions.Size] = deal([]);
      end
      if ~isfield(regions,'reversed');
        [regions.reversed] = deal(0);
      end
      
      
      for i = 1:length(regions)
        region = regions(i);
        bb = floor(double(region.BoundingBox));
        
        if isempty(region.Size)
          region.Size = [region.MajorAxisLength,region.MinorAxisLength];
        end
        
        % save some features 
        
        bb(1:2) = max(bb(1:2),0);
        bb2 = bb;
        bb2(1:2) = max(floor(bb(1:2)-bb(3:4)/2),0);
        bb2(3:4) = bb(3:4)*2;
        % border effect
        bb2(3) = min(bb2(3) + bb2(1),size(Iframe,2))-bb2(1);
        bb2(4) = min(bb2(4) + bb2(2),size(Iframe,1))-bb2(2);
        region.BoundingBox2x = bb2;
        
        region.FilledImage =   Iframe(bb(2)+1:bb(2)+bb(4),bb(1)+1:(bb(1)+bb(3)));
        region.FilledImage2x =   Iframe(bb2(2)+1:(bb2(2)+bb2(4)),bb2(1)+1:(bb2(1)+bb2(3)));
        
        if ~isempty(Cframe)
          region.FilledImageCol = Cframe(bb(2)+1:(bb(2)+bb(4)),bb(1)+1:(bb(1)+bb(3)),:);
          region.FilledImageCol2x = Cframe(bb2(2)+1:(bb2(2)+bb2(4)),bb2(1)+1:(bb2(1)+bb2(3)),:);
          
        else % needed?
          region.FilledImageCol = region.FilledImage;
          region.FilledImageCol2x = region.FilledImage2x;
        end
        
        
        % enlarge mask
        region.Image2x = logical(zeros(bb2([4,3])));
        offset = bb(1:2)-bb2(1:2);
        region.Image2x(offset(2)+1:offset(2)+bb(4),offset(1)+1:offset(1)+bb(3)) = region.Image;
        
        
        % NOTE: close the holes of the msk
        %region.Image = self.a_closeHoles(region.Image);
        region.Image2x= self.a_closeHoles(region.Image2x); % better borders
        
        
        %% compute MSER features
        region.MSERregions = []; % field is expected by some functions...
        region.MSERregionsOffset = [];

        if self.computeMSER &&  bb(3)*bb(4)>self.computeMSERthres*self.fishwidth*self.fishlength
          xy.helper.verbose('%d>%f1.0\r',bb(3)*bb(4),self.computeMSERthres*self.fishwidth*self.fishlength)
          region = self.a_computeMSERregions(region,bb2);
        end
        
        %% add 
        regionext = cat(1,regionext,region);
      end
    end
    
    function goodmsk = getGoodMsk(self,rp);

      nspots = length(rp);
      goodmsk = true(nspots,1);
      
      % reduce according to size (width + height in pixels)
      if nspots
        bbox =  cat(1,rp.BoundingBox);
        area = cat(1,rp.Area);
        extent = cat(1,rp.MajorAxisLength) + cat(1,rp.MinorAxisLength);
        width = min(bbox(:,[3,4]),[],2);
        delidx =  extent<self.minextent | extent>self.maxextent | self.minArea>area | self.minWidth>width;
        goodmsk(delidx) = false;
      end

    end
    
    function rp = splitRegions(self,rp,Iframe, Cframe)
    % split regions based on the MSER detection
      
      nspots = length(rp); 

      %% split regions
      splitif = zeros(nspots,1);
      newrplist = [];
      for i = 1:nspots

        if length(rp(i).MSERregions)>1
          splitif(i) = 1;

          featim = rp(i).FilledImage2x;
          points = rp(i).MSERregions;
          
          newspots = [];
          newspots.Connectivity = 4; % dummy
          newspots.ImageSize = size(Iframe);
          newspots.PixelIdxList = {};
          for ii = 1:length(points)
            idx = xy.helper.s2i(size(Iframe),bsxfun(@plus,double(points(ii).PixelList(:,[2,1])),rp(i).MSERregionsOffset([2,1])-1));
            newspots.PixelIdxList{1,ii} = idx;
          end
          newspots.NumObjects = length(newspots.PixelIdxList);

          % check for overlapping
          overlap = zeros(newspots.NumObjects,newspots.NumObjects);
          for ii = 1:newspots.NumObjects
            pl1 = newspots.PixelIdxList{ii};
            for jj = ii+1:newspots.NumObjects
              pl2 = newspots.PixelIdxList{jj};
              idx = intersect(pl1,pl2);
              
              overlap(ii,jj) = length(idx)/min(length(pl2),length(pl1));
              if overlap(ii,jj)>self.overlapthres
                newspots.PixelIdxList{jj} = setdiff(pl2,idx);
                newspots.PixelIdxList{ii} = setdiff(pl1,idx);
              else
                % no deleting
              end
              
            end
            
          end

          % do not handle regions that have massive overlap. Take
          % only smaller one in this case
          [ii,jj] = find(overlap>self.overlapthres);
          delsp = zeros( newspots.NumObjects,1);
          for k = 1:length(ii)
            if length(newspots.PixelIdxList{ii(k)})>length(newspots.PixelIdxList{jj(k)})
              delsp(ii(k)) = 1;
            else
              delsp(jj(k)) = 1;
            end
          end
          newspots.PixelIdxList(delsp==1) =[];
          newspots.NumObjects = length(newspots.PixelIdxList);

          newrp = regionprops(newspots,Iframe,[self.rprops,{'Image','MajorAxisLength','MinorAxisLength'}]);
          newrp(~cat(1,newrp.Area)) = [];
          if ~isempty(newrp)
            newrp = self.getMoreFeatures(newrp,Iframe, Cframe);
            newrplist = cat(1,newrplist, newrp);
          end
          
          
        end
        
      end


      try % the strucutres might have different fields so first
          % have a try
        rp = [rp(~splitif);newrplist];
      catch 
        f = fieldnames(rp)';
        rp = rp(~splitif);
        for i = 1:length(newrplist)
          rpnew = cell2struct(cell(length(f),1),f);
          for ff = fieldnames(newrplist)'
            rpnew.(ff{1}) = newrplist(i).(ff{1});
          end
          rp = [rp;rpnew];
        end
      end
      
     
    end
    
    
    
    
    function segments = getFishFeatures(self,segments);
      
      plotif = self.plotif; % for debugging

      
      if plotif 
        figure(1);
      end

      fixedwidth = self.fishwidth*3;
      fixedheight = self.fishlength;

      smallerwidth = self.featurewidth;
      smallerheight = self.featureheight;

      debendsmoothing = 5; % needs to be odd
      avgheadonset = 2;%5

      for i = 1:length(segments)
        seg = segments(i);

        segments(i).FishFeature = [];
        segments(i).bendingStdValue = [];

        
    % $$$         if ~isempty(seg.MSERregions)
    % $$$           oimg = seg.FilledImage2x;
    % $$$           if self.colorfeature && isfield(seg,'FilledImageCol2x')
    % $$$             oimg_col = seg.FilledImageCol2x;
    % $$$           end
    % $$$           
    % $$$           ori = -seg.MSERregions(1).Orientation+ 90; %! why take
    % $$$                                                      %the first?
    % $$$           lst = seg.MSERregions(1).PixelList;
    % $$$           if ~lst(end)
    % $$$             % somehow zeros at the end ? WHY ?
    % $$$             lst(~any(lst,2),:) =  [];
    % $$$           end
    % $$$           omsk = zeros(size(oimg));
    % $$$           omsk(xy.helper.s2i(size(oimg),lst(:,[2,1]))) = 1;
    % $$$         else
    % $$$           % take normal regions
        
        oimg = seg.FilledImage2x;
        omsk =  seg.Image2x; % Note: Image has no closed borders
                             % (only Image2x)
        if self.colorfeature && isfield(seg,'FilledImageCol2x')
          oimg_col = seg.FilledImageCol2x;
        end
        

        ori = -seg.Orientation + 90;
    % $$$         end

        % need single for the features
        if isa(oimg,'uint8')
          oimg = single(oimg)/255;
        end
        
        if self.colorfeature && isa(oimg_col,'uint8')
          oimg_col = single(oimg_col)/255;
        end
        
        
        
        if plotif
          subplot(2,length(segments),i)
          if self.colorfeature
            tmp = min(max(oimg_col*255,0),255);
          else
            tmp = min(max(round(255*((oimg-min(oimg(:)))/diff(xy.helper.minmax(oimg(:)')))),0),255);
            tmp = repmat(tmp,[1,1,3]);
          end
          % delete green channel within mask
          tmp2 = tmp(:,:,2);
          tmp2(omsk) = 50;
          tmp(:,:,2) = tmp2;
          image(uint8(tmp));
          
          daspect([1,1,1]);
          axis off;
        end
        

        % regenerate mask to get mean background
        mback = mean(oimg(~omsk));
        if self.colorfeature
          tmp = mean(oimg_col,3);
          mback_col = mean(tmp(~omsk));
        end
        
        
        % rotate fish to axes
        if self.colorfeature
          [oimg,msk,oimg_col] = self.a_imrotate(ori,oimg,mback,omsk,oimg_col,mback_col);
        else
          [oimg,msk] = self.a_imrotate(ori,oimg,mback,omsk,[],[]);
        end

        
        % get updated bounding box of the fish
        [bb,center] = self.a_getMaxAreaRegion(msk);


        % get direction of movement (assuming the fish get's "thinner")
        y =sum(msk,2);
        inds = find(y);
        X = [ones(size(y(inds))),(1:length(inds))'];
        w = X\y(inds);
        if w(2)>0 % probably backwards
          msk = flipud(msk);
          oimg = flipud(oimg);
          center(2) = size(oimg,1)-center(2) + 1;
          if self.colorfeature
            for ii = 1:size(oimg_col,3)
              oimg_col(:,:,ii) = flipud(oimg_col(:,:,ii));
            end
          end
          
        end
        
        % maybe use BB to check whether fixed width is good choice
        if  fixedheight/2 < bb(4)/2
          %xy.helper.verbose('WARNING: Fishheight probably too small')
        end
        
        startidx = max(find(any(msk,2),1,'first')-2,1);
        oimg = oimg(startidx:end,:);
        oimg(max(end-startidx+2,1):end,:) = mback;
        center(2) = center(2)-startidx+1;

        if self.colorfeature
          oimg_col = oimg_col(startidx:end,:,:);
          oimg_col(max(end-startidx+2,1):end,:,:) = mback_col;
        end

        
        % select fixed window
        fbb = [center,fixedwidth,fixedheight];
        fbb([1,2]) = max(round(fbb([1,2])-fbb([3,4])/2),1);

        foimg = mback*ones(fixedheight,fixedwidth);;
        indy = fbb(2):min(fbb(2)+fbb(4)-1,size(oimg,1));
        indx = fbb(1):min(fbb(1)+fbb(3)-1,size(oimg,2));
        foimg(1:length(indy),1:length(indx)) = oimg(indy,indx);

        if self.colorfeature
          foimg_col = mback_col*ones(fixedheight,fixedwidth,3);
          foimg_col(1:length(indy),1:length(indx),:) = oimg_col(indy,indx,:);
        end
        
        % compute center of mass to "de-bend" the fish
        x= 1:fixedwidth;
        bw = round(debendsmoothing); % some smooting
        z = 1-foimg; % assume black fish on white !!
        z = max(bsxfun(@minus,z,mean(z,2)),0);
        com = z*x'./(sum(z,2));
        com(isnan(com)) = fixedwidth/2;
        com = conv(com,ones(bw,1)/bw,'valid');
        com = [com(1:floor(bw/2));com;com(end-floor(bw/2)+1:end)];

        
        
        
        if ~self.interpif
          [indy,indx] = ndgrid(single(1:fixedheight),single(1:fixedwidth));
          if self.colorfeature        
            doimg_col = zeros(size(foimg_col));
          end

          indx = mod(bsxfun(@plus,indx,round(com-fixedwidth/2))-1,fixedwidth)+1;
          idx = xy.helper.s2i(size(foimg),[indy(:),indx(:)]);
          doimg = reshape(foimg(idx),size(foimg));
          if self.colorfeature
            for ii = 1:size(doimg_col,3)
              tmp = foimg_col(:,:,ii);
              doimg_col(:,:,ii) =  reshape(tmp(idx),size(tmp));
            end
          end
          
        else
          indxI = max(min(bsxfun(@plus,single(1:fixedwidth),com-fixedwidth/2),fixedwidth),1);
          indy = repmat(single(1:fixedheight)',[1,size(indxI,2)]);
          doimg = self.a_interp(foimg,com-fixedwidth/2,mback,'Linear');
          if self.colorfeature
            doimg_col = self.a_interp(foimg_col,com-fixedwidth/2,mback_col,'Linear');
          end
        end



        if self.readjustposition
          % readjust position to start 
          mline = mean(doimg(1:floor(fixedheight/2),:),2);
          %headidx2 = find(cumsum(abs(diff(mline))<std(diff(mline(avgheadonset+1:end))))>3,1)-1;
          bw1 = 3;
          [~,headidx] = max(diff(conv(mline,ones(bw1,1)/bw1,'valid'),2));
          headidx = headidx + 1 + (bw1-1)/2;
          
          startidx = headidx - avgheadonset;
          if startidx>1
            doimg = doimg(startidx:end,:);
            doimg(end-startidx+2:end,:) = mback;
            
            if self.colorfeature
              doimg_col = doimg_col(startidx:end,:,:);
              doimg_col(end-startidx+2:end,:,:) = mback_col;
            end
            
          end
        end
        
        

        if plotif>2
          figure(3);
          clf;
          plot(mline,'.')
          hold on;
          plot(headidx,mline(headidx),'or');
          plot(headidx2,mline(headidx2),'xb');
          drawnow;
        end
        
        % take middle part (the fish)
        istart = floor(fixedwidth/2 - smallerwidth/2);
        jstart = 1;
        %final_oimg = zeros(smallerheight,smallerwidth);
        idx1 = jstart:jstart+smallerheight-1;
        idx2 = istart:istart+smallerwidth-1;
        final_oimg = doimg(min(idx1,end),min(idx2,end));

        if self.colorfeature
          final_oimg_col = doimg_col(min(idx1,end),min(idx2,end),:);
          segments(i).FishFeature = final_oimg_col;
        else
          segments(i).FishFeature = final_oimg;
        end
        
        %final_oimg = (final_oimg(:,1:smallerwidth/2) + final_oimg(:,end:-1:smallerwidth/2+1))/2;
        %final_oimg = (final_oimg-mean(final_oimg(:)))/std(final_oimg(:));
        %x = dct2(x);

        segments(i).bendingStdValue = std(com(1:min(smallerheight,end)));
        
        if plotif
          figure(1)

          subplot(2,length(segments)*2,2*(i-1)+length(segments)*2+1);
          
          if isfield(seg,'FishFeatureC')
            imagesc(seg.FishFeatureCRemap);
            title('C')
          else
            
            if self.colorfeature
              imagesc(foimg_col);
            else
              imagesc(foimg);
            end
            %daspect([1,1,1]);
            hold on ;
            plot(com,1:fixedheight,'r','linewidth',2)
          
          end
          axis off;

          subplot(2,length(segments)*2,2*(i-1)+length(segments)*2+2);
          if self.colorfeature
            imagesc(final_oimg_col);
          else
            imagesc(final_oimg);
          end
          %daspect([1,1,1]);
          axis off;
          
        end
      end
      

      if plotif
        drawnow;
        pause(0.1);
      end
      
    end
    
    
    function segm = detect(self, bwimg, Iframe, Cframe)
      
    % get the spots from the binary image
      rp = self.a_getRegions(bwimg,Iframe,[self.rprops,{'Image'}]);
      
      goodmsk = self.getGoodMsk(rp);
      
      %fprintf('%d/%d\r',sum(goodmsk),length(goodmsk));
      % compute more features
      rp = self.getMoreFeatures(rp(goodmsk),Iframe, Cframe);
      rp = self.splitRegions(rp,Iframe, Cframe);
      %goodmsk = self.getGoodMsk(rp);
      %rp = rp(goodmsk);
      rp = self.getFishFeatures(rp);

      
      segm = rp;
    end
    
  end


  methods
    
    function self = FishBlobAnalysis(varargin) % constructor
      
      self = self@handle();
     
      nargs = length(varargin);
      if nargs==1 && isstruct(varargin{1})

        for f = fieldnames(varargin{1})'
          if isprop(self,f{1})
            try
              self.(f{1}) = varargin{1}.(f{1});
            end
            
          end
        end
        if isfield(varargin{1},'blob')
          for f = fieldnames(varargin{1}.blob)'
            if isprop(self,f{1})
              try
                self.(f{1}) = varargin{1}.blob.(f{1});
              end
            end
          end
        end
      elseif nargs>0 && mod(nargs,2)
        error('expected arguments of the type ("PropName",pvalue)');
      else
        for i = 1:2:nargs
          if ~isprop(self,varargin{i})
            error('expected arguments of the type ("PropName",pvalue)');
          else
            self.(varargin{i}) = varargin{i+1};
          end
        end
      end
      
      
      if isempty(self.se)
        self.se = strel('disk',1);
      end

      initialize(self);
    end
    
    function featureSize = getFeatureSize(self);
      featureSize  = [self.featureheight,self.featurewidth];
    end
    
    function initialize(self,verboseif)
      if ~exist('verboseif','var')
        verboseif = 1;
      end
      if isempty(self.fishwidth)
        self.fishwidth = 20; % random guess. Can be set afterwords;
      end
      if isempty(self.fishlength)
        self.fishlength = 100; % random guess. Can be set afterwords;
      end
      if isempty(self.headprop)
        self.headprop = 0.6; % random guess. Can be set afterwords;
      end
      
      setFishSize(self,self.fishlength,self.fishwidth,self.headprop,verboseif);
    end
    
    
    function setFishSize(self,fishlength,fishwidth,headprop,verboseif)
    % adjust paremeters according to the fishsize
      
      if ~exist('verboseif','var')
        verboseif = 1;
      end
      if ~exist('headprop','var')
        headprop = self.headprob;
      end
      self.headprop = headprop;

      if ~isempty(fishlength)
        self.fishlength = fishlength;
      end
      if ~isempty(fishwidth)
        self.fishwidth = fishwidth;        
      end

      if ~isempty(self.fishlength) && ~isempty(self.fishwidth) 
        if verboseif
          xy.helper.verbose('Initiated %s with fish features.',upper(class(self)))
          xy.helper.verbose('Assumed approx. fish dimensions in pixels: %d x %d',self.fishlength, self.fishwidth)
        end
      else
        error('Provide valid fishsizes');
      end


      self.featurewidth = self.fishwidth;
      self.featureheight = floor(self.fishlength*headprop);

      self.maxextent = 5*(self.fishlength+self.fishwidth);
      self.minextent = 0.05*(self.fishlength+self.fishwidth);
      self.minWidth = 0.05*self.fishwidth;
      self.minArea  = 0.05*self.fishlength*self.fishwidth;
      self.maxArea  = 5*self.fishlength*self.fishwidth;
      self.a_init();
    end
      
    function segm = stepBlob(self,bwimg, Iframe, Cframe)
    % SEGM = STEP(SELF,BWIMG, IFRAME, CFRAME).  IFRAME is an intensity
    % based frame (i.e. one channel) and CFRAME can be color frame
    % and is optional. 
      
      if nargin<4
        Cframe = [];
      end
      
      % detect objects
      segm = detect(self,bwimg, Iframe, Cframe);

      % save for plotting
      self.segm = segm;
     end

    
  
    
  end

end

      

    