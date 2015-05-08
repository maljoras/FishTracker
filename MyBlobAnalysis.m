classdef MyBlobAnalysis < handle;
  
  properties 
    
    rprops = {'Centroid','MeanIntensity','BoundingBox','Orientation','Area'};
    minArea = 50;
    minextent = 50;
    maxextent = 1000;
    minmaxintensity = [0.1,0.5];
    nhistbins = 25;
    minMSERDistance = 3; % in pixels
    overlapthres = 0.8;
  
    fishwidth = 10;% approximate value in pixels
    fishlength = 50; 

    featurewidth = [];
    featureheight = [];
    colorfeature = true;
  end

  properties (SetAccess = private)
    mask = [];
    segm = [];
    frame = [];
    cframe = [];
  end
  
  methods (Access=private)
  
    function regionext = getMoreFeatures(self,regions,Iframe,Cframe)
      
      if isempty(regions)
        regionext = regions;
        return
      end
      
      regionext = [];
      mm = self.minmaxintensity;
      scale = diff(mm)/(self.nhistbins-1);
      cframe = reshape(Cframe,[],size(Cframe,3));

      for i = 1:length(regions)
        region = regions(i);
        

% $$$         msk = region.Image;
% $$$         fmsk = find(region.Image);
% $$$         dist = min(floor(bwdist(~msk)),self.maxdist);
% $$$         dist1 = dist(:);
% $$$         pv = max(min(region.PixelValues,mm(2)),mm(1))-mm(1);
% $$$         pv = floor(pv/scale)+1; % discretize
% $$$         hst = accumarray([dist1(fmsk),pv(:)],1,[self.maxdist,self.nhistbins],[],0);
% $$$         region.BWDist =  dist;        
% $$$         region.P2d =  hst/sum(hst(:)); % to prob

        % save some features 
        bb = floor(double(region.BoundingBox));
        bb(1:2) = max(bb(1:2),1);
        region.FilledImage =   Iframe(bb(2):bb(2)+bb(4)-1,bb(1):(bb(1)+bb(3)-1));
        bb2 = bb;
        bb2(1:2) = floor(max(bb(1:2)-bb(3:4)/2,1));
        bb2(3:4) = bb(3:4)*2;
        region.FilledImage2x =   Iframe(bb2(2):min(bb2(2)+bb2(4)-1,end),bb2(1):min((bb2(1)+bb2(3)-1),end));
        region.FilledImageCol = Cframe(bb(2):(bb(2)+bb(4)-1),bb(1):(bb(1)+bb(3)-1),:);
        region.FilledImageCol2x = Cframe(bb2(2):min(bb2(2)+bb2(4)-1,end),bb2(1):min((bb2(1)+bb2(3)-1),end),:);

        
        %% compute MSER features
        featim = region.FilledImage2x;
        range = [self.minArea,max(ceil(region.Area*0.75),self.minArea+10)];
        points= detectMSERFeatures(featim,'RegionAreaRange',range,'MaxAreaVariation',0.6);
        orgpoints = points;
        
        
        delinds = zeros(points.Count,1);
        for j = 1:length(points)
          if any(points(j).Location < bb2(3:4)/4) | any(points(j).Location > 3*bb2(3:4)/4)
            delinds(j) =  1;
          end
        end
        points = points(~delinds);

        if points.Count>1
          % delete very near points
          p = tril(squareform(pdist(points.Location)));
          x = 1;
          delinds = zeros(points.Count,1);
          while ~isempty(x)
            [x,y] = find(p<self.minMSERDistance & p>0,1);
            if ~isempty(x)
              if length(points(x).PixelList)>length(points(y).PixelList)
                z = y;
              else
                z = x;
              end
              delinds(z) = 1;
              p(:,y) = 0; % y still deleted
            end
          end
          points = points(~delinds);
        end


        region.MSERregions = points;
        region.MSERregions = [];
        for ii = 1:points.Count
          for f = {'Location','Orientation','Axes','PixelList'}
            region.MSERregions(ii,1).(f{1}) = points(ii).(f{1}); % re-map to structure
          end
        end

        region.MSERregionsOffset = bb2(1:2);
        
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
        extent = bbox(:,3) + bbox(:,4);
        delidx =  extent<self.minextent | extent>self.maxextent | self.minArea>area;
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
            idx = s2i(size(Iframe),bsxfun(@plus,double(points(ii).PixelList(:,[2,1])),rp(i).MSERregionsOffset([2,1])-1));
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
              if overlap(ii,jj)<self.overlapthres
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
          
          newrp = regionprops(newspots,Iframe,[self.rprops,{'Image','PixelValues','PixelIdxList'}]);
          newrp = self.getMoreFeatures(newrp,Iframe, Cframe);

% $$$           % use old MSER regions
% $$$           s = 0;
% $$$           for j = find(~delsp)
% $$$             s = s+1;
% $$$             newrp(s).MSERregions =  rp(i).MSERregions(j);
% $$$           end
% $$$                     
          newrplist = cat(1,newrplist, newrp);
          
        end
        
        
      end


      rp = [rp(~splitif);newrplist];
      
      
    end
    
    
    function segments = getFishFeatures(self,segments);
      
      plotif = 0; % for debugging
      interpif = 0; % much better effect it seems
      
      if plotif 
        figure(1);
      end

      fixedwidth = self.fishwidth*2;
      fixedheight = self.fishlength;

      smallerwidth = self.featurewidth;
      smallerheight = self.featureheight;

      debendsmoothing = 5; % needs to be odd
      avgheadonset = 5;

      for i = 1:length(segments)
        seg = segments(i);

        segments(i).fishFeature = [];
        segments(i).bendingStdValue = [];


        if ~isempty(seg.MSERregions)
          oimg = seg.FilledImage2x;
          if self.colorfeature
            oimg_col = double(seg.FilledImageCol2x);
          end
          
          ori = -seg.MSERregions(1).Orientation/2/pi*360+ 90;
          lst = seg.MSERregions(1).PixelList;
          msk = zeros(size(oimg));
          msk(s2i(size(oimg),lst(:,[2,1]))) = 1;
        else
          % take normal regions
          %oimg = seg.FilledImage;
          %msk =  seg.Image;

          oimg = seg.FilledImage2x;
          if self.colorfeature
            oimg_col = double(seg.FilledImageCol2x);
          end
          
          msk = oimg<mean(oimg(:)); % very crude
          ori = -seg.Orientation + 90;
        end
        
                
        if plotif
          subplot(2,length(segments),i)
          imagescbw(oimg);
          daspect([1,1,1]);
          axis off;
        end
        

        % regenerate mask to get mean background
        mback = mean(oimg(~msk));
        if self.colorfeature
          tmp = mean(oimg_col,3);
          mback_col = mean(tmp(~msk));
        end
        
        
        % rotate fish to axes
        oimg = 1-max(imrotate(1-oimg,ori,'bilinear','loose'),1-mback);
        msk = imrotate(msk,ori,'nearest','loose');
        
        if self.colorfeature
          oimg_col= imrotate(oimg_col,ori,'bilinear','loose');
          oimg_col(~oimg_col) = mback_col;
        end
        
        
        % get updated bounding box of the fish
        rp1 = regionprops(msk);
        [~,midx] = max(cat(1,rp1.Area));
        center = rp1(midx).Centroid;
        bb = max(floor(rp1(midx).BoundingBox),1);
        
        
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
          %warning('Fishheight probably too small')
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
        [indy,indx] = ndgrid(single(1:fixedheight),single(1:fixedwidth));
        
        if self.colorfeature        
          doimg_col = zeros(size(foimg_col));
        end
        
        if ~interpif
          indx = mod(bsxfun(@plus,indx,round(com-fixedwidth/2))-1,fixedwidth)+1;
          idx = s2i(size(foimg),[indy(:),indx(:)]);
          doimg = reshape(foimg(idx),size(foimg));
          if self.colorfeature
            for ii = 1:size(doimg_col,3)
              tmp = foimg_col(:,:,ii);
              doimg_col(:,:,ii) =  reshape(tmp(idx),size(tmp));
            end
          end
          
        else
          indxI = max(min(bsxfun(@plus,indx,com-fixedwidth/2),fixedwidth),1);
          doimg = interp2(indx,indy,foimg,indxI,indy,'linear',mback);
          
          if self.colorfeature
            for ii = 1:size(doimg_col,3)
              tmp = foimg_col(:,:,ii);
              doimg_col(:,:,ii) = interp2(indx,indy,tmp,indxI,indy,'linear',mback_col);
            end
          end
        end


        
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
          segments(i).fishFeature = final_oimg_col;
        else
          segments(i).fishFeature = final_oimg;
        end
        
        %final_oimg = (final_oimg(:,1:smallerwidth/2) + final_oimg(:,end:-1:smallerwidth/2+1))/2;
        %final_oimg = (final_oimg-mean(final_oimg(:)))/std(final_oimg(:));
        %x = dct2(x);

        segments(i).bendingStdValue = std(com(1:min(smallerheight,end)));
        
        if plotif
          figure(1)
          subplot(2,length(segments)*2,2*(i-1)+length(segments)*2+1);
          if self.colorfeature
            imagesc(foimg_col);
          else
            imagesc(foimg);
          end
          
          daspect([1,1,1]);
          axis off;
          hold on ;
          plot(com,1:fixedheight,'r','linewidth',2)


          subplot(2,length(segments)*2,2*(i-1)+length(segments)*2+2);
          if self.colorfeature
            imagesc(final_oimg_col);
          else
            imagesc(final_oimg);
          end
          daspect([1,1,1]);
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
      spots=bwconncomp(bwimg);
      rp = regionprops(spots,Iframe,[self.rprops,{'Image','PixelValues','PixelIdxList'}]);

      goodmsk = self.getGoodMsk(rp);

      % compute more features
      rp = self.getMoreFeatures(rp(goodmsk),Iframe, Cframe);
      rp = self.splitRegions(rp,Iframe, Cframe);
      rp = self.getFishFeatures(rp);

      segm = rp;
    end
    
  end


  methods
    
    function self = MyBlobAnalysis(varargin) % constructor
      
      self = self@handle();

      
      self.featurewidth = self.fishwidth;
      self.featureheight = floor(self.fishlength*0.6);

      
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
      


      verbose('Initiated BlobAnalyzer with fish features.')
      verbose('Assumed approx. fish dimensions in pixels: %d x %d',self.fishlength, self.fishwidth)
    end
    
    
      
      
    function segm = step(self,mask, Iframe, Cframe)
      
      if nargin<4
        Cframe = Iframe;
      end
      
      % detect objects
      segm = detect(self,mask, Iframe, Cframe);

      % save
      self.segm = segm;
      self.mask = mask;
      self.frame = Iframe;
      self.cframe = Cframe;
     end

    
    function plot(self)

      if isempty(self.segm)
        error('Call "step" first');
      end
      
      assert(any(strcmp(self.rprops,'BoundingBox')));
      assert(any(strcmp(self.rprops,'Centroid')));
      
      cla;
      imagesc(self.frame);
      hold on;
      plot(self.segm.Centroid(:,1),self.segm.Centroid(:,2),'xr','markersize',10,'linewidth',2);
      
      if ~isempty(self.segm.BoundingBox)
        z = bbox2points(self.segm.BoundingBox);
        plot(squeeze(z(:,1,:)),squeeze(z(:,2,:)),'r');
      end
      
      drawnow;
    end

    
  end

end

      

    