classdef FishBlobAnalysisMatlab < FishBlobAnalysis;
  
  
  methods

    function  [oimg,msk,oimg_col] = a_imrotate(self,ori,oimg,mback,omsk,omsk_back,oimg_col,mback_col)
    % rotate fish to axes
      oimg = 1-max(imrotate(1-oimg,ori,'bilinear','loose'),1-mback);
      msk = imrotate(omsk,ori,'nearest','loose');
      
      if ~isempty(oimg_col)
        oimg_col= 1-max(imrotate(1-oimg_col,ori,'bilinear','loose'),1-mback_col);
        oimg_col(~oimg_col) = mback_col;
      end
    end
    
    function rp = a_getRegions(self,bwimg,Iframe,rprops);
      bwimg = bwareaopen(bwimg,self.minArea); % needed ?
      spots=bwconncomp(bwimg,8);
      rp = regionprops(spots,Iframe,[self.rprops,{'Image'}]);
    end
    
    function [bb,center] = a_getMaxAreaRegion(self,bwimg);
      rp = regionprops(bwimg);
      [~,midx] = max(cat(1,rp.Area));
      bb = max(floor(rp(midx).BoundingBox),1);
      center = rp(midx).Centroid;
    end
    
    function doimg = a_interp(self,foimg,xshift,mback,type)

      [h,w] = size(foimg);
      indw = single(1:w);
      indh = single(1:h);
      indxI = max(min(bsxfun(@plus,indw,xshift),w),1);
      [indy,indx] = ndgrid(indh,indw);

      if size(foimg,3)==3 % colorimage
        doimg = zeros(size(foimg));
        for ii = 1:size(doimg,3)
          tmp = foimg(:,:,ii);
          doimg(:,:,ii) = interp2(indx,indy,tmp,indxI,indy,type,mback);
        end
      else
        doimg2 = interp2(indx,indy,foimg,indxI,indy,type,mback);
      end
    end

    
    function [region] = a_computeMSERregions(self,region,bb2)
    % not yet implemented. USE MATLAB VERSION FOR NOW

      featim = region.FilledImage2x;
      p= detectMSERFeatures(featim,'RegionAreaRange',range,'MaxAreaVariation',0.6);
      loc = p.Location;
      pixelList = p.PixelList;
      ori = p.Orientation/2/pi*360;


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
              d(:,z) = 0; % y still deleted
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
      end
      
    end

    
    function [region] = a_computeMSERregions(self,region,bb2);

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
      

      %region.MSERregions = points;

      for ii = 1:points.Count
        for f = {'Location','Orientation','Axes','PixelList'}
          region.MSERregions(ii,1).(f{1}) = points(ii).(f{1}); % re-map to structure
        end
      end
      
      region.MSERregionsOffset = bb2(1:2);

    end
    
    function self = FishBlobAnalysisMatlab(varargin) % constructor
      
      self = self@FishBlobAnalysis(varargin{:});
    end
  
  end

end

      

    