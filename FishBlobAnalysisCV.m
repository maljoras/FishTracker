classdef FishBlobAnalysisCV < FishBlobAnalysis;
  
  
  properties 
    mser = [];
  end
  

  
  
  methods
    
    function a_init(self);
      % initialize mser
      self.mser = [];
      area = 1.5*self.fishwidth*self.fishlength;
      range =  [self.minArea,max(area,self.minArea+10)];
      self.mser  =  cv.MSER('MinArea',range(1),'MaxArea',range(2),'MaxVariation',0.6);
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
    % not yet implemented. USE MATLAB VERSION FOR NOW

      featim = region.FilledImage2x;
      if ~isa(featim,'uint8')
        featim = uint8(featim*255);
      end

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
      end
      
      
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
          plot(pcv, 'showPixelList', true, 'showEllipses', false);
        end
        title('openCV')
        keyboard
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
    
    %% NEED TO BE UPDATED TO REMAP!!!!
    function doimg = a_interp(self,foimg,xshift,mback,type)

      [h,w] = size(foimg);
      indw = single(0:w-1);
      indh = single(0:h-1)';
      indxI = max(min(bsxfun(@plus,indw,xshift),w),1);
      indy = repmat(indh,[1,w]);

      % color should work, too.
      doimg = cv.remap(foimg,indxI,indy,'BorderType','Constant','BorderValue',mback,'Interpolation',type);

    end
    %%

    
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
        rp(s,1).MajorAxisLength = max(info.size);
        rp(s,1).MinorAxisLength = min(info.size);
        rp(s,1).Orientation = -info.angle + 90;
        smallmask = bwimg(bb(2)+1:bb(2)+bb(4),bb(1)+1:(bb(1)+bb(3)));
        rp(s,1).Image =  smallmask;

        
        % COULD USE THIS FUNCTION: cv.getRectSubPix % but maybe overkill
        
        %[idxi,idxj] = find(smallmask);
        %startidx = bb([2,1]);
        %rp(s,1).PixelIdxList = (idxi + startidx(1)) + (idxj + (startidx(2)-1))*sizei;
        %imregion = Iframe(bb(2)+1:bb(2)+bb(4),bb(1)+1:(bb(1)+bb(3)));
        %if isuint
        %  rp(i).FilledImage = single(imregion);
        %else
        %  rp(i).FilledImage = imregion;
        %end
        %rp(i).MeanIntensity = mean(rp(i).FilledImage(rp(i).Image));
      end
    end
    
    
      
    function self = FishBlobAnalysisCV(varargin) % constructor
      
      self = self@FishBlobAnalysis(varargin{:});
    

    end
  
  end

end

      

    