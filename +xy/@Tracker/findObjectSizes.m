function [nObjects,objectSize] =findObjectSizes(self,minAxisWidth)

  if nargin==1
    minAxisWidth = 3;
  end

  xy.helper.verbose('Detect approx body sizes (minAxisWidth=%g)...',minAxisWidth);

  self.videoHandler.reset();
  self.videoHandler.originalif = true;
  
  if numel(self.useScaledFormat)==3 || numel(self.useScaledFormat)==6;
    % THIS IS JUST AN INITIAL GUESS
    scale = self.useScaledFormat(1:3);
    if numel(self.useScaledFormat)==6
      delta = self.useScaledFormat(4:6);
    else
      mu = [0.5,0.5,0.5];
      delta = -scale.*mu + 0.5/3;
    end
    self.videoHandler.setToScaledFormat(scale,delta);
    xy.helper.verbose('Scale initially RGB format with %1.1f %1.1f %1.1f',...
                        scale(1),scale(2),scale(3))    
  elseif numel(self.useScaledFormat)~=1
    error('useScaledFormat: Either provide scale RGB + delta values or set to true/false');
  end
  
  n = min(floor(self.videoHandler.history/2),floor(self.videoHandler.timeRange(2)*self.videoHandler.frameRate));
  n = max(min(n,500),10);
  self.videoHandler.initialize(0);
  s = 0;

  if self.displayif && self.opts.display.bodySearchResults
    figure;
  end
  
  for i = 1:n
    if ~self.videoHandler.hasFrame()
      break
    end
    
    [segm] = self.videoHandler.step();
    if ~(mod(i,50)) 
      xy.helper.verbose('%1.1f%%\r',i/n*100); % some output
    end
    
    if self.displayif && self.opts.display.bodySearchResults && ~mod(i,5)
      imagesc(self.videoHandler.getCurrentBWImg);
      drawnow;
    end
    
    if isempty(segm)
      continue;
    end
    s = s+1;
    idx = [segm.MinorAxisLength]>minAxisWidth;
    count(s) = length(segm(idx));
    width(s) = median([segm(idx).MinorAxisLength]);
    height(s) = median([segm(idx).MajorAxisLength]);
  end

  if ~s
    error('Something wrong with the blobAnalyzer ... cannot find objects.');
  end

  nObjects = median(count(1:s));
  objectSize = ceil([quantile(height(1:s),0.9),quantile(width(1:s),0.9)]); % body are bending
  xy.helper.verbose('Detected %d objects of size %1.0fx%1.0f pixels.',nObjects, objectSize);

  if self.displayif && self.opts.display.bodySearchResults
    clf;
    bwmsk = self.videoHandler.getCurrentBWImg();
    imagesc(bwmsk);
    title(sprintf('Detected %d objects of size %1.0fx%1.0f pixels.',nObjects, objectSize));
  end


  if any(self.useScaledFormat)

      % get good/better color conversion
      
      bwmsks = {};
      cframes = {};
      
      for i =1:10 % based on not just one frame to get better estimates
        [segm] = self.videoHandler.step();
        bwmsks{i} =self.videoHandler.getCurrentBWImg();
        cframes{i} = im2double(self.videoHandler.getCurrentFrame());
      end
      [scale,delta] = xy.helper.getColorConversion(bwmsks,cframes);

    
      
    if ~isempty(scale)
      self.videoHandler.setToScaledFormat(scale,delta);
      self.videoHandler.resetBkg(); 
      self.videoHandler.computeSegments = false;
      % re-generate some background
      xy.helper.verbose('Scale to scaled RGB format: %1.1f %1.1f %1.1f',...
                          scale(1),scale(2),scale(3))    

      xy.helper.verbose('Regenerate background.')
      for i =1:(n/2)
        [~,~] = self.videoHandler.step();    
        fprintf('%1.1f%%\r',i/n*2*100); % some output
      end
      self.videoHandler.computeSegments = true;
    else
      self.videoHandler.setToRGBFormat();
      self.videoHandler.frameFormat = ['GRAY' self.videoHandler.frameFormat(end)];
    end
  end

  self.videoHandler.originalif = false;
  self.videoHandler.reset();
end
