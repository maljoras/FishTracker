function [nObjects,objectSize] =findObjectSizes(self,minAxisWidth)

  if nargin==1
    minAxisWidth = 4;
  end

  verbose('Detect approx fish sizes (minAxisWidth=%g)...',minAxisWidth);

  self.videoHandler.reset();
  self.videoHandler.originalif = true;
  n = min(self.videoHandler.history,floor(self.videoHandler.timeRange(2)*self.videoHandler.frameRate));
  self.videoHandler.initialize(0);
  s = 0;
  for i = 1:n

    [segm] = self.videoHandler.step();
    fprintf('%1.1f%%\r',i/n*100); % some output

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
    error('Something wrong with the bolbAnalyzer ... cannot find objects.');
  end

  nObjects = median(count);
  objectSize = ceil([quantile(height,0.9),quantile(width,0.9)]); % fish are bending
  verbose('Detected %d fish of size %1.0fx%1.0f pixels.',nObjects, objectSize);

  if self.displayif && self.opts.display.fishSearchResults
    figure;
    bwmsk = self.videoHandler.getCurrentBWImg();
    imagesc(bwmsk);
    title(sprintf('Detected %d fish of size %1.0fx%1.0f pixels.',nObjects, objectSize));
  end

  if self.useScaledFormat
    % get good color conversion
    bwmsks = {};
    cframes = {};
    for i =1:10 % based on not just one frame to get better estimates
      [segm] = self.videoHandler.step();
      bwmsks{i} =self.videoHandler.getCurrentBWImg();
      cframes{i} = im2double(self.videoHandler.getCurrentFrame());
    end
    
    [scale,delta] = getColorConversion(bwmsks,cframes);
    
    if ~isempty(scale)
      self.videoHandler.setToScaledFormat(scale,delta);
      self.videoHandler.resetBkg(); 
      self.videoHandler.computeSegments = false;
      % re-generate some background
      verbose('Set to scaled format. Regenerate background..')
      for i =1:(n/2)
        [~,frame] = self.videoHandler.step();    
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
