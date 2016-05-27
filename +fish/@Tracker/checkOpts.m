function checkOpts(self)
% checks (and modifies potentially) the given options. 

  
  
  %% look for objects 
  if isempty(self.nfish) || isempty(self.fishlength) || isempty(self.fishwidth) || self.useScaledFormat  
    [nfish,fishSize] = self.findObjectSizes();
    
    
    if (nfish>25 || ~nfish) && isempty(self.nfish)% too many... something wrong
      fish.helper.verbose('WARNING: The fish size and number cannot be determined')
      if self.displayif && self.opts.display.fishSearchResults
        self.nfish = fish.helper.chooseNFish(vid,1); % only if interactively
        fishSize = [100,20]; % wild guess;
        close(gcf);
      else
        error('Please manual provide fishlength, fishwidth and nfish');
      end
    end

    
    if isempty(self.opts.fishlength) % otherwise already set by hand
      self.fishlength = fishSize(1);
    end
    if isempty(self.opts.fishwidth) 
      self.fishwidth = fishSize(2); 
    end
    if isempty(self.opts.nfish) 
      self.nfish = nfish;
    end

  end  
  assert(self.fishlength>self.fishwidth);      
  % overwrite the given options
  self.fishlength = self.fishlength;

  
  %% classifier
  if self.opts.classifier.nlfd
    if self.opts.classifier.nlfd<0
      nlfd = self.nfish+1;
      fish.helper.verbose('Set nlfd to nfish+1 [%d]',nlfd);
      self.opts.classifier.nlfd = nlfd;
      self.opts.classifier.npca = max(self.opts.classifier.npca,nlfd+1);
      
    elseif self.opts.classifier.nlfd<self.nfish+1;
      nlfd = self.nfish+1;
      fish.helper.verbose('Nlfd too small. Increased nlfd to number of fish+1 [%d]',nlfd);
      self.opts.classifier.nlfd = nlfd;
      self.opts.classifier.npca = max(self.opts.classifier.npca,nlfd+1);
      
    end
  end
end


