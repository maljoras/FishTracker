function checkOpts(self)
% checks (and modifies potentially) the given options. 

  
 
  
  
  assert(self.fishlength>self.fishwidth);      
  % overwrite the given options
  self.fishlength = self.fishlength;

  
  %% classifier
  if self.opts.classifier.nlfd
    if self.opts.classifier.nlfd<0
      nlfd = max(self.nfish+1,10);
      fish.helper.verbose('Set nlfd [%d]',nlfd);
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


