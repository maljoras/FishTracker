function checkOpts(self)
% checks (and modifies potentially) the given options. 

  assert(self.bodylength>self.bodywidth,'Assertion error: bodylength>bodywidth');      
  % overwrite the given options
  self.bodylength = self.bodylength;

  
  if self.opts.tracks.useDagResults
    if self.bodylength<80
      warning(['Using DAG results for small bodysizes might result ' ...
               'in inferior results. Set default to SWB.']);
      self.opts.tracks.useDagResults = 0;
    end
    if self.nindiv>10
      warning(['Using DAG results for large number of individials might result ' ...
               'in inferior results. Set default to SWB.']);
      self.opts.tracks.useDagResults = 0;
    end
  end
  
  if self.nindiv>10 && self.bodylength<80 &&   self.opts.classifier.allSwitchProbThres<0.5
    warning(['Many small individuals. Consider increasing ' ...
             '"classifier.allSwitchProbThres"']);
    self.opts.classifier.allSwitchProbThres = 0.6;
  end
  
    
    
  %% classifier
  if self.opts.classifier.nlfd
    if self.opts.classifier.nlfd<0
      nlfd = max(self.nindiv+1,10);
      xy.helper.verbose('Set nlfd [%d]',nlfd);
      self.opts.classifier.nlfd = nlfd;
      self.opts.classifier.npca = max(self.opts.classifier.npca,nlfd+1);
      
    elseif self.opts.classifier.nlfd<self.nindiv+1;
      nlfd = self.nindiv+1;
      xy.helper.verbose('Nlfd too small. Increased nlfd to number of individuals+1 [%d]',nlfd);
      self.opts.classifier.nlfd = nlfd;
      self.opts.classifier.npca = max(self.opts.classifier.npca,nlfd+1);
      
    end
  end
end


