classdef PresenterTrackTextureRegions < xy.stimulus.PresenterTrackTextures;
  
  properties
    xRegions = [0.5];
    regSizeFactorScale = [0.5,10];
  end

  properties(SetAccess=private,GetAccess=private)
  
    xRegHelperUpper = [];
    xRegHelperLower = [];
  end
  
  
  methods
  
    
    function reset(self)
      reset@xy.stimulus.PresenterTrackTextures(self);
    
      self.xRegHelperUpper = [self.xRegions(:)',1];
      self.xRegHelperLower = [0,self.xRegions(:)'];
    end

    
    function [stmSizeFactor] = assignStmSizeFactor(self,x,y,t,identityIds)
      
      rmsk = bsxfun(@ge,self.xRegHelperUpper,x(:));
      rmsk = rmsk & bsxfun(@le,self.xRegHelperLower,x(:));
      [msk,regidx] = max(rmsk,[],2);

      szfac = self.stmSizeFactor(min(identityIds,end));
      rsf = self.regSizeFactorScale(:);
      szfac = szfac(:).*rsf(regidx);
      stmSizeFactor = self.grand(szfac,self.stmSizeFactorCV);     
      stmSizeFactor = stmSizeFactor.*msk;
    end

  end
  
end




