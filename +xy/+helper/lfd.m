function [FV,eVar,Proj,r] = lfd(Z,iidx,nbody,npca,whiteif)
% [FV,EVAR,PROJ,RSTRUC] = LFD(Z,IIDX,NBODY,NPCA) computes the LInear Fisher
% Discrimant directions of Z (axyTer initial PCA with NPCA components,)
% using generalized eigenvectors. IIDX is has a unique number for each
% individual (like class labels). NPCA=0 shuts the PCA preprocessing
% off.
% LFD(..,WHITEIF) whites the data axyTer PCA (default 0). 

  
  %first do pca on Z
  NPCADIM = npca;
  N = size(Z,1);
  Nc = length(unique(iidx));

  if NPCADIM>0
    npcadim = min(N-Nc,NPCADIM);
    [PC, eVar, Proj, m] = xy.helper.pca1(Z,npcadim);
    X = Proj;
  else
    X = Z;
  end

  r.totalVar = sum(var(Z));
  r.retainedVar = sum(var(X))/r.totalVar;

  if exist('whiteif','var') && whiteif
    % AXYTER Pca
    X = whiten(X); % use function from matlab central; 
  end
    
    
    
  
  noclassnorm=1;% all classes are weighted the same regardless the
                   % number of samples

  [u dummy ui]= unique(iidx);  
  Nc = length(u);
  Ni = accumarray(ui,ui,[],@numel);
  

  mX = zeros(Nc,size(X,2));
  for i = 1:size(X,2)
    mX(:,i) = accumarray(ui,X(:,i),[],@mean);
  end
  
  
  %within class and between class covariance
  mu = mean(X,1);
  x = bsxfun(@minus,mX,mu); 

  Sb = 0;
  for k = 1:Nc % just do it in a loop...
    if noclassnorm 
      Sb= Sb + x(k,:)'*x(k,:); 
    else
      Sb= Sb + Ni(k)*x(k,:)'*x(k,:); 
    end
    
  end

  Sw = 0;
  for i = 1:Nc
    j = find(iidx==u(i));
    x = bsxfun(@minus,X(j,:),mX(i,:));
    for k = 1:size(x,1)
      if noclassnorm
        Sw = Sw + 1/size(x,1)*x(k,:)'*x(k,:);
      else
        Sw = Sw + x(k,:)'*x(k,:);
      end
      
    end
  end
  
  % generalized eigenvalue problem
  [V,D] = eig(Sb,Sw);

  % norm V !!
  for i = 1:size(V,2);
    V(:,i) = V(:,i)/norm(V(:,i));
  end

  d = diag(D);
  [sd,sidx] = sort(d,'descend');
  
  
  r.D = sd;
  r.V = V(:,sidx);
  eVar = sd(1:nbody);
  r.Sb = Sb;
  r.Sw = Sw;
  
  % save other stuff
  r.mX = mX;
  r.mu = mu;

  r.Nc = Nc;
  r.Ni = Ni;
  
  if npca>0
    r.m = m';
    r.PC = PC;
  end
  
  
  
  Proj = bsxfun(@minus,X,mu)*r.V(:,1:nbody);
  FV = r.V(:,1:nbody);
  
