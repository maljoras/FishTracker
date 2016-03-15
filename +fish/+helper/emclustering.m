function [pcl,a,mu,Sigma] = emclustering(X,k,plotif)
% [PCL,A,MU,SIGMA] = EMCLUSTERING(X,K) clusters data X  (n x d) while assuming K clusters with a
% Guassian mixture model

 n = size(X,1);
 d = size(X,2);

 maxiter = 100; 
 tol = 1e-10;
 reg = 0.001;
 plotif = exist('plotif','var');
 lfdif = 1; % for plotting only


 maxd = 10;
 if d>maxd
   % preprocess;
   [prepc, ~,X, prempc] = pca1(X,maxd); % on
   d = maxd;
 end

 % model parameters
 prob = rand(n,k);
 %Sigma = repmat(std(X(:))/k*eye(d),[1,1,k]);
 %mu = repmat(mean(X,1)',[1,k]);
 a = ones(1,k);
 oldtol = 0;
  
 % iteration step
 for i = 1:maxiter

   
   proba =  bsxfun(@times,prob,a);    
   pcl =  bsxfun(@rdivide,proba,sum(proba,2));     
   
   % new a
   a = mean(pcl,1);
   a = a/sum(a);

   % new mean
   sumcl = sum(pcl,1);
   mu = bsxfun(@rdivide,X' * pcl,sumcl);
   
   %var (already use new mu)
   for j = 1:k
     xm = bsxfun(@minus,X,mu(:,j)');
     Sigma(:,:,j) = xm'*(bsxfun(@times,xm,pcl(:,j)))/sumcl(j)+ reg*eye(d);
   end
   
   % compute tolerance
   newtol =0;
   newtol = mean(Sigma(:)) + mean(mu(:)) + mean(a);
   toldiff = abs(newtol-oldtol);
   oldtol = newtol;

   if plotif
     % plot first two dimensions
     clf;
     
     map = jet(k);
     if lfdif
       [~,cl] = max(pcl,[],2);
       [pc,~,proj,rstruc] = lfd(X,cl,2,0);
       mpc = rstruc.mu(:);
     else
       [pc, ~, proj, mpc] = pca1(X,2); % only for plotting
     end
     
     for j = 1:k
       plot(proj(cl==j,1),proj(cl==j,2),'.','Markersize',10,'linewidth',2, ...
            'color',map(j,:));
       hold on;
       plotGauss(pc'*bsxfun(@minus,mu(:,j),mpc),pc'*Sigma(:,:,j)*pc,'color',map(j,:));
     end
     drawnow;
   end
   

   if toldiff<tol
     break
   end
   
   % compute the class probs according to weighted distance
   for j = 1:k
     xm = bsxfun(@minus,X,mu(:,j)');
     xms = xm*inv(Sigma(:,:,j));
     prob(:,j) = exp(-sum(xms.*xm,2)/2)/sqrt(det(Sigma(:,:,j)))/(2*pi)^(d/2);
   end

   
 end
 
   
   
   
 
 
 
 function h = plotGauss(mu,Sigma,varargin);
% adapted from Mark A. Paskin
  k = 1;
  n = 100;
  [V, D] = eig(Sigma);
  
  % Compute the points on the surface of the ellipse.
  t = linspace(0, 2*pi, n);
  u = [cos(t); sin(t)];
  w = (k * V * sqrt(D)) * u;
  z = repmat(mu, [1 n]) + w;
  h = plot(z(1, :), z(2, :),varargin{:});