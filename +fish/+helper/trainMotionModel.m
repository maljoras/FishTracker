function trainMotionModel(self,fishIds,timeRange);

  if ~exist('fishIds','var') || isempty(fishIds)
    fishIds = 1:self.nfish;
  end

  if ~exist('timeRange','var') || isempty(timeRange)
    timeRange = self.timerange;
  end

  assert(isa(self,'fish.Tracker'));

  centerLine = self.deleteInvisible('centerLine',timeRange);
  centerLine = centerLine(:,fishIds,:,:);
  velocity = self.deleteInvisible('velocity',timeRange);
  velocity = velocity(:,fishIds,:);
  

  

  orgX = permute(centerLine,[1,4,3,2]);
  headpos = orgX(1:end-1,1,:,:);
  X = bsxfun(@minus,orgX(1:end-1,:,:,:),headpos);

  % rotate in direction to velocity
  Xr = zeros(size(X));
  vori = atan2(velocity(:,:,1),velocity(:,:,2));
  vori = permute(vori(1:end-1,:),[1,3,4,2]);
  
  Xr(:,:,1,:) = bsxfun(@times,cos(vori),X(:,:,1,:)) - bsxfun(@times,sin(vori),X(:,:,2,:));
  Xr(:,:,2,:) = bsxfun(@times,cos(vori),X(:,:,1,:)) + bsxfun(@times,sin(vori),X(:,:,2,:));
  
  
  % build training set
  Xr(:,1,:,:) = permute(velocity(1:end-1,:,:),[1,4,3,2]); 
  Y = bsxfun(@minus,orgX(2:end,:,:,:),headpos);

  
  keyboard
  
  % make model
  sz = size(Y);
  predX = zeros(sz);
  
  for i = 1:size(X,4);
    x = X(:,:,:,i);
    y = Y(:,:,:,i);
    x1 = cat(2,x(:,:),ones(size(x,1),1));
    [beta{i},sig{i},~,~,llh] = mvregress(x1,y(:,:));

    predy = reshape(x1*beta{i},sz(1:3));
    predY(:,:,:,i) = predy;

  end
  predY = bsxfun(@plus,predY,headpos);
  clf;
  
  
  res = self.getTrackingResults();
  t = res.t(:,1);
  ptr = [10,20];
  idx = find(t>=ptr(1) & t<ptr(2));
  idx = idx(1:end);
  
  orgX = orgX(2:end,:,:,:);
  
  cmap = jet(self.nfish);
  for i = 1:length(fishIds)
    col = cmap(fishIds(i),:);
    
    plot(orgX(idx,:,1,i)',orgX(idx,:,2,i)','color',col,'linewidth',2);
    hold on;        
    plot(orgX(idx,1,1,i),orgX(idx,1,2,i),'o','color',col,'linewidth',1,'markersize',4);
    
    % prediction
    plot(predY(idx,:,1,i)',predY(idx,:,2,i)','color','k','linewidth',1);
    hold on;        
    plot(predY(idx,1,1,i),predY(idx,1,2,i),'s','color','k','linewidth',1,'markersize',4);
  end

  for i = 1:100;
    cla;
    plot(predY(i,:,1,1),predY(i,:,2,1),'r');
    hold on;
    plot(orgX(i+1,:,1,1),orgX(i+1,:,2,1),'b');
    plot(orgX(i,:,1,1),orgX(i,:,2,1),'k');
    drawnow;
    pause;
  end
  
  %save(fileName,'X','Y');
  
  