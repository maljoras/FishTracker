function plotWrongAssignments(self,idx)
% plotting helper function
  clf;
  imagesc(self.videoHandler.getCurrentFrame());
  hold on;

  detectionIdx = self.assignments(:, 2);
  trackIdx = self.assignments(:, 1);

  centroids = self.centroids(detectionIdx, :);
  latestcentroids = cat(1,self.tracks.centroid);
  plot(latestcentroids(:,1),latestcentroids(:,2),'<k');
  plot(centroids(:,1),centroids(:,2),'or');
  plot(self.centroids(:,1),self.centroids(:,2),'xm');

  plot(latestcentroids(idx,1),latestcentroids(idx,2),'.w');


  bboxes = self.bboxes(detectionIdx, :);
  latestbboxes = cat(1,self.tracks.bbox);

  for i = 1:size(latestbboxes,1)
    rectangle('position',latestbboxes(i,:),'edgecolor','k');
  end
  for i = 1:size(bboxes,1)
    rectangle('position',bboxes(i,:),'edgecolor','r');
  end
end

