% comparse the FT and ID tracker

LOAD = 1;
COMPUTE = 0;
PLOT = 1;

if LOAD
  

  addpath('~/work/progs/toolboxes/mexopencv');
  addpath('~/work/progs/bnu/fish/FishTracker/helper');
  vid = '/home/malte/Videos/5Zebrafish_nocover_22min.avi';
  id = load('~/work/projects/zebra/trackering/trajectories.mat');
  idres.pos = permute(id.trajectories,[1,3,2]);
end

if COMPUTE
  ft = FishTracker(vid,'displayif',0);
  ft.track();
  ftres = ft.getTrackingResults();
  ft.save('compareFTID');
end



if PLOT
  nfish = ft.nfish;
  ftres = ft.getTrackingResults();
  dist = zeros(nfish);

  for i = 1:nfish
    for j = 1:nfish
      dist(i,j) = nanmean(sqrt(sum((ftres.pos(1:1000,:,i) - idres.pos(1:1000,:,j)).^2,2)));
    end
  end
  assignments = assignDetectionsToTracks(dist,1e3);

  idpos = idres.pos(1:end-1,:,:);
  ftpos = ftres.pos(:,:, assignments(assignments(:,1),2));
  
  
  %distances
  clf;
  r1 = nfish;
  r2 = 1;
  t = cat(1,ftres.tracks(:,1).t);
  for i = 1:nfish
    subplot(r1,r2,i);
    plot(t,sqrt(sum((ftpos(:,:,i) - idpos(:,:,i)).^2,2)))
    title(i);
  end
  
end
