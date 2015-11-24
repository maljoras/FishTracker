

fnames = {'/data/videos/othervideos/ABfishforAppa/fishA.avi','/data/videos/othervideos/ABfishforAppa/Bfish.avi','/data/videos/othervideos/ABfishforAppa/AandBfish.avi'};
nfish = [1,1,2];
tr = [0,1200]; % 20 minutes

for i = 1:length(fnames)
  ft =  FishTracker(fnames{i},'nfish',nfish(i),'displayif',0,'blob.colorfeature',true,'useScaledFormat',1,'useMex',1,'fishlength',110,'fishwidth',25);
  ft.saveFields = {ft.saveFields{:},'segment.MinorAxisLength','segment.MajorAxisLength', 'segment.FishFeatureC','segment.FishFeature','segment.FishFeatureCRemap'};

  ft.timerange = tr;
  
  ft.track();
  save([fnames{i}(1:end-3) 'mat'],'ft','-v7.3')
end
