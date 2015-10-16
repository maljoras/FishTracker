

fnames = {'~/data/zebra/videos/ABfish/fishA.avi',['~/data/zebra/' ...
                    'videos/ABfish/Bfish.avi'],'~/data/zebra/videos/ABfish/AandBfish.avi'};
nfish = [1,1,2];
tr = [0,600]; % 10 minutes

for i = 1:length(fnames)
  ft =  FishTracker(fnames{i},'nfish',nfish(i),'displayif',0,'detector.thres',-0.044,'fishlength',110,'fishwidth',35);
  ft.saveFields = {ft.saveFields{:} 'segment.FilledImageCol2x','segment.fishFeature','segment.Image2x'};

  ft.timerange = tr;
  
  ft.track();
  save([fnames{i}(1:end-3) 'mat'],'ft','-v7.3')
end
