classdef Tracker
% dummy for loading the old fishTracker objects into new xyTracker format
  
  
  methods(Static)
    
    function obj = loadobj(S);
      S.opts.nindiv = S.nfish;
      S.opts.bodylength = S.fishlength;
      S.opts.bodywidth = S.fishwidth;
      S.opts = rmfield(S.opts,{'fishwidth','fishlength','nfish'});
      
      S.nindiv = S.nfish;
      S.identityClassifier = S.fishClassifier;
      S.bodylength = S.fishlength;
      S.bodywidth = S.fishwidth;
      S.uniqueIdentityFrames = S.uniqueFishFrames;
      S.identityId2TrackId = S.fishId2TrackId;
      S = rmfield(S,{'fishwidth','fishlength','nfish','fishClassifier','uniqueFishFrames','fishId2TrackId'});      

      S.opts.tracks = renameStructField(S.opts.tracks,'probThresForFish','probThresForIdentity');
      % rename fishId fields
      if isfield(S.tracks,'fishId')
        [S.tracks.identityId] = deal(S.tracks.fishId);
        S.tracks = rmfield(S.tracks,'fishId');
      end
      
      if isfield(S.tracks,'predFishId')
        [S.tracks.predIdentityId] = deal(S.tracks.predFishId);
        S.tracks = rmfield(S.tracks,'predFishId');
      end

      for f = {'savedTracks','saveFieldsOut','saveFieldsIn','savedTracksFull'}
        if isfield(S.(f{1}),'fishId')
          S.(f{1}) = renameStructField(S.(f{1}),'fishId','identityId');
        end
        if isfield(S.(f{1}),'predFishId')
          S.(f{1}) = renameStructField(S.(f{1}),'predFishId','predIdentityId');
        end
      end
      
      obj = xy.Tracker.loadobj(S);
      obj.getTrackingResults([],[],1); % force generate new res;
    end
  end
end
