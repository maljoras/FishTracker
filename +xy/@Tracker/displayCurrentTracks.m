function varargout = displayCurrentTracks(self)
% UFRAME = DISPLAYCURRENTTRACKS(SELF) generates a uframe from current
% tracking results.

  cjet = jet(self.nextId);

  minVisibleCount = 2;
  if self.opts.display.BWImg
    uframe = getCurrentBWImg(self.videoHandler);
  else
    uframe = getCurrentFrame(self.videoHandler);
  end
  
  if ~isa(uframe,'uint8')
    uframe = uint8(uframe*255);
  end

  if size(uframe,3)==1
    uframe = repmat(uframe,[1,1,3]);
  end


  if ~isempty(self.tracks)
    
    % Noisy detections tend to result in short-lived tracks.
    % Only display tracks that have been visible for more than
    % a minimum number of frames.
    reliableTrackInds = ...
        [self.tracks(:).totalVisibleCount] > minVisibleCount;
    reliableTracks = self.tracks(reliableTrackInds);
    
    % Display the objects. If an object has not been detected
    % in this frame, display its predicted bounding box.
    if ~isempty(reliableTracks)
      % Get bounding boxes.
      bboxes = cat(1, reliableTracks.bbox);
      
      % Get ids.
      ids = int32([reliableTracks(:).id]);

      % Create labels for objects indicating the ones for
      % which we display the predicted rather than the actual
      % location.
      labels = cellstr(int2str(ids'));
      predictedTrackInds = [reliableTracks(:).consecutiveInvisibleCount] > 0;
      isPredicted = cell(size(labels));
      isPredicted(predictedTrackInds) = {' predicted'};

      for i = 1:length(reliableTracks)
        if ~predictedTrackInds(i) & reliableTracks(i).assignmentCost>1
          isPredicted{i} = sprintf('  %1.0d',round(reliableTracks(i).assignmentCost));
        end
      end
      

      labels = strcat(labels, isPredicted);
      cols = jet(self.nanimals);
      cols_grey = [0.5,0.5,0.5;cols];
      
      
      if length(reliableTracks)==self.nanimals && self.isInitClassifier

        ids = [reliableTracks(:).identityId];

        pids = cat(1,reliableTracks.predIdentityId);

        
        clabels = uint8(cols(ids,:)*255);
        pclabels = uint8(cols(pids,:)*255);
        % Draw the objects on the frame.
        uframe = insertShape(uframe, 'FilledRectangle', bboxes,'Color',pclabels,'Opacity',0.2);
        uframe = insertObjectAnnotation(uframe, 'rectangle', bboxes, labels,'Color',clabels);

      else
        % Draw the objects on the frame inm grey 
        uframe = insertObjectAnnotation(uframe, 'rectangle', bboxes, labels,'Color', uint8([0.5,0.5,0.5]*255));
      end
      
      center = cat(1, reliableTracks.location);
      center = cat(2,center,max(self.fishlength/2,self.fishwidth)*ones(size(center,1),1));
      
      crossedTrackIdStrs = arrayfun(@(x)num2str(sort(x.crossedTrackIds)), reliableTracks,'uni',0);
      [u,idxct,idxu] = unique(crossedTrackIdStrs);
      for i =1:length(u)
        crossId = reliableTracks(idxct(i)).crossedTrackIds;
        [~,cross] =  ismember(crossId,[reliableTracks.id]);
        cross(~cross) = [];
        
        uframe = insertObjectAnnotation(uframe, 'circle', center(cross,:), 'Crossing!','Color', uint8(cols(i,:)*255));
      end
      
      
      %% insert some  features
      if self.opts.display.level>1
        for i = 1:length(reliableTracks)
          
          id = reliableTracks(i).identityId;
          if isnan(id)
            id = 1;
          end
          
          segm = reliableTracks(i).segment;
          
          
          if isfield(segm,'MSERregions') && ~isempty(segm.MSERregions)
            for j = 1:length(segm.MSERregions)
              pos = segm.MSERregionsOffset + segm.MSERregions(j).Location;
              uframe = insertMarker(uframe, pos, '*', 'color', uint8(cols(id,:)*255), 'size', 4);
            end
          end
          pos = reliableTracks(i).location;
          uframe = insertMarker(uframe, pos, 'o', 'color', uint8(cols(id,:)*255), 'size', 5);
          
          if ~isempty(reliableTracks(i).centerLine) && ~any(isnan(reliableTracks(i).centerLine(:)))
            uframe = insertMarker(uframe, reliableTracks(i).centerLine, 's', 'color', uint8([1,1,1]*255), 'size', 2);
          end
          
        end
      end
      
      
      
      if self.opts.display.level>2
        %% insert more markers
        if length(self.tracks)==self.nanimals
          howmany = 25;
          idx = max(self.currentFrame-howmany,1):self.currentFrame;
          trackpos = self.pos(:,:,idx);
          f2t = self.identityId2TrackId(idx,:);
          delidx = find(any(any(isnan(trackpos),1),2));
          trackpos(:,:,delidx) = [];
          f2t(delidx,:) = [];
          cli = NaN(length(idx),self.nanimals);
          if ~isempty(self.tracks(1).classProbHistory)

            for iii = 1:length(self.tracks)
              [classProbs, w]= self.tracks(iii).classProbHistory.getData(length(idx));
              [~,classIdx] = max(classProbs,[],2);
              msk = any(isnan(classProbs),2) | w < self.tracks(iii).classProbHistory.reasonableThres;
              classIdx(msk) = 0;
              cli(end-size(classIdx)+1:end,iii) = classIdx;
            end

          end
          cli(delidx,:) = [];
          cli = cat(2,ones(size(cli,1),1),cli+1);
          
          trackIds = [self.tracks.id];
          [~,f2i] = ismember(f2t,trackIds);

          if ~isempty(trackpos)
            for ii = 1:self.nanimals
              idx1 = f2i(:,ii)+1;
              inds = xy.helper.s2i(size(cli),[(1:size(cli,1))',idx1]);
              inds2 = cli(inds);
              pos1 = squeeze(trackpos(:,ii,~isnan(inds2)))';
              cols1 = uint8(cols_grey(inds2(~isnan(inds2)),:)*255);
              cols2 = uint8(cols(ii,:)*255);
              uframe = insertMarker(uframe, pos1, 'o', 'color', cols2 , 'size', 3);
              uframe = insertMarker(uframe, pos1, 'x', 'color', cols1 , 'size', 2);
            end
          end
        end
      end
      
    end
  end
  if ~nargout
    self.videoPlayer.step(uframe);
    self.stepVideoWriter(uframe);
  else
    varargout{1} = uframe;
  end
end
