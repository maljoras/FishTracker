function plotCrossings_(self,section,updateIdx)


  figure(2);
  if section==1
    clf;
  end

  cols = jet(self.nfish);
  [r1,r2] = fish.helper.getsubplotnumber(length(updateIdx));

  s = 0;
  for trackIndex = updateIdx
    s = s+1;
    a = subplot(r1,r2,s);
    set(a,'ColorOrder',cols,'ColorOrderIndex',1);
    hold on;
    
    t1 = self.tracks(trackIndex).firstFrameOfCrossing;
    t2 = self.tracks(trackIndex).lastFrameOfCrossing;
    t0 = max(t1-30,1);
    t = self.currentFrame;
    
    if section==1
      ids = self.tracks(trackIndex).crossedTrackIds;
      
      %first plot all traces
      plot(squeeze(self.pos(1,:,t0:t))',squeeze(self.pos(2,:,t0:t))','--','Linewidth',1);
      hold on; 
      
      % plot those that are in the current crossing
      fids = find(ismember(self.fishId2TrackId(t2,:),ids));
      set(a,'ColorOrder',cols(fids,:));set(a,'ColorOrderIndex',1)
      plot(squeeze(self.pos(1,fids,t1:t2))',squeeze(self.pos(2,fids,t1:t2))','o','Linewidth',2);
      
    else
      %section 2
      fid = self.tracks(trackIndex).fishId;
      plot(squeeze(self.pos(1,fid,t0:t))',squeeze(self.pos(2,fid,t0:t))','-','Linewidth',2,'color',cols(fid,:));
      frameSize = self.videoHandler.frameSize;
      axis ij;
      xlim([0,frameSize(2)])
      ylim([0,frameSize(1)]);
      

    end
    
  end
  if section==2
    drawnow;
    
    self.daGraph.plotTraceAssignments([self.tracks.fishId],self.assignments);
  end
end
