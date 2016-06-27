

LOAD = 0;
COMPUTE = 1;
PLOT = 1;

if LOAD 
 load ~/data/zebra/videos/textures/texturesinregionsF4-1-1.mat
 ft.getTrackingResults([],[],1); % force to generate new DAG
end


if COMPUTE

  timeRange = [100,200];  
  dagif = 1;

  res = ft.getTrackingResults(timeRange,dagif);
  turn = ft.getTurningStats(res);

  % we first have to match the current fishID to the stimFishID at the time of
  % stimulus presentations
  stmInfo = res.tracks.stmInfo; % this is sorted accroding to gloval estimate
                                % of fishIDs (might include some ID switches)
  stmFishId = stmInfo(:,:,ft.stimulusPresenter.IDX_FISHID);


  
  % msk includes boundary off (same as ~isnan(bbox(1)))
  stmmsk = stmInfo(:,:,ft.stimulusPresenter.IDX_MSK); 

  
  stmstate = stmInfo(:,:,ft.stimulusPresenter.IDX_STATE); 
  %sf = stmInfo(:,:,ft.stimulusPresenter.IDX_SIZEFACTOR);
  %stmstate = stmstate.*(sf>0); % stmmsk sets sf to zero

  % NOTE: some of the stmstates correspond to sf==0. These can be used for
  % control comparison (same length statistics)
  
  % do not use stmmsk because occasonially has boundary effect 
  
  for i = 1:length(turn)
    stmtidx = find(diff(stmstate(:,i))==1) + 1;
    turn(i).stm.tidx = stmtidx;
  
    for f = {'IDX_SIZEFACTOR','IDX_SHIFTORI','IDX_SHIFT','IDX_XY','IDX_FISHID'}
      turn(i).stm.(lower(f{1}(5:end))) = stmInfo(stmtidx,i,ft.stimulusPresenter.(f{1}));
    end
    
  end

end
