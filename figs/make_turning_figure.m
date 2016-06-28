

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
  
  
  velocity = [zeros(1,2,ft.nfish);diff(ft.interpolateInvisible(res,'pos',3))];

  
  % msk includes boundary off (same as ~isnan(bbox(1)))
  stmmsk = stmInfo(:,:,ft.stimulusPresenter.IDX_MSK); 

  
  stmstate = stmInfo(:,:,ft.stimulusPresenter.IDX_STATE); 
  %sf = stmInfo(:,:,ft.stimulusPresenter.IDX_SIZEFACTOR);
  %stmstate = stmstate.*(sf>0); % stmmsk sets sf to zero

  % NOTE: some of the stmstates correspond to sf==0. These can be used for
  % control comparison (same length statistics)
  
  % do not use stmmsk because occasonially has boundary effect 
  
  for i = 1:length(turn)
    stmstatei = stmstate(:,i);
    stmtidx_start = find(diff([0;stmstatei])==1);
    stmtidx_stop = find(diff([stmstatei;0])==-1);
    turn(i).stm.tidx = stmtidx_start;

    nstm = length(stmtidx_start);
    turn(i).stm.stoptidx = stmtidx_stop;

    accummsk = zeros(size(stmstate,1),1);
    accummsk(stmtidx_start) = 1;
    accummsk = cumsum(accummsk);
    accummsk(~accummsk) = nstm+1;

    accummsk_nogap = accummsk;
    turn(i).stm.accummsk_nogap = accummsk_nogap;

    accummsk_gap  =  accummsk;
    accummsk_gap(logical(stmstatei)) = nstm+1;
    turn(i).stm.accummsk_gap = accummsk_gap;
    
    accummsk(~stmstatei) = nstm+1;
    turn(i).stm.accummsk = accummsk;

    for f = {'IDX_SIZEFACTOR','IDX_SHIFTORI','IDX_SHIFT','IDX_XY','IDX_FISHID'}
      turn(i).stm.(lower(f{1}(5:end))) = stmInfo(stmtidx_start,i,ft.stimulusPresenter.(f{1}));
    end
    turn(i).stm.acc_fishId = accumarray(accummsk,stmFishId(:,i)==i,[nstm+1,1],@mean);
    turn(i).stm.acc_len = accumarray(accummsk,1,[nstm+1,1],@sum);

    
    turnmsk = zeros(size(accummsk));
    turnmsk(turn(i).tidx) = 1;
    turn(i).stm.acc_nturns = accumarray(accummsk,turnmsk,[nstm+1,1],@sum);

    turnmsk(turn(i).tidx) = 1:length(turn(i).tidx);
    turnmsk(~turnmsk) = NaN;
    turn(i).stm.acc_firstturn_idx = accumarray(accummsk,turnmsk,[nstm+1,1],@nanmin);
    turn(i).stm.acc_lastturn_idx = accumarray(accummsk,turnmsk,[nstm+1,1],@nanmax);
    turn(i).stm.acc_gapfirstturn_idx = accumarray(accummsk_gap,turnmsk,[nstm+1,1],@nanmin);
    turn(i).stm.acc_gaplastturn_idx = accumarray(accummsk_gap,turnmsk,[nstm+1,1],@nanmax);

    turn(i).stm.acc_nogapfirstturn_idx = accumarray(accummsk_nogap,turnmsk,[nstm+1,1],@nanmin);
    turn(i).stm.acc_nogaplastturn_idx = accumarray(accummsk_nogap,turnmsk,[nstm+1,1],@nanmax);
  end

end
