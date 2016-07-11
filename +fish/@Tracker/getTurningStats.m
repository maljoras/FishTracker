function turn = getTurningStats(self,res)
% TURNSTATS = GETTURNINGSTATS(SELF,RES) gets some stats on the
% turnings from RES

  deltat = max(floor(self.avgTimeScale),1); % frames/BL
  deltat = deltat + mod(deltat-1,2);

  turn = self.getTurningPoints(res);
  velocity = self.getResField(res,'velocity',1);
  vel = sqrt(sum(velocity.^2,3));
  
  pos = self.getResField(res,'pos',1);
  locvel = permute(bsxfun(@rdivide,sqrt(sum(diff(pos).^2,2)),diff(res.tabs)),[1,3,2]);
  locvel(end+1,:) = NaN;
  locvel(locvel>quantile(locvel(:),0.99))=NaN;
  locvel = fish.helper.movavg(locvel,max(round(deltat/2),1));
  

  mxt = deltat*10;
  for i = 1:size(velocity,2)
    tidx = turn(i).tidx;

    turn(i).vel0 = vel(tidx,i);
    turn(i).vel1 = vel(min(tidx+deltat,end),i);
    turn(i).vel2 = vel(min(tidx+2*deltat,end),i);
    turn(i).acc1 = turn(i).vel1-turn(i).vel0;
    turn(i).acc2 = turn(i).vel2-turn(i).vel1;
    
    turn(i).locvel0 = locvel(tidx,i);
    turn(i).locvel1 = locvel(min(tidx+deltat,end),i);
    turn(i).locvel2 = locvel(min(tidx+2*deltat,end),i);

    turn(i).len = sqrt(diff(turn(i).x).^2 + diff(turn(i).y).^2);
    turn(i).len(end+1) = NaN;
    turn(i).dori = angle(exp(1i*diff(turn(i).ori)));
    turn(i).dori(end+1) = NaN;
      
    accmsk = ones(size(vel,1),1);
    accmsk(1:tidx(1)-1) = 0;
    accmsk(tidx(2:end)) = -diff(tidx)+1;
    accmsk = cumsum(accmsk);
    accmsk(~accmsk | accmsk>mxt) = mxt+1;
    
    
    % avg turn-triggered veocity % WRONG ! HAVE TO ACCOUNT FOR NON-EQUAL DT!!
    turn(i).attv = accumarray(accmsk,vel(:,i),[mxt+1,1],@nanmean);
    turn(i).attv = turn(i).attv(1:end-1);
    turn(i).sattv = accumarray(accmsk,vel(:,i),[mxt+1,1],@fish.helper.stderr);
    turn(i).sattv = turn(i).sattv(1:end-1);
    
     
    turn(i).attlocv = accumarray(accmsk,locvel(:,i),[mxt+1,1],@nanmean);
    turn(i).attlocv = turn(i).attlocv(1:end-1);
    turn(i).sattlocv = accumarray(accmsk,locvel(:,i),[mxt+1,1],@fish.helper.stderr);
    turn(i).sattlocv = turn(i).sattlocv(1:end-1);

  end
  