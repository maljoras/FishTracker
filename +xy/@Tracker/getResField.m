function rfield = getResField(self,res,field,delinvif)
% RFIELD = GETRESFIELD(SELF,RESmFIELD,DELINVIF)

  
  if nargin<4
    delinvif = 0;
  end
  
  if delinvif
    rfield = self.deleteInvisible(res,field);
  else
    if strcmp(field,'pos')
      rfield = res.pos;
    elseif isfield(res.tracks,field)
      rfield = double(res.tracks.(field));
    else
      error(sprintf('Cannot find requested field "%s"',upper(field)));
    end
  end
  

