
function  save(self,savename,savepath,vname)
% SELF.SAVE([SAVENAME,SAVEPATH,VNAME]) saves the object

  [defname,exists] = self.getDefaultFileName();
  
  if ~exist('savename','var') || isempty(savename)
    [~,savename] = fileparts(defname);
  end
  if ~exist('savepath','var') || isempty(savepath)
    savepath = fileparts(savename);
    if isempty(savepath)
      [savepath] = fileparts(defname);
    end
  end
  [~,savename] = fileparts(savename);
  
  if ~exist('vname','var') || isempty(vname)
    vname = 'ft';
  end
  eval([vname '=self;']);

  if ~isempty(savepath)
    fname = [savepath filesep savename '.mat'];
  else
    fname = [savename '.mat'];
  end

  % do not overwrite
  fname1 = fname;
  s = 0;
  while exist(fname1,'file')
    s = s+1;
    [a,b,c] = fileparts(fname);
    fname1 = [a filesep b '-' num2str(s) c];
  end
  fname = fname1;
  
  save(fname,vname,'-v7.3');
  
  
  fish.helper.verbose('saved variable %s to %s',vname,fname);
end      
