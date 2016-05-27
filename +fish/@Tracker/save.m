
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
  fname = [savepath filesep savename '.mat'];
  save([savepath filesep savename],vname,'-v7.3');

  fish.helper.verbose('saved variable %s to %s',vname,fname);
end      
