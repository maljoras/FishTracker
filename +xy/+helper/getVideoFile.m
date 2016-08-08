    function vid = getVideoFile(fname)

      if exist('fname','var') && ~isempty(fname)
        [~,a,b] = fileparts(fname);
        fname = [a,b];
      else
        fname = '';
      end

      [filename, pathname]  = uigetfile({'*.avi';'*.mp4';'*.mpg';'*.*'}, ...
                                        sprintf('Pick movie file %s for tracking',fname));
      vid = fullfile(pathname,filename);
      if ~exist(vid)
        error('Please select file');
      end
      xy.helper.verbose('Selected file %s',vid);
    end
