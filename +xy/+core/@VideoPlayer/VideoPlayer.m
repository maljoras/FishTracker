classdef VideoPlayer < handle;
  
  
  
  properties 
    fig = [];
    figno  = 100;
  end
  

  
    
  methods
  
    
    function self = VideoPlayer() 
    % constructor
      self = self@handle();
      
      open(self);
    end
    
    function release(self)
      close(self.fig);
    end
    
    function bool =isOpen(self)
      bool = ishandle(self.fig);
    end

    function open(self)
      self.fig = figure(self.figno);
      drawnow;
    end
    
    function show(self)
      if ~isOpen(self)
        open(self);
      end
    end
    
    
    function step(self,frame)
      self.show();
      
      figure(self.fig)
      clf;
      image(frame);
      axis off;
      drawnow;
    end
    
    
  end
  
end
