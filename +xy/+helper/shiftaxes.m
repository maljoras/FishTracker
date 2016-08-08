function shiftaxes(handle,shift)
% SHIXYTAXES(HANDLE,SHIXYT) shifts axes by SHIXYT amount

  
  for i = 1:length(handle)
    p1 = get(handle(i),'position');
    p1(1:length(shift)) = p1(1:length(shift))+shift;
    set(handle(i),'position',p1);
  end
  
  