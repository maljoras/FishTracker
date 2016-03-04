
function nfish = chooseNFish(fname,nfish)
  if hasOpenCV
    vid = FishVideoReader(fname);
  else
    vid = FishVideoReaderMATLAB(fname);
  end
  
  frame = readFrame(vid);
  
  f = figure;
  p1 = get(f,'position');
  %p1(4) = p1(3)/2;
  set(f,'name','Choose the number of fish:','WindowStyle','modal');

  a = subplot(1,2,1);
  image(frame);
  axis off
  handles = [];
  handles.cancel= 0;
  daspect(a,[1,1,1]);
  p1 = get(a,'position');
  p1(3) = p1(3)*2;
  p1(1) = 0.05;
  set(a,'pos',p1);

  
  title('Choose the number of fish:');
  
  gap = 0.02;
  bw = 0.07;
  p2 = [p1(1)+p1(3)+gap, p1(2), 1-p1(3)-p1(1)-2*gap, p1(4)];
  handles.sld = uicontrol('Style', 'listbox', 'Min',1,'Max',1,'Value',nfish,...
                          'units','normalized','Position', p2,...
                          'String',num2cell([1:50]),'Fontsize',18,...
                  'Callback', @subCallbackList); 

  handles.btn = uicontrol('Style', 'pushbutton', 'String', 'OK',...
        'units','normalized','Position',[p2(1),p2(2)-gap-bw,p2(3),bw],...
        'Callback', 'uiresume(gcbf)');       
  handles.btn2 = uicontrol('Style', 'pushbutton', 'String', 'Do Not Track',...
        'units','normalized','Position',[p2(1)-p2(3)-gap,p2(2)-gap-bw,p2(3),bw],...
        'Callback', @subCallbackCancel);       

  set(f,'UserData',handles);

  uiwait(f); 
  
  if ~ishandle(f)
    error('Window closed');
  end
  
  handles = get(f,'UserData');
  if handles.cancel
    nfish= [];
  else
    nfish = handles.sld.Value;
  end
  close(f);
  
  
end

function subCallbackList(s,c)
  f=gcbf;
  f.UserData.txt.String = num2str(floor(s.Value));
end

function subCallbackCancel(s,c)
  f=gcbf;
  f.UserData.cancel = 1;
  uiresume(f);
end
