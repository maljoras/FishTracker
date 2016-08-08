function [a varargout]=labelsubplot(figureh,labels,fs);
% labels all children of a given figure handle by putting a new axis with
% a label on top of it. Gives back the handles of the new axes.
    
  if ~exist('figureh','var') || isempty(figureh)
    figureh = gcf;
  end
  
  
  allch = get(figureh,'children');
  
  if isempty(allch)
    a =[];
    return;
  end
  
  
  if ~exist('labels','var') || isempty(labels)
    labels = 'ABCDEFGHIJKLMNOBQRSTUVWXYZ';
  end
  
  if ~exist('fs','var') || isempty(fs)
    fs = 16;
  end
  

  s = 0;
  for i = 1:length(allch)
    if ~strcmp(get(allch(i),'Type'),'axes') || ~isempty(get(allch(i),'Tag')), continue;end
    s =s+1;
    ch(s) = allch(i);
    p(s,:) = get(ch(s),'position');
  end

  %from lexyT to right, top to down
  [dummy,inds] = sortrows([-p(:,2),p(:,1)]);
  
  ch = ch(inds);
  s = 0;
  pold = [-1,-1];
  for i = 1:length(ch)
    p1 = get(ch(i),'position');
    if all(pold(1:2)==p1(1:2))
      continue
    end
    das = get(ch(i),'dataaspect');
    s = s+1;
    pold = p1;
    %delete(ch(i))
    a(s) =axes;
    set(a(s),'position',p1,'color','none');
    axis(a(s),'off');
    
    if das(1)==1 && das(2)==1
      xl = diff(get(ch(i),'xlim'));
      yl = diff(get(ch(i),'ylim'));
      daspect(a(s),[1,xl/yl,1])
    end
    
    if iscell(labels)
      label = labels{s};
    else
      label = labels(s);
    end
      
    axes(a(s))
    d(s) = text(-0.15,1,label,'fontweight','bold','fontsize', ...
                fs,'verticalalignment','bottom','horizontalal','right');
%    d(s) = text(-0.15,1,a(s),['\textsf{' labels(s) '}'],'fontweight','bold','fontsize', ...
%                fs,'verticalalignment','bottom','interpreter','latex');

   set(a(s),'Tag','LabelAxes')
  end
  
  
  %put new axes behind old
  set(figureh,'children',circshift(get(figureh,'children'),[-length(a),0]));

  varargout{1} = d;