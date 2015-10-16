function plotSegments(seg)

  
  if ~length(seg)
    return;
  end

  clf; 
  [r1,r2] = getsubplotnumber(length(seg));
  
  for i = 1:length(seg)
    
    subsubplot(r1,r2,i,2,1,1);
    x = 0:size(seg(i).FilledImage,2)-1;
    y = 0:size(seg(i).FilledImage,1)-1;
    imagesc(x,y,seg(i).FilledImage);
    daspect([1,1,1])
    axis off;
    hold on; 
    c = cat(1,seg(i).CenterLine{:});
    if ~isempty(c)
      scatter(c(:,1),c(:,2),100,'w','filled');
    end
    
    subsubplot(r1,r2,i,2,1,2);
    imagesc(seg(i).RotImage);
    daspect([1,1,1])
    axis off;
  end
  
    