function overlap = getBBoxOverlap(bbtracks,bbsegs,margin)

  overlap = zeros(size(bbtracks,1),size(bbsegs,1));
  for i = 1:size(bbtracks,1)
    
    % better judge overlap from bounding box 
    bbtrack = bbtracks(i,:);
    bbtrack(1:2) = bbtrack(1:2) - margin/2;
    bbtrack(3:4) = bbtrack(3:4) + margin;
    
    mnx = min(bbsegs(:,3),bbtrack(3));
    mny = min(bbsegs(:,4),bbtrack(4));
    
    bBoxXOverlap = min(min(max(bbtrack(1) + bbtrack(3) - bbsegs(:,1),0), ...
                           max(bbsegs(:,1) + bbsegs(:,3) - bbtrack(1),0)),mnx);
    bBoxYOverlap = min(min(max(bbtrack(2) + bbtrack(4) - bbsegs(:,2),0), ...
                           max(bbsegs(:,2) + bbsegs(:,4) - bbtrack(2),0)),mny);
    
    overlap(i,:) = bBoxXOverlap.*bBoxYOverlap./mnx./mny;
  end
end
