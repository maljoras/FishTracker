function res = selectedFishIds(self,res,fishIds)
%RES = SELECTEDFISHIDS(SELF,RES,FISHIDS)

  for f = fieldnames(res.tracks)'
    S = [];
    S.type = '()';
    if size(res.tracks.(f{1}),2)==self.nfish
      sz= size(res.tracks.(f{1}));      
      S.subs = cell(1,length(sz));
      [S.subs{:}] = deal(':');
      S.subs{2} = fishIds;
      res.tracks.(f{1}) = subsref(res.tracks.(f{1}),S);
    end
  end
  
  if isfield(res,'pos')
    res.pos = res.pos(:,:,fishIds);
  end

end
