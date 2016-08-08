function res = selectedIdentityIds(self,res,identityIds)
%RES = SELECTEDIDENTITYIDS(SELF,RES,IDENTITYIDS)

  for f = fieldnames(res.tracks)'
    S = [];
    S.type = '()';
    if size(res.tracks.(f{1}),2)==self.nanimals
      sz= size(res.tracks.(f{1}));      
      S.subs = cell(1,length(sz));
      [S.subs{:}] = deal(':');
      S.subs{2} = identityIds;
      res.tracks.(f{1}) = subsref(res.tracks.(f{1}),S);
    end
  end
  
  if isfield(res,'pos')
    res.pos = res.pos(:,:,identityIds);
  end

end
