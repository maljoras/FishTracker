
function setProps(self,opts)

  if nargin<2
    opts = self.opts;
  end

  for f = fieldnames(self.opts)'
    if ~any(strcmp(f{1},{'tracks','classifier','blob','detector','reader','stimulus'}))  && isprop(self,f{1})
      self.(f{1}) = opts.(f{1});
    end
  end
end
