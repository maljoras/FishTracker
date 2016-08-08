function outstr = allfieldnames(p)
% ALLFIELDNAMES(S) returns a cell array of all fieldnames from a struc with
% multiple layers in a form that eval can be used to set the fields. 
%
% MJR 12.5.2006
  
  s = 0;    
  outstr = [];
  
  if isstruct(p) && ~isempty(p) && length(p)==1
    
    f = fieldnames(p);

    for i = 1:length(f)

      daughternodes = xy.helper.allfieldnames(p.(f{i}));
      
      if ~isempty(daughternodes)
        for j = 1:length(daughternodes)
          s = s+1;  
          if strcmp(daughternodes{j}(1),'{') || strcmp(daughternodes{j}(1),'(')
            outstr{s} = [ f{i}  daughternodes{j}];
          else
            outstr{s} = [ f{i} '.' daughternodes{j}];
          end
        end
      else
        s = s+1;      
        outstr{s} = f{i};
      end
      
    end

  elseif iscell(p) && ~isempty(p)
    %same for cells

    for i = 1:length(p(:))
      if ~isstruct(p{i})
        continue
      end

      daughternodes = xy.helper.allfieldnames(p{i});
        
      if ~isempty(daughternodes)
        for j = 1:length(daughternodes)
          s = s+1;      
          outstr{s} = [ '{' num2str(i) '}.' daughternodes{j}];
        end
      else
        s = s+1;      
        outstr{s} = [];
      end
    
    end

  elseif length(p)>1 && ~isempty(p)
    %same for arrays

    for i = 1:length(p(:))
      if ~isstruct(p(i))
        continue
      end

      daughternodes = xy.helper.allfieldnames(p(i));
        
      if ~isempty(daughternodes)
        for j = 1:length(daughternodes)
          s = s+1;    

          if isvector(p)
            str = num2str(i);
          else
            c = {};
            [c{1:length(size(p))}] = ind2sub(size(p),i);
            str = [];
            for iii = 1:length(c)-1
              str = [str,num2str(c{iii}) ','];
            end
            str = [str, num2str(c{end})];
          end
          outstr{s} = [ '(' str ').' daughternodes{j}];
        end
      else
        s = s+1;      
        outstr{s} = [];
      end
    
    end

  else
    outstr = [];
  end
  

