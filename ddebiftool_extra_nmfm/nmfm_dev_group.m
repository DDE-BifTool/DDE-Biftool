function groups=nmfm_dev_group(funs)
%% Group history functions into groups of equals
% output cell array of index sets where each cell contains indices with
% equal funs.
%
% $Id: nmfm_dev_group.m 309 2018-10-28 19:02:42Z jansieber $
%%
nf=length(funs);
remaining=1:nf;
groups=cell(1,nf);
count=0;
%% for test purposes (such that funs can be complex numbers or dev functions)
if isnumeric(funs)
    iseq=@(x,y)all(x(:)==y(:));
else
    iseq=@nmfm_dev_equal;
end
%% collect into groups
while ~isempty(remaining)
    count=count+1;
    groups{count}=remaining(1);
    remaining(1)=[];
    for i=length(remaining):-1:1
        if iseq(funs(groups{count}(1)),funs(remaining(i)))
            groups{count}(end+1)=remaining(i);
            remaining(i)=[];
        end
    end
end
groups=groups(1:count);
end
 