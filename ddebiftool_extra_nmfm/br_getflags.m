function spoints = br_getflags(branch,kind)
%% Collect list of special point numbers in branch
%
% $Id: br_getflags.m 309 2018-10-28 19:02:42Z jansieber $
%
%%
%#ok<*AGROW>
uninitialized=arrayfun(@(x)~ischar(x.flag),branch.point);
if any(uninitialized)
    warning('br_getflags:ini',...
        'uninitialized flags at points %s',num2str(find(uninitialized)));
end
%% convert uninitialized flags to ''
if nargin>1
    spoints=find(arrayfun(@(x)strcmp(x.flag,kind),branch.point));
else
    branch.point=arrayfun(@(x)setfield(x,'flag',[x.flag,'']),branch.point); %#ok<SFLD>
    ll = length(branch.point);
    spoints = [];
    num = [];
    for i = 1:ll
        f = branch.point(i).flag;
        ind = bif2num(f);
        if ind > 0
            if ind>length(num)
                num(ind)=0;
            end
            num(ind) = num(ind) + 1;
            spoints(ind, num(ind)) = i;
        end
    end
end
end

