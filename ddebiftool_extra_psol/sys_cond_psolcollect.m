function [r,J]=sys_cond_psolcollect(p,pref,conds)
%% collect a cell array of extra conditions
% $Id: sys_cond_psolcollect.m 369 2019-08-27 00:07:02Z jansieber $
%%
for i=length(conds):-1:1
    [res{i},Jac{i}]=conds{i}(p,pref);
    res{i}=res{i}(:)';
    Jac{i}=Jac{i}(:)';
end
r=[res{:}]';
J=[Jac{:}]';
end
