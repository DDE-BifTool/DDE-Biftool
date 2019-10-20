function [r,J]=dde_sys_cond_collect(funcs,p,pref,conds)
%% collect a cell array of extra conditions
% $Id$
%%
for i=length(conds):-1:1
    [res{i},Jac{i}]=conds{i}(p,pref);
    res{i}=res{i}(:)';
    Jac{i}=Jac{i}(:)';
end
[ru,Ju]=sys_cond_user(p,pref,funcs);
r=[ru;res{:}]';
J=[Ju;Jac{:}]';
end
function [r,J]=sys_cond_user(point,pref,funcs)
%% call user defined extra conditions and embed derivative into extended point
% for psol extensions only at the moment
% $Id$
%%
if isfield(funcs,'sys_cond_reference') && funcs.sys_cond_reference
    [r,J]=funcs.sys_cond(point,pref);
else
    [r,J]=funcs.sys_cond(point);
end
%% check if this is extended system (works for coll/psol)
if ~isfield(funcs,'get_comp') ||~isfield(funcs,'userfuncs')
    return
else
    r0=r;
    J0=J;
end
userpoint=funcs.get_comp(point,'solution');
userfuncs=funcs.userfuncs;
if isfield(userfuncs,'sys_cond_reference') && userfuncs.sys_cond_reference
    userref=funcs.get_comp(pref,'solution');
    [r,userJ]=userfuncs.sys_cond(userpoint,userref);
else
    [r,userJ]=userfuncs.sys_cond(userpoint);
end
J=dde_apply({'dde_',point.kind,'_extendblanks'},userJ,point);
r=[r0;r];
J=[J0;J];
end
