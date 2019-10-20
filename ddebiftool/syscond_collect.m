function [r,J]=syscond_collect(p,pref,conds,funcs)
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
if ~isfield(funcs,'get_comp')
    return
else
    r0=r;
    J0=J;
end
userpoint=get_comp(point,'solution');
if isfield(userfuncs,'sys_cond_reference') && userfuncs.sys_cond_reference
    userref=get_comp(pref,'solution');
    [r,userJ]=userfuncs.sys_cond(userpoint,userref);
else
    [r,userJ]=userfuncs.sys_cond(userpoint);
end
nuserpar=length(userpoint.parameter);
Jtemplate=p_axpy(0,point,[]);
J=repmat(Jtemplate,length(userJ),1);
%% append artificial parameters and components as zero in derivative
for i=1:length(userJ)
    J(i).parameter(1:nuserpar)=userJ(i).parameter;
    J(i).profile(1:ind.dim,:)=userJ(i).profile;
    J(i).period=userJ(i).period;
end
r=[r0;r];
J=[J0;J];
end
