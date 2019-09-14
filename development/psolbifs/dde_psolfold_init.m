function pfoldini=dde_psolfold_init(funcs,point,free_par,method)
%% crude initial guess for fold of periodic orbits
%
% $Id: dde_psolfold_init.m 362 2019-07-14 15:49:40Z jansieber $
%
pfoldini=dde_psolfold_create('point',point,'free_par',free_par);
J=dde_psol_jac_res(funcs,point,free_par,method.point);
nullvecs=dde_svdspaces_lr(J,1);
vpoint=dde_point_from_x(nullvecs,point,pfoldini.free_par);
normv2=dde_coll_profile_dot(vpoint,vpoint)+...
    vpoint.period^2/point.period^2+...
    sum(vpoint.parameter(pfoldini.free_par).^2);
normv=sqrt(normv2);
pfoldini.profile_var=vpoint.profile/normv;
pfoldini.parameter_var=[vpoint.period/pfoldini.period,...
    vpoint.parameter(pfoldini.free_par)]/normv;
if ~funcs.tp_del
    pfoldini.parameter_delay=[0,point.parameter(funcs.sys_tau())];
else
    pfoldini.parameter_delay=[];
end    
pfoldini.profile_var=vpoint.profile/normv;
end