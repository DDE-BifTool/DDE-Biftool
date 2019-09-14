function pfoldini=POfoldInit(funcs,point,method,ip,orig_free_par)
%% crude initial guess for fold of periodic orbits
%
% $Id: POfoldInit.m 369 2019-08-27 00:07:02Z jansieber $
%
pfoldini=point;
[f,x0]=p_correc_setup(funcs,point,orig_free_par,method.point,...
    'remesh_flag',0,'previous',point,'output','J');
J=f(x0);
%dde_psol_jac_res(funcs,pfoldini,orig_free_par,method.point);
nullvecs=dde_svdspaces_lr(J,1);
vpoint=dde_point_from_x(nullvecs,pfoldini,orig_free_par);
normv2=p_dot(vpoint,vpoint,'free_par_ind',orig_free_par,'period',true);
normv=sqrt(normv2);
vpoint=p_axpy(1/normv,vpoint,[]);
pfoldini.profile=[pfoldini.profile;vpoint.profile];
pfoldini.parameter(1:ip.nuserpar)=pfoldini.parameter;
pfoldini.parameter(ip.beta)=vpoint.period;
pfoldini.parameter(ip.period)=pfoldini.period;
if isfield(ip,'nullparind')
    pfoldini.parameter(ip.nullparind(:,2))=vpoint.parameter(ip.nullparind(:,1));
end
if ~funcs.tp_del
    pfoldini.parameter(ip.ext_tau)=point.parameter(funcs.sys_tau());
end
end