function data=dde_psol2torus(funcs,info,varargin)
%% constructor for torus bifurcation from psol
% is equal to SetupTorusBifurcation but without correction and with square
% Jacobian
% $Id: dde_psol2torus.m 344 2019-05-12 00:08:16Z jansieber $
%% 
default={'hjac',1e-4,'nremove',1,'biftype','torus','closest',[]};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
point=info.point(1);
%% initialize branch of torus bifurcations (trbranch) and pass on optional args
info=replace_branch_pars(info,[],pass_on);
%% set up numbering and values of additional parameters
point=dde_psol_create(point);
dim=size(point.profile,1);    % dimension of original problem
npar=length(point.parameter); % number of original system parameters
omega_ind=npar+1;             % location of add. parameter omega
period_ind=npar+2;            % location of add. parameter (equal to period)
%% set up functions of extended system
get_comp=@(p,component)extract_from_tr(p,component,options.biftype);
if ~funcs.tp_del 
    %% constant delays
    tau_ind=funcs.sys_tau();
    sys_rhs=@(x,p)sys_rhs_TorusBif(x,p(1:npar),...
        p(omega_ind),p(period_ind),[0,p(tau_ind)],...
        funcs.sys_rhs,dim,funcs.sys_deri);
    sys_deri=@(x,p,nx,np,v)...
        sys_deri_TorusBif(x,p,nx,np,v,options.hjac,omega_ind,period_ind,dim,...
        funcs,struct('sys_rhs',sys_rhs));
    trfuncs=set_funcs('sys_rhs',sys_rhs,'sys_tau',@()tau_ind,...
        'sys_deri',sys_deri,'x_vectorized',funcs.x_vectorized);
else
    %% state dependent delay
    n_orig_tau=funcs.sys_ntau();  % number of state-dependent delays
    % additional delays needed for extended system:
    xtau_ind=tauSD_ext_ind(n_orig_tau);
    sys_rhs=@(x,p)sys_rhs_SD_TorusBif(x,p(1:npar),...
        p(omega_ind),p(period_ind),...
        funcs.sys_rhs,funcs.sys_tau,funcs.sys_deri,funcs.sys_dtau,dim,xtau_ind);
    sys_ntau=@()n_orig_tau*(n_orig_tau+1);
    sys_tau=@(itau,x,p)sys_tau_SD_PObif(itau,x,p(1:npar),funcs.sys_tau,dim,xtau_ind);
    %sys_dtau=@(itau,x,p,nx,np)sys_dtau_SD_PObif(itau,x,p,nx,np,funcs.sys_dtau,...
    %    dim,xtau_ind);
    trfuncs=set_funcs('sys_rhs',sys_rhs,'sys_tau',sys_tau,'sys_ntau',sys_ntau,...
        'x_vectorized',funcs.x_vectorized,'hjac',options.hjac);
end
trfuncs.sys_cond=@(p)sys_cond_TorusBif(p,funcs.sys_cond,dim,period_ind,...
    get_comp,false);
trfuncs.get_comp=get_comp;
trfuncs.kind=options.biftype;
trfuncs.userfuncs=funcs;
%% required amendments of structures for extended system
info.parameter.free=[info.parameter.free,omega_ind,period_ind];
info.method.point.extra_condition=1;
%% create initial guess for correction
info.point=TorusInit(funcs,point,info.method,options.nremove,options.closest);
data.info=info;
data.funcs=trfuncs;
end
