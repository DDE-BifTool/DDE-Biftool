function [trfuncs,trbranch] = init_NS_br(funcs,psol1,psol2,varargin)
%% auxiliary function used to setup up Neimark-Sacker branch
% used by the functions init_ZH, init_HT and init_HOHO
%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: init_NS_br.m 314 2019-01-24 14:28:23Z mmbosschaert $
%%

%% process options
default={'contpar',[],'sys_deri',1e-4,'correc',true,'dir',[],'step',1e-3,...
    'hjac',1e-4,'nremove',1};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
if isfield(funcs,'sys_deri_provided') && funcs.sys_deri_provided
    options.sys_deri=funcs.sys_deri;
end
if funcs.tp_del && isfield(funcs,'sys_dtau_provided') && funcs.sys_dtau_provided
    options.sys_dtau=funcs.sys_dtau;
end
branch=df_brnch(funcs,get_free_pars(),'psol');
branch.point=psol1; % add first cycle
%% initialize branch of torus bifurcations (trbranch) and pass on optional args
trbranch=replace_branch_pars(branch,options.contpar,pass_on);
%% set up numbering and values of additional parameters
point=trbranch.point;
if isfield(point,'stability')
    point=rmfield(point,'stability');
end
dim=size(point.profile,1);    % dimension of original problem
npar=length(point.parameter); % number of original system parameters
omega_ind=npar+1;             % location of add. parameter omega
period_ind=npar+2;            % location of add. parameter (equal to period)
%% set up functions of extended system
trfuncs=funcs;
trfuncs.get_comp=@(p,component)extract_from_tr(p,component);

%% constant delays
tau_ind=funcs.sys_tau();
trfuncs.sys_rhs=@(x,p)sys_rhs_TorusBif(x,p(1:npar),...
    p(omega_ind),p(period_ind),[0,p(tau_ind)],...
    funcs.sys_rhs,dim,options.sys_deri);
trfuncs.sys_deri=@(x,p,nx,np,v)...
    sys_deri_TorusBif(x,p,nx,np,v,options.hjac,omega_ind,period_ind,dim,...
    funcs,struct('sys_rhs',trfuncs.sys_rhs));

trfuncs.sys_cond=@(p)sys_cond_TorusBif(p,funcs.sys_cond,dim,period_ind,...
    trfuncs.get_comp);
%% required amendments of structures for extended system
trbranch.parameter.free=[trbranch.parameter.free,omega_ind,period_ind];
trbranch.method.point.extra_condition=1;
%% create initial guess for correction
trini1=TorusInit(funcs,psol1,branch.method,options.nremove,[]);
trini2=TorusInit(funcs,psol2,branch.method,options.nremove,[]);

trbranch.point=trini1;
trbranch.point(2)=trini2;
end
