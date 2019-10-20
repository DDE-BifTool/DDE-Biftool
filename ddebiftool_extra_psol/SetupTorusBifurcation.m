%% SetupTorusBifurcation - Initialize continuation of torus or period doubling bifurcations
%%
function [trfuncs,trbranch,suc]=SetupTorusBifurcation(funcs,branch,ind,varargin)
%% Inputs
%
% * |funcs|: functions used for DDE
% * |branch|: branch of psols along which bifurcation was discovered
%
% optional inputs
%
% * |contpar| (integer default |[]|):  set of continuation parameters (if []
%   then free parameters of input branch are used)
% * |sys_deri| (default |1e-4|): used for finite differencing when approximating
%   jacobian of rhs, will be replaced by funcs.sys_deri if funcs.sys_deri is
%   provided by user
% * |correc| (logical, default true): apply |p_correc| to first points on
%    torus branch
% * |dir| (integer, default |[]|): which parameter to vary initially along
%    torus branch (trbranch has only single point if |dir| is empty)
% * |step| (real, default |1e-3|): size of initial step if |dir| is non-empty
% * |hjac| (default |1e-4|) deviation for numerical derivatives if needed
%
% all other named arguments are passed on to trbranch.method.continuation,
% trbranch.method.point and trbranch.parameter
%% Outputs
%
% * |pfuncs|: functions used for extended DDE
% * |trbranch|: bifurcation branch with first point (or two points)
% * |suc|: flag whether corection was successful
%
%
% <html>
% $Id: SetupTorusBifurcation.m 374 2019-09-14 14:02:58Z jansieber $
% </html>
%
%% process options
default={'contpar',[],'sys_deri',1e-4,'sys_dtau',1e-4,'correc',true,'dir',[],'step',1e-3,...
    'hjac',1e-4,'nremove',1,'biftype','torus','closest',[],...
    'jacobian_nonsquare',true,'extra_cond',{},...
    'stop_1_2',false};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
branch.point=branch.point(ind);
%% If branch is psol bifurcation remove extra args
[funcs,branch,addremove]=PsolFromPsolbif(funcs,branch);
options.nremove=[options.nremove(:).',addremove(:).'];
%% initialize branch of torus bifurcations (trbranch) and pass on optional args
trbranch=replace_branch_pars(branch,options.contpar,pass_on);
trbranch.method.point.preprocess='dde_jac2square_preprocess';
trbranch.method.point.postprocess='dde_TorusBifurcation_postprocess';
%% set up numbering and values of additional parameters
point=dde_psol_create(trbranch.point);
ip.dim=size(point.profile,1);    % dimension of original problem
ip.nuserpar=length(point.parameter); % number of original system parameters
ip.omega=ip.nuserpar+1;             % location of add. parameter omega
ip.period=ip.omega+1;            % location of add. parameter (equal to period)
if isfield(funcs,'sys_deri_provided') && funcs.sys_deri_provided
    options.sys_deri=funcs.sys_deri;
end
if funcs.tp_del && isfield(funcs,'sys_dtau_provided') && funcs.sys_dtau_provided
    options.sys_dtau=funcs.sys_dtau;
end
%% set up functions of extended system
get_comp=@(p,component)extract_from_tr(p,component,options.biftype);
user_lhs_num=funcs.lhs_matrix(ip.dim);
lhs_num=kron(eye(3),user_lhs_num);
fun_args={'x_vectorized',funcs.x_vectorized,...
    'lhs_matrix',lhs_num};
if ~funcs.tp_del 
    %% constant delays
    ip.orig_tau=funcs.sys_tau();
    ip.orig_ntau=length(ip.orig_tau);
    sys_rhs=@(x,p)sys_rhs_TorusBif(x,p,ip,funcs,options.sys_deri,user_lhs_num);
    sys_deri=@(x,p,nx,np,v)...
        sys_deri_TorusBif(x,p,nx,np,v,options.hjac,ip,funcs,struct('sys_rhs',sys_rhs));
    trfuncs=set_funcs('sys_rhs',sys_rhs,'sys_tau',@()ip.orig_tau,...
        'sys_deri',sys_deri,fun_args{:});
else
    %% state dependent delay
    ip.orig_ntau=funcs.sys_ntau();  % number of state-dependent delays
    % additional delays needed for extended system:
    sys_rhs=@(x,p)sys_rhs_SD_TorusBif(x,p,ip,funcs,options.sys_deri,options.sys_dtau);
    sys_ntau=@()ip.orig_ntau*2;
    sys_tau=@(itau,x,p)sys_tau_SD_PObif(itau,x,p,ip,funcs);
    %sys_dtau=@(itau,x,p,nx,np)sys_dtau_SD_PObif(itau,x,p,nx,np,funcs.sys_dtau,...
    %    dim,xtau_ind);
    trfuncs=set_funcs('sys_rhs',sys_rhs,'sys_tau',sys_tau,'sys_ntau',sys_ntau,...
        fun_args{:},'hjac',options.hjac);
    trfuncs.delayed_derivs=[zeros(1,ip.orig_ntau),ones(1,ip.orig_ntau)];
end
trfuncs.sys_cond_reference=true;
%if options.jacobian_nonsquare
%    torus_extra_cond={};
%else
%    torus_extra_cond={@(p,pref)sys_cond_TorusBif(p,pref,ip)};
%end
trfuncs.get_comp=get_comp;
trfuncs.kind=options.biftype;
trfuncs.userfuncs=funcs;
trfuncs.ip=ip;
trfuncs.sys_cond=@(p,pref)dde_sys_cond_collect(trfuncs,p,pref,...
    {... %@(p,pref)sys_cond_TorusBif(p,pref,ip,options.jacobian_nonsquare),...
    @(p,pref)sys_cond_fixperiod(p,ip.period),...
    options.extra_cond{:}}); %#ok<CCAT>
%% required amendments of structures for extended system
free_par=[trbranch.parameter.free,ip.omega,ip.period];
[dum,is]=unique(free_par); %#ok<ASGLU>
trbranch.parameter.free=free_par(sort(is));
trbranch.method.point.extra_condition=1;
%% create initial guess for correction
trini0=TorusInit(funcs,point,trbranch.method,options.nremove,options.closest);
[trbranch,suc]=correct_ini(trfuncs,trbranch,trini0,...
    options.dir,options.step,options.correc);
if strcmp(options.biftype,'torus') && options.stop_1_2
    trbranch.parameter.max_bound=cat(1,trbranch.parameter.max_bound,...
        [ip.omega,1]);
    trbranch.parameter.min_bound=cat(1,trbranch.parameter.min_bound,...
        [ip.omega,-1]);
end
end
