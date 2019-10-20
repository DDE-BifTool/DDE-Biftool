%% SetupPOfold - Initialize continuation of folds of periodic orbits
%%
function [pfuncs,pbranch,suc]=SetupPOfold(funcs,branch,ind,varargin)
%% Inputs
%
%  * |funcs|: functions used for DDE
%  * |branch|: branch of psols along which fold was discovered
%  * |ind|: index in points array that is closest to fold (for initial
%      guess)
%
% optional inputs
%
% * |contpar| (integer default |[]|): set of continuation parameters (if []
%   then free parameters of input branch are used)
% * |sys_deri| (default |1e-4|): used for finite differencing when approximating
%   jacobian of rhs, will be replaced by funcs.sys_deri if funcs.sys_deri is
%   provided by user
% * |correc| (logical, default true): apply |p_correc| to first points on fold
%   branch
% * |dir| (integer, default |[]|): which parameter to vary initially along fold
%   branch (pbranch has only single point if |dir| is empty)
% * |step| (real, default |1e-3|): size of initial step if |dir| is non-empty
% * |pitchfork| (logical, default false) if true then initialization assumes
%   that detected "fold" is pitchfork bifurcation (not safe)
% * |hjac| (default |1e-4|) deviation for numerical derivatives if needed
%
% All other named arguments are passed on to pbranch.method.continuation,
% pbranch.method.point and pbranch.parameter
%% Outputs
%
% * |pfuncs|: functions used for extended DDE
% * |pbranch|: fold branch with first point (or two points)
% * |suc|: flag whether corection was successful
%
%
% <html>
% $Id: SetupPOfold.m 374 2019-09-14 14:02:58Z jansieber $
% </html>
%
%% process options
default={'contpar',[],'sys_deri',1e-4,'sys_dtau',1e-4,'correc',true,'dir',[],...
    'step',1e-3,'nextrapar',0,'hjac',1e-4,'df_deriv',false,...
    'nullparind',zeros(0,1),'extra_cond',{}};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
branch.point=branch.point(ind); % remove all points but approx fold
%% If branch is psol bifurcation remove extra args
[funcs,branch]=PsolFromPsolbif(funcs,branch);
if isfield(funcs,'sys_deri_provided') && funcs.sys_deri_provided
    options.sys_deri=funcs.sys_deri;
end
if funcs.tp_del && isfield(funcs,'sys_dtau_provided') && funcs.sys_dtau_provided
    options.sys_dtau=funcs.sys_dtau;
end
%% initialize branch of folds (pbranch) and pass on optional args
pbranch=replace_branch_pars(branch,options.contpar,pass_on);
pbranch.method.point.extra_condition=1;
%% set up numbering and values of additional parameters and artificial delays
point=dde_psol_create(pbranch.point);
%% dimension of original problem
if isfield(point,'profile')
    ip.dim=size(point.profile,1);    
elseif isfield(point,'x')
    ip.dim=length(point.x);
end
ip.extdim=ip.dim;
ip.nuserpar=length(point.parameter); % number of original system parameters
%% extend problem
ip.beta=ip.nuserpar+1;       % location of add. parameter beta
ip.period=ip.nuserpar+2;     % location of copy of period
ip.nextrapar=options.nextrapar;
nnull=length(options.nullparind);
ip.nullparind=NaN(nnull,2); % location of parameters for nullspace vector
ip.nullparind(:,1)=options.nullparind;
ip.nullparind(:,2)=ip.period+(1:nnull)';
get_comp=@(p,component)extract_from_POfold(p,component,ip);
user_lhs_num=funcs.lhs_matrix(ip.dim);
lhs_num=kron(eye(2),user_lhs_num);
fun_args={'x_vectorized',funcs.x_vectorized,...
    'lhs_matrix',lhs_num};
if ~funcs.tp_del % constant delay
    ip.orig_tau=funcs.sys_tau(); % system delays
    % additional delays needed for extended system:
    ip.ext_tau=ip.period+nnull+(1:length(ip.orig_tau));
    ip.orig_ntau=length(ip.orig_tau);
    delay_cond={@(p,pref)sys_cond_extradelay(p,ip.orig_tau,ip.orig_tau,ip.ext_tau)};
    %% set up functions of extended system
    sys_rhs=@(x,p)sys_rhs_POfold(x,p,ip,funcs,options.sys_deri);
    sys_tau=@()[ip.orig_tau,ip.ext_tau];
    %sys_deri=@(x,p,nx,np,v)...
    %    sys_deri_POfold(x,p,nx,np,v,options.hjac,...
    %    beta_ind,period_ind,dim,xtau_ind,funcs,sys_rhs,options.df_deriv);
    pfuncs=set_funcs('sys_rhs',sys_rhs,'sys_tau',sys_tau,...%disabled forn now: 'sys_deri',sys_deri,...
        fun_args{:});
else % state-dependent delay
    % indices of additional delays needed for extended system:
    ip.orig_tau=[];
    ip.orig_ntau=funcs.sys_ntau();  % number of state-dependent delays
    ip.ext_tau=[];
    delay_cond={};
    %% set up functions of extended system
    sys_rhs=@(x,p)sys_rhs_SD_POfold(x,p,ip,funcs,options.sys_deri,options.sys_dtau);
    sys_tau=@(itau,x,p)sys_tau_SD_PObif(itau,x,p,ip,funcs);
    sys_ntau=@()ip.orig_ntau*2;
    pfuncs=set_funcs('sys_rhs',sys_rhs,'sys_tau',sys_tau,...
        'sys_ntau',sys_ntau,...
        fun_args{:},'hjac',options.hjac);
end
%% required amendments of structures for extended system
pfuncs.delayed_derivs=[zeros(1,ip.orig_ntau),ones(1,ip.orig_ntau)];
free_par=[pbranch.parameter.free,...
    ip.beta,ip.period,ip.nullparind(:,2)',ip.ext_tau];
[dum,is]=unique(free_par); %#ok<ASGLU>
pbranch.parameter.free=free_par(sort(is));
%% extended delays are for derivatives
pfuncs.sys_cond_reference=true;
pfuncs.get_comp=get_comp;
pfuncs.kind='POfold';
pfuncs.userfuncs=funcs;
pfuncs.ip=ip;
pfuncs.sys_cond=@(p,pref)dde_sys_cond_collect(pfuncs,p,pref,...
    {@(p,pref)sys_cond_POfold(p,ip),...
    @(p,pref)sys_cond_fixperiod(p,ip.period),...
    delay_cond{:},...
    options.extra_cond{:}}); %#ok<CCAT>
%% create initial guess for correction
pfoldini0=POfoldInit(funcs,point,branch.method,ip,ip.nullparind(:,1)');
%% for branch points, preprocess to add dummy parameters
if pfuncs.ip.nextrapar>0
    pbranch.method.point.preprocess='dde_BP_preprocess';
end
%% correct initial guess and find 2nd point along branch if desired
[pbranch,suc]=correct_ini(pfuncs,pbranch,pfoldini0,...
    options.dir,options.step,options.correc);
end
%
