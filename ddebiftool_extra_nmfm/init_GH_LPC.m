function [LPCfuncs,LPCbranch,psol1,psol2] = init_GH_LPC(funcs,GH,varargin)
%% Initialize branch for continuing the Limit point of cycles curve
% emanating from the the generalized-Hopf point.
%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: init_GH_LPC.m 309 2018-10-28 19:02:42Z jansieber $
%%
%% process options
default={'contpar',[],'sys_deri',1e-4,'sys_dtau',1e-4,'correc',true,'dir',[],...
    'step',1e-3,'pitchfork',false,'hjac',1e-4,'df_deriv',false};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
if isfield(funcs,'sys_deri_provided') && funcs.sys_deri_provided
    options.sys_deri=funcs.sys_deri;
end
if funcs.tp_del && isfield(funcs,'sys_dtau_provided') && funcs.sys_dtau_provided
    options.sys_dtau=funcs.sys_dtau;
end
%% approximate first two cycles
psol1=cycle_approx(GH,options.radius);
psol2=cycle_approx(GH,options.radius*1.01);

%% Initialize Limit point of cycles branch
branch=df_brnch(funcs,get_free_pars(),'psol');
branch.point=psol1; % add first cycle
%% initialize branch of folds (pbranch) and pass on optional args
pbranch=replace_branch_pars(branch,options.contpar,pass_on);
%% set up numbering and values of additional parameters and artificial delays
if isfield(psol1,'stability')
    psol1=rmfield(psol1,'stability');
end
dim=size(psol1.profile,1);    % dimension of original problem
npar=length(psol1.parameter); % number of original system parameters
beta_ind=npar+1;              % location of add. parameter beta
period_ind=npar+2;            % location of copy of period
pfuncs=funcs;
pfuncs.get_comp=@(p,component)extract_from_POfold(p,component,npar);

orig_tau_ind=funcs.sys_tau(); % system delays
% additional delays needed for extended system:
[ext_tau,xtau_ind,relations]=generate_tau_ext(psol1.parameter,orig_tau_ind,...
    'ext_tau_base',period_ind);
ext_tau_ind=period_ind+(1:length(ext_tau));
%% set up functions of extended system
pfuncs.sys_rhs=@(x,p)sys_rhs_POfold(x,p(1:npar),p(beta_ind),p(period_ind),...
    [0,p(orig_tau_ind)],funcs.sys_rhs,dim,xtau_ind,options.sys_deri);
pfuncs.sys_tau=@()[orig_tau_ind,ext_tau_ind];
pfuncs.sys_deri=@(x,p,nx,np,v)...
    sys_deri_POfold(x,p,nx,np,v,options.hjac,...
    beta_ind,period_ind,dim,xtau_ind,funcs,pfuncs,options.df_deriv);
%% required amendments of structures for extended system
pbranch.parameter.free=[pbranch.parameter.free,beta_ind,period_ind,ext_tau_ind];
%%
pfuncs.sys_cond=@(p)sys_cond_POfold(p,funcs.sys_cond,dim,beta_ind,...
    period_ind,pfuncs.get_comp,relations);
%% create initial guess for correction
pfoldini1=POfoldInit(funcs,psol1,branch.method.point,ext_tau,options.pitchfork);
%% correct initial guess and find 2nd point along branch if desired
pbranch.method.point.extra_condition=1;
pfoldini2=POfoldInit(funcs,psol2,branch.method.point,ext_tau,options.pitchfork);

pbranch.point=pfoldini1;
pbranch.point(2)=pfoldini2;

LPCfuncs=pfuncs;
LPCbranch=pbranch;

end

function psol = cycle_approx(GH,eps)
%% normal form coefficients
PHI=GH.nmfm.PHI;
H1001=GH.nmfm.H1001;
K01=GH.nmfm.K01;
H0001=GH.nmfm.H0001;
H11=GH.nmfm.H11;
H20=GH.nmfm.H20;
H30=GH.nmfm.H30;
H21=GH.nmfm.H21;
L2=GH.nmfm.L2;
c1=GH.nmfm.c1;
c2=GH.nmfm.c2;
b12=GH.nmfm.b12;

%% approxiamtion to the cycle
T=2*pi/(GH.omega+(imag(c1)-2*real(c2)*b12)*eps^2);
nphase=1;n=length(GH.x);
x0 = zeros(n,80*nphase+1);
for i=1:80*nphase+1
    zz = -exp(sqrt(-1.0)*(2*pi*(i-1)/(80*nphase)+pi/2));
    x0(:,i) = GH.x+2*real(zz*PHI(:,1))*eps+(H11(:,1)-2*L2*H0001(:,1)+real(zz^2*H20(:,1)))*eps^2+...
        (-4*L2*real(zz*H1001(:,1))+1/3*real(zz^3*H30(:,1))+real(zz*H21(:,1)))*eps^3;
end
KK=GH.parameter(get_free_pars())'-2*L2*K01*eps^2;

%% setup periodic orbit
psol.kind='psol';
psol.parameter=GH.parameter;
psol.parameter(get_free_pars())=real(KK)';
psol.mesh=linspace(0,1,80*nphase+1);
psol.degree=4;
psol.profile=real(x0);
psol.period=real(T);
end
