%% SetupRWFold - Initialize continuation of folds of relative equilibria
%%
function [pfuncs,pbranch,suc]=SetupRWFold(funcs,branch,ind,varargin)
%% Inputs
%
% * |funcs|: functions used for DDE
% * |branch|: branch of psols along which fold was discovered
% * |ind|: number of point close to fold
%
%% Outputs
%
% * |pfuncs|: functions used for extended DDE
% * |pbranch|: fold branch with first point (or two points)
% * |suc|: flag whether corection was successful
%
%% Optional inputs
%
% * |contpar| (integer default |[]|): set of continuation parameters  
%   (if empty free parameters in argument branch are used)
% * |hbif| (default |1e-3|): used for finite differencing when approximating
%   linearized system,
% * |correc| (logical, default |true|): apply |p_correc| to first points on fold
%   branch
% * |dir| (integer, default |[]|): which parameter to vary initially along fold
%   branch (|pbranch| has only single point if dir is empty)
% * |step| (real, default |1e-3|): size of initial step if dir is non-empty
% * |hjac| (default |1e-6|) deviation for numerical derivatives if needed
%
% all other named arguments are passed on to pbranch.method.continuation,
% pbranch.method.point and pbranch.parameter
%
% $Id: SetupRWFold.m 357 2019-06-30 23:59:31Z jansieber $
%

%% process options
default={'contpar',[],'hbif',1e-3,'correc',true,'dir',[],...
    'step',1e-3,'hjac',1e-6};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
branch.point=branch.point(ind); % remove all points but approx fold
% initialize branch of folds (pbranch)
pbranch=replace_branch_pars(branch,options.contpar,pass_on);
%% set up numbering and values of additional parameters and artificial delays
point=pbranch.point;
if isfield(point,'stability')
    point=rmfield(point,'stability');
end
dim=size(point.x,1);            % dimension of original problem
ind_omega=length(point.parameter);   % number of original system parameters
ind_rho=ind_omega+1;                 % location of extra parameter rho
%% set up functions of extended system
pfuncs=funcs;
pfuncs.userfuncs=funcs;
pfuncs.get_comp=@(p,component)extract_from_RWfold(p,component,ind_omega);
pfuncs.sys_rhs=@(x,p)sys_rhs_RWFold(x,p,ind_omega,ind_rho,funcs.sys_rhs,dim,funcs.sys_dirderi{1});
pfuncs.orig_cond=funcs.sys_cond;
pfuncs.orig_cond_reference=funcs.sys_cond_reference;
pfuncs.sys_cond=@(p,pref)sys_cond_RWFold(p,pfuncs,dim,ind_rho,pref);
pfuncs=pfuncs.add_deriv(pfuncs,@(x,p,dx,dp)pfuncs.sys_rhs(x,p),0);
pfuncs.kind='RWFold';
%% required amendments of structures for extended system
pbranch.parameter.free=[pbranch.parameter.free,ind_rho];
pbranch.method.point.extra_condition=1;
%% create initial guess for correction
pfoldini0=RWfoldInit(funcs,point,branch.parameter.free,pbranch.method.point);
%% correct initial guess and find 2nd point along branch if desired
[pbranch,suc]=correct_ini(pfuncs,pbranch,pfoldini0,...
    options.dir,options.step,options.correc);
end

%% crude initial guess
function pfoldini=RWfoldInit(funcs,point,free_par_ind,mth)
pfoldini=point;
J=p_correc_rhs(funcs,mth,point,free_par_ind,'pref',point,'output','J');
%J=stst_jac(funcs,point.x,point.parameter,free_par_ind);
%[rdum,Jcond]=funcs.sys_cond(point); %#ok<ASGLU>
%J=[J;Jcond.x',Jcond.parameter(free_par_ind)];
[U,S,V]=svd(J); %#ok<ASGLU>
nullvecs=V(:,end);
v=nullvecs(1:length(point.x));
rho=nullvecs(end);
vpoint=p_axpy(0,point,[]);
vpoint.x=v;
vpoint.parameter=rho;
normv=sqrt(v'*v+rho^2);
rho=rho/normv;
vpoint.x=vpoint.x/normv;
pfoldini.x=[pfoldini.x;vpoint.x];
pfoldini.parameter=[pfoldini.parameter,rho];
end
%% extract components from pfold
function result_array=extract_from_RWfold(pfold_array,component,npar)
dim=size(pfold_array(1).x,1)/2;
for i=1:length(pfold_array)
    pfold=pfold_array(i);
    switch component
        case 'kind'
            result='RWfold';
        case {'solution','solution_for_stability'}
            result=pfold;
            result.x=result.x(1:dim,:);
            result.parameter=result.parameter(1:npar);
        case 'nullvector'
            result=pfold;
            result.x=result.x(dim+1:end,:);
            result.parameter=result.parameter(npar+1);
    end
    result_array(i)=result; %#ok<AGROW>
end
end
