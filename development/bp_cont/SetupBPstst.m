%% SetupBPstst - Initialize continuation of folds of relative equilibria
%%
function [pfuncs,pbranch,suc]=SetupBPstst(funcs,branch,ind,cond,varargin)
%% Inputs
%
% * |funcs|: functions used for DDE
% * |branch|: branch of psols along which fold was discovered
% * |ind|: number of point close to fold
% * |cond|: extra condition(s) enforcing symmetry or invariance
%
%% Outputs
%
% * |pfuncs|: functions used for extended DDE
% * |pbranch|: fold branch with first point (or two points)
% * |suc|: flag whether corection was successful
%
%% Optional inputs
%
% * |invariance_cond_reference| (logical default |false|): if symmetry
% condition requires second argument (previous point)
%
% all other named arguments are passed on to SetupFold
%
% $Id$
%

%% process options
default={'contpar',[],'invariance_cond_reference',false,...
    'branchpoint',true};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
branch.point=branch.point(ind); % remove all points but approx fold
% initialize branch of folds (pbranch)
pbranch=replace_branch_pars(branch,options.contpar,pass_on);
%% set up numbering and values of additional parameters and artificial delays
point=pbranch.point;
if isfield(point,'stability')
    point=rmfield(point,'stability');
end
if isnumeric(cond) && ismatrix(cond)
    extra_cond=@(p,pref)sym_cond(p,cond);
elseif options.invariance_cond_reference
    extra_cond=cond;
else
    extra_cond=@(p,pref)cond(p);
end
ip.dim=size(point.x,1);            % dimension of original problem
ip.nuserpar=length(point.parameter);   % number of original system parameters
ip.nextrapar=length(extra_cond(point,point));
%% set up functions of extended system
pfuncs=funcs;
pfuncs.userfuncs=funcs;
pfuncs.ip=ip;
pfuncs.sys_cond_reference=true;
pfuncs.sys_cond=@(p,pref)dde_sys_cond_collect(pfuncs,p,pref,{extra_cond});
%% required amendments of structures for extended system
pbranch.method.point.extra_condition=1;
pbranch.method.point.preprocess='dde_BP_preprocess';
%% create initial guess for correction
[pbranch,suc]=SetupStstBifurcation(pfuncs,pbranch,1,'fold',...
    'branchpoint',false,'contpar',options.contpar,pass_on{:});
end

%%
function [r,J]=sym_cond(p,M)
ncond=size(M,1);
J=repmat(p_axpy(0,p,[]),ncond,1);
r=M*p.x;
for i=length(J):-1:1
    J.x=M(i,:).';
end
end