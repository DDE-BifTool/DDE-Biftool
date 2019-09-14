%% Initialize continuation of Fold bifurcations of periodic orbits
function [pfbranch,suc]=SetupPsolfold(funcs,branch,ind,varargin)
%% Inputs
%
% * |funcs|: functions used for DDE
% * |branch|: branch of stst along which Fold was discovered
% * |ind|: index of approximate Fold point
%
% Important optional inputs (name-value pairs)
% outputs
% foldbranch: Fold branch with first point (or two points)
% suc: flag whether corection was successful
%
% Important name-value pair inputs
%
% * |'contpar'| (integers default |[]|): index of continuation parameters
% * |'basepar'| (integers default |[]|): index of square problem
% parameters: if the basic problem has extra conditions and a corresponding
% set of free parameters then these need to be duplicated for the fold
% problem.
% * |'correc'| (logical, default true): apply |p_correc| to first points on fold
%   branch
% * |'dir'| (integer, default |[]|): which parameter to vary initially along fold
%   branch (foldbranch has only single point if dir is empty)
% * |'step'| (real, default |1e-3|): size of initial step if dir is non-empty
%
% All other named arguments are passed on to fields of |foldbranch|.
%% Outputs
% 
% * |foldbranch|: Fold branch with first point (or two points)
% * |suc|: flag whether corection was successful
%
% Parameter limits for etc are inherited from branch, unless overridden by
% optional input arguments.
%
% $Id: SetupPsolfold.m 375 2019-09-14 14:26:50Z jansieber $
%
%% process options
default={'contpar',[],'basepar',[],'correc',true,'dir',[],'step',1e-3};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
% initialize branch of bifurcations (bifbranch)
bifbranch=branch;
bifbranch.method=dde_psolfold_method(bifbranch.method);
bifbranch=replace_branch_pars(bifbranch,options.contpar,pass_on);
point=branch.point(ind);
pfoldini=dde_psolfold_init(funcs,point,options.basepar,bifbranch.method);
%% correct and add 2nd point if desired
[pfbranch,suc]=correct_ini(funcs,bifbranch,pfoldini,...
    options.dir,options.step,options.correc);
end