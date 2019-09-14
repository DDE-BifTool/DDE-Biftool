%% Initialize continuation of Fold bifurcations
%%
function [foldbranch,suc]=SetupFold(funcs,branch,ind,varargin)
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
% * |'contpar'| (integers default |[]|): index of additional continuation parameters
%  (in addition to free pars of branch
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
% $Id: SetupFold.m 309 2018-10-28 19:02:42Z jansieber $
%
%% process options
%% wrapper around SetupStstBifurcation to ensure backward compatibility
[foldbranch,suc]=SetupStstBifurcation(funcs,branch,ind,'fold',varargin{:});
end
