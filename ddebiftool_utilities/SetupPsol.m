%% Start branch of periodic orbits from Hopf point or bifurcation of periodic orbit branch
%%
function [per,suc]=SetupPsol(funcs,branch,ind,varargin)
%% Inputs
% 
% * |funcs|: structure with functions provided by user
% * |branch|: branch of |'stst'| steady state solutions, |'hopf'|
% solutions from which one wants to branch off, or |'psol'| solutions that are bifurcations of psols
% * |ind|: index in |'point'| field of |branch| where we want to branch off
% 
% Important optional inputs (name-value pairs)
%
% (for branching off from Hopf)
%
% * |'degree'|: degree used for collocation
% * |'intervals'|: number of collocation intervals (overall mesh size is
% degree x intervals + 1
% * |'hopfcorrection'|: whether Hopf point still needs to be corrected
%
% (general)
%
% * |'contpar'|: index of continuation parameters if different from
% |'branch'|, (this should be set manually when the original branch was a
% two-parameter branch)
% * |'corpar'|: parameters left free for initial correction (if different
% from |contpars|)
%
% All other optional inputs are passed on to fields of out branch |per|
%% Outputs
%
% * |per|: branch of periodic orbits with desired settings and two initial
% corrected points
% * |suc|: flag indicating success
% 
% $Id: SetupPsol.m 309 2018-10-28 19:02:42Z jansieber $
%
%%
default={'contpar',[],'corpar',[]};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
if isempty(options.contpar)
    options.contpar=branch.parameter.free;
end
if isempty(options.corpar)
    % assume that initial correction have to be made in all continuation
    % parameters (to accomodate step condition)
    options.corpar=options.contpar;
end
kind=branch.point(ind).kind;
[per,suc]=dde_apply({'SetupPsolFrom_',kind,''},funcs,branch,ind,...
    'contpar',options.contpar,'corpar',options.corpar,pass_on{:});
end
