%% Generate steady-state branch with initial point(s)
%%
function [br,suc]=SetupStst(funcs,varargin)
%% Inputs
%
% |funcs|: structure with functions provided by user
%
% All other inputs are name-value pairs. Important inputs:
%
% * |'parameter'|, |'x'|: passed on to initial point structure
% * |'contpar'|: index of continuation parameter(s)
% * |'corpar'|: index of parameter(s) used for initial correction (if
% different from |'contpar'|)
% * |'dir'|: if more than one index is given in |'contpar'| then the
% parameter with index |'dir'| is varied in the first step
% * |'step'|: length of initial step (can be negative, default |0.01|)
%
% All other optional inputs are passed on to fields of out branch |per|
%% Outputs
%
% * |br|: branch of steady states with desired settings and two initial corrected points
% * |suc|: flag indicating success
%
% $Id: SetupStst.m 357 2019-06-30 23:59:31Z jansieber $
%
%%
default={'step',0.01,'contpar',[],'corpar',[],'dir',[],'correc',true};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
if isempty(options.corpar)
    % assume that initial correction have to be made in all continuation parameters but
    % first
    options.corpar=options.contpar(2:end);
end
% create branch br of steady state solutions traversing through
% point changing par in direction dir
point=dde_stst_create(pass_on{:});
[point,pass_on]=dde_set_options(point,pass_on,'pass_on');
%% define and initialize branch
br=df_brnch(funcs,options.contpar,'stst');
br=replace_branch_pars(br,options.contpar,pass_on);
if isempty(options.dir)
    dir=br.parameter.free;
else
    dir=options.dir;
end
[br,suc]=correct_ini(funcs,br,point,dir,options.step,options.correc);
end
