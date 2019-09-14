%% create branch of psols from dde23 sol structure
%
% $Id: branch_from_sol.m 126 2016-09-05 22:37:08Z jansieber $
%
% inputs:
%
% * |funcs|: functions defining right-hand side (created with set_funcs)
% * |sol|: solutions structure output by dde23 (settled to rotation)
% * |free_par_ind|: index of parameter to be varied along branch initially
% * |parameters|: initial values for parameters
%
% Important optional inputs (name-value pairs, all passed on to )
%
% * |'degree'|: (default 3) degree used for collocation
% * |'intervals'|: (default 20) number of collocation intervals (overall mesh size is
% degree x intervals + 1
% * |'stepcond'|: (default []) passed on to p_correc
% * |'step'|: (default 1e-2) initial step between first two points
% * |'corpar'|: parameters left free for initial correction (if different
% from |free_par_ind|)
%
%%
function branch=branch_from_sol(funcs,sol,free_par_ind,parameters,varargin)
default={'degree',3,'intervals',20,'corpar',[],'step',1e-2,'indperiod',[],'stepcond',[]};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
if ~isfield(sol,'xe') || length(sol.xe)<2
    error('branch_from_sol:cross_section',...
        'solution structure does not contain two Poincare cross section events');
end
period=sol.xe(end)-sol.xe(end-1);
per0=struct('kind','psol','parameter',parameters,'mesh',[],'degree',options.degree,...
    'profile',[],'period',period);
if ~isempty(options.indperiod)
    per0.parameter(options.indperiod)=period;
end
per0.mesh=linspace(0,1,options.degree*options.intervals+1);
solmesh=sol.xe(end-1)+per0.mesh*per0.period;
per0.profile=deval(sol,solmesh);
per0.profile=per0.profile-floor(per0.profile(1)/(2*pi))*2*pi;
per1=p_remesh(per0,options.degree,options.intervals);
branch=df_brnch(funcs,free_par_ind,'psol');
branch=replace_branch_pars(branch,free_par_ind,pass_on);
corpar=union(options.corpar,options.indperiod);
[per2,suc]=p_correc(funcs,per1,corpar,options.stepcond,branch.method.point);
if ~suc
    error('branch_from_sol:correction',...
        '1st correction failed: suc=%d',suc);
end
per3=per2;
per3.parameter(free_par_ind(1))=per2.parameter(free_par_ind(1))+options.step;
[per3,suc]=p_correc(funcs,per3,corpar,options.stepcond,branch.method.point);
if ~suc
    error('branch_from_sol:correction',...
        '2nd correction failed: suc=%d',suc);
end
branch.point=[per2,per3];
end

