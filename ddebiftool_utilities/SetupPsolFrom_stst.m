%% Branch off at Hopf point (point is either of kind hopf or stst)
%%
function [per,suc]=SetupPsolFrom_stst(funcs,ststbranch,ind,varargin)
%% Inputs
% 
% * |funcs|: structure with functions provided by user
% * |ststbranch|: branch of |'stst'| steady state solutions or |'hopf'|
% solutions from which one wants to branch off
% * |ind|: index in |'point'| field of |ststbranch| which is close to Hopf
% point where we want to branch off
% 
% Important optional inputs (name-value pairs)
%
% * |'degree'|: degree used for collocation
% * |'intervals'|: number of collocation intervals (overall mesh size is
% degree x intervals + 1
% * |'hopfcorrection'|: whether Hopf point still needs to be corrected
% * |'contpar'|: index of continuation parameters if different from
% |'ststbranch'|
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
% $Id: SetupPsolFrom_stst.m 369 2019-08-27 00:07:02Z jansieber $
%
%%
default={'radius',0.01,'contpar',[],'corpar',[],...
    'degree',3,'intervals',20,'hopfcorrection',true};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
% create branch per of periodic solutions starting from an
% approximate Hopf point num on a branch br of steady states
if isempty(options.contpar)
    options.contpar=ststbranch.parameter.free;
end
if isempty(options.corpar)
    % assume that initial correction have to be made in all continuation
    % parameters (to accomodate step condition)
    options.corpar=options.contpar;
end
hopf=ststbranch.point(ind);
if ~strcmp(hopf.kind,'hopf')
    if ~isfield(hopf,'stability') || isempty(hopf.stability)
        hopf.stability=p_stabil(funcs,hopf,ststbranch.method.stability);
    end
    hopf=p_tohopf(funcs,hopf);
end
hm=df_mthod(funcs,'hopf');
hm.point=dde_set_options(hm.point,pass_on,'pass_on');
if options.hopfcorrection
    [hopf,suc]=p_correc(funcs,hopf,options.corpar,[],hm.point);
    if suc==0
        per=[];
        return;
    end
end
[deg_psol,step_cond]=p_topsol(funcs,hopf,...
    options.radius,options.degree,options.intervals,pass_on{:});
per=df_brnch(funcs,options.contpar,'psol');
per.parameter=ststbranch.parameter;
per.parameter.free=options.contpar;

per=replace_branch_pars(per,options.contpar,pass_on);
[psol,suc]=p_correc(funcs,deg_psol,options.contpar,step_cond,per.method.point,1, ...
    deg_psol);
if suc==0
    per=[];
    return;
end
per.point=deg_psol;
per.point.profile=repmat(hopf.x,[1,size(per.point.profile,2)]);
per.point(2)=psol;
end
