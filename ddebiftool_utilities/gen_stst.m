function [br,suc]=gen_stst(funcs,varargin)
%% generate initial points of stst branch
%
% $Id: gen_stst.m 309 2018-10-28 19:02:42Z jansieber $
%
default={'step',0.01,'contpar',[],'corpar',[]};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
if isempty(options.corpar)
    % assume that initial correction have to be made in all parameters but
    % first
    options.corpar=options.contpar(2:end);
end
% create branch br of steady state solutions traversing through
% point changing par in direction dir
br=[];
point=struct('kind','stst','parameter',[],'x',[],'stability',[]);
[point,pass_on]=dde_set_options(point,pass_on,'pass_on');
%first step: correct point
if isfield(pass_on,'newheuristics_tests') && pass_on.newheuristics_tests==0
    mth=df_mthod(funcs,'stst',0);
else
    mth=df_mthod(funcs,'stst',0);
end
[mth.point,pass_on]=dde_set_options(mth.point,pass_on,'pass_on');
[mth.continuation,pass_on]=dde_set_options(mth.continuation,pass_on,'pass_on');
[mth.stability,pass_on]=dde_set_options(mth.stability,pass_on,'pass_on');

[point,suc]=p_correc(funcs,point,options.corpar,[],mth.point);
if ~suc
    return
end
if isempty(options.contpar)
    br=point;
    return
end
%define and initialize branch
br=df_brnch(funcs,options.contpar,'stst');
br.parameter=dde_set_options(br.parameter,pass_on,'pass_on');
br.method=mth;

br.point=point;

%find second point on the branch
p2=point;
p2.parameter(options.contpar(1))=p2.parameter(options.contpar(1))+options.step;
[p2,suc]=p_correc(funcs,p2,options.corpar,[],mth.point);
if ~suc
    return
end
br.point(2)=p2;
end
