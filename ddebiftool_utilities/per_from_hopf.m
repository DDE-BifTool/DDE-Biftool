function [per,suc]=per_from_hopf(funcs,point,varargin)
%% branch off at Hopf point (point is either of kind hopf or stst)
%
% $Id: per_from_hopf.m 309 2018-10-28 19:02:42Z jansieber $
%
%%
default={'radius',0.01,'contpar',[],...
    'degree',3,'intervals',20,'hopfcorrection',true};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
% create branch per of periodic solutions starting from an
% approximate Hopf point num on a branch br of steady states
hopf=point;
if ~strcmp(hopf.kind,'hopf')
    hopf=p_tohopf(funcs,hopf);
end;
hm=df_mthod(funcs,'hopf');
hm.point=dde_set_options(hm.point,pass_on,'pass_on');
if options.hopfcorrection
    [hopf,suc]=p_correc(funcs,hopf,options.contpar,[],hm.point);
    if suc==0
        per=[];
        return;
    end;
end
[deg_psol,step_cond]=p_topsol(funcs,hopf,...
    options.radius,options.degree,options.intervals);
pm=df_mthod(funcs,'psol');
[pm.point,pass_on]=dde_set_options(pm.point,pass_on,'pass_on');
[pm.continuation,pass_on]=dde_set_options(pm.continuation,pass_on,'pass_on');
[psol,suc]=p_correc(funcs,deg_psol,options.contpar,step_cond,pm.point,1, ...
    deg_psol);
if suc==0
    per=[];
    return;
end;
per=df_brnch(funcs,options.contpar,'psol');
per.method=pm;
per.parameter=dde_set_options(per.parameter,pass_on,'pass_on');
per.point=deg_psol;
per.point.profile=repmat(mean(per.point.profile,2),[1,size(per.point.profile,2)]);
per.point(2)=psol;
end
