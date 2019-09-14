function data=dde_stst2psol(funcs,info,varargin)
%% constructor of periodic orbit from Hopf bifurcation
% $Id: dde_stst2psol.m 345 2019-05-12 00:11:15Z jansieber $
%%
default={'radius',0.01,'degree',3,'intervals',20};
options=dde_set_options(default,varargin,'pass_on');
hopf=info.point(1);
if ~strcmp(hopf.kind,'hopf')
    if ~isfield(hopf,'stability') || isempty(hopf.stability)
        hopf.stability=p_stabil(funcs,hopf,info.method.stability);
    end
    hopf=p_tohopf(funcs,hopf);
end
[deg_psol,step_cond]=p_topsol(funcs,hopf,...
    options.radius,options.degree,options.intervals);
info.method=df_mthod('psol');
info.method=replace_method_pars(info.method,varargin{:});
data.info=info;
data.info.point=deg_psol;
data.info.tangent=p_axpy(-1,step_cond,[]);
data.funcs=funcs;
end