function F=nmfm_deriv_define(funcs,point,varargin)
default={'free_pars',[]};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
F.B=@(v1,v2)nmfm_mfderiv(funcs,point,[v1,v2],pass_on{:});
F.C=@(v1,v2,v3)nmfm_mfderiv(funcs,point,[v1,v2,v3],pass_on{:});
F.D=@(v1,v2,v3,v4)nmfm_mfderiv(funcs,point,[v1,v2,v3,v4],pass_on{:});
F.E=@(v1,v2,v3,v4,v5)nmfm_mfderiv(funcs,point,[v1,v2,v3,v4,v5],pass_on{:});
if ~isempty(options.free_pars)
    par0=point.parameter'*0;
    par=@(vals)dde_array_insert(par0,options.free_pars,vals);
    dp=@(alpha)nmfm_dev_fun([0*point.x;par(alpha)]);
    F.J1=[nmfm_mfderiv(funcs,point,dp([1;0]),pass_on{:}),...
        nmfm_mfderiv(funcs,point,dp([0;1]),pass_on{:})];
    F.A1=@(v1,v2) nmfm_mfderiv(funcs,point,[v1,dp(v2)],pass_on{:});
    F.B1=@(v1,v2,v3)  nmfm_mfderiv(funcs,point,[v1,v2,dp(v3)],pass_on{:});
    F.C1=@(v1,v2,v3,v4) nmfm_mfderiv(funcs,point,[v1,v2,v3,dp(v4)],pass_on{:});
end
end