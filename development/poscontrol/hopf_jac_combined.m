function [J,res]=hopf_jac_combined(funcs,point,free_par_ind,varargin)
default={'reference',point,'extra_condition',repmat(struct(),1,0)};
options=dde_set_options(default,varargin,'pass_on');
c=options.reference.v;
c=c'/norm(c);
[J,res]=hopf_jac(funcs,point.x,point.omega,point.v,...
    point.parameter,free_par_ind,c);
end

