%% Codimension-2 normal form computation (wrapping for finite differences)
%% Input
%
% * funcs: problem functions
% * pt: point that is (approx) codim-2 point
% * nmfm_compute: function of type newpoint=nmfm_type(funcs,point) where
% type is hoho, zeho, BT, cusp etc
%% Output
%
% * nf: point with normal form
% * nflow: if numerical finite differences are used then computation is
% done twice, once with higher order, once with lower order, this output is
% the result with lower order. Use the difference between nf and nflow to
% estimate the error
%
%%
function [nf,nflow]=nmfm_wrap_findiff(funcs,pt,nmfm_compute)
if funcs.sys_mfderi_provided
    fin_diff=false;
else
    fin_diff=true;
    funcs1=funcs;
    funcs2=set_funcs(funcs1,'sys_mfderi',{'output',2});
end
if fin_diff
    newpoint=nmfm_compute(funcs1,pt);
    nf=newpoint;
    newpoint=nmfm_compute(funcs2,pt);
    nflow=newpoint;
else
    newpoint=nmfm_compute(funcs,pt);
    nf=newpoint;
    nflow=nf;
end
end
