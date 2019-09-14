function [J,res]=dde_fold_jac_res(funcs,pt,free_par,varargin)
%% Jacobian and residual of nonlinear system for fold
% function [J,res]=dde_fold_jac_res(funcs,pt,free_par,...)
% INPUT:
%   funcs problem functions
%	pt current fold solution guess
%	free_par free parameter numbers
%	pref: reference point
% OUTPUT: 
%	J jacobian in R^(2n+1+s x 2n+p)
%	res residual in R^(2n+1+s x 1)
 
% (c) DDE-BIFTOOL v. 2.00, 23/11/2001
%
% $Id: dde_fold_jac_res.m 348 2019-06-19 13:09:01Z jansieber $
%
%%
[ind,len]=dde_ind_from_point(pt,free_par);
dfx=dde_stst_linearize_rhs(funcs,pt);
m=size(dfx,3);
dfxsum=sum(dfx,3);
om=ones(m,1);
xx=pt.x(:,om);
res0=funcs.sys_rhs(xx,pt.parameter);...
res=[res0;...
    dfxsum*pt.v;...
    pt.v'*pt.v-1];
%% Jacobian
dfp=dde_stst_linearize_rhs(funcs,pt,'free_par',free_par);
[dfxxv,dfxpv]=sum_dfdxxv(funcs,xx,pt.parameter,pt.v,free_par);
nf=length(res0);
J=zeros(length(res),len);
J(1:nf,[ind.x(:);ind.parameter(:)])=              [ dfxsum,          dfp  ];
J(nf+(1:nf),[ind.x(:);ind.v(:);ind.parameter(:)])=[ dfxxv, dfxsum,   dfxpv];
J(end,ind.v)=                                              2*pt.v';
end
%%
function [dxx,dxp]=sum_dfdxxv(funcs,xx,p,v,free_par)
[n,m]=size(xx);
np=length(free_par);
dxp=zeros(n,np);
dxx=zeros(n,n);
for i=1:m
    for k=1:m
        dxx=dxx+funcs.sys_deri(xx,p,[i-1,k-1],[],v);
    end
    for k=1:np
        dxp(:,k)=dxp(:,k)+funcs.sys_deri(xx,p,i-1,free_par(k),[])*v;
    end
end
end
