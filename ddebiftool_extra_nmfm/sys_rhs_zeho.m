function res=sys_rhs_zeho(xx,par,funcs)
%% r.h.s. for Zero-Hopf bifurcation
%
% $Id: sys_rhs_zeho.m 309 2018-10-28 19:02:42Z jansieber $
%%
n=size(xx,1)/4;
user_rhs=funcs.sys_rhs;
q0=xx(n+(1:n),1);
q1r=xx(2*n+(1:n),1);
q1i=xx(3*n+(1:n),1);
q1=q1r+1i*q1i;
omega=par(end);
D0=ch_matrix(funcs,xx(1:n,1),par(1:end-1),0);
D_om=ch_matrix(funcs,xx(1:n,1),par(1:end-1),1i*omega);
res(1:n,1)=user_rhs(xx(1:n,:),par(1:end-1));
res(n+1:2*n,1)=D0*q0;
rz=D_om*q1;
res(2*n+1:4*n,1)=[real(rz);imag(rz)];
end
