function res=sys_rhs_BT(xx,par,funcs)
%% r.h.s. for Takens-Bogdanov bifurcation
% sys_cond_BT needs to be appended
%
% $Id: sys_rhs_BT.m 309 2018-10-28 19:02:42Z jansieber $
%%
n=size(xx,1)/3;
user_rhs=funcs.sys_rhs;
user_deri=funcs.sys_deri;
if ~funcs.tp_del
    tau=[0,par(funcs.sys_tau())];
else
    tau=zeros(1,funcs.sys_ntau());
    for j=2:length(tau)
        tau(j)=funcs.sys_tau(j-1,xx(1:n,ones(1,j-1)),par);
    end
end
m=length(tau)-1;
q0=xx(n+(1:n),1);
q1=xx(2*n+(1:n),1);
D=zeros(n);
dD=eye(n);
for k=0:m
    B=user_deri(xx(1:n,:),par,k,[],[]);
    D=D-B;
    dD=dD+tau(k+1)*B;
end
res(1:n,1)=user_rhs(xx(1:n,:),par);
res(n+1:2*n,1)=D*q0;
res(2*n+1:3*n,1)=dD*q0+D*q1;
end
