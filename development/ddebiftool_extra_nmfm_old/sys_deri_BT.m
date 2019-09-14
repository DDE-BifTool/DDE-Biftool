function J=sys_deri_BT(xx,par,nx,np,v,funcs) %#ok<INUSL>
%% Jacobian for Takens-Bogdanov bifurcation (fitting to sys_rhs_BT)
% (see sys_rhs_BT for formula)
%
% $Id: sys_deri_BT.m 109 2015-08-31 23:45:11Z jansieber $
%%
n=size(xx,1)/3;
user_deri=funcs.sys_deri;
ind_tau=funcs.sys_tau();
tau=[0,par(ind_tau)];
m=length(tau)-1;
x=xx(1:n,:);
q0=xx(n+1:2*n,1);
q1=xx(2*n+1:3*n,1);
D=zeros(n);
dD=eye(n);
df_dxk=zeros(n,n,m+1);
for k=0:m
    df_dxk(:,:,k+1)=user_deri(x,par,k,[],[]);
    D=D-df_dxk(:,:,k+1);
    dD=dD+tau(k+1)*df_dxk(:,:,k+1);
end
if length(nx)==1 && isempty(np) % deriv wrt x(:,nx+1)
    J=zeros(3*n,3*n);
    J(1:n,1:n)=df_dxk(:,:,nx+1);
    if nx==0
        J(n+(1:n),n+(1:n))=D;
        J(2*n+(1:n),n+(1:n))=dD;
        J(2*n+(1:n),2*n+(1:n))=D;
    end
    for k=0:m
        dfdx_nxk_q0=user_deri(x,par,[k,nx],[],q0);
        dfdx_nxk_q1=user_deri(x,par,[k,nx],[],q1);
        J(n+(1:n),1:n)=J(n+(1:n),1:n)-dfdx_nxk_q0;
        J(2*n+(1:n),1:n)=J(2*n+(1:n),1:n)-dfdx_nxk_q1;
        J(2*n+(1:n),1:n)=J(2*n+(1:n),1:n)+tau(k+1)*dfdx_nxk_q0;
    end
elseif isempty(nx) && ~isempty(np) % deriv wrt par
    J=zeros(3*n,1);
    J(1:n)=user_deri(x,par,[],np,[]);
    for k=0:m
        dfdxp=user_deri(x,par,k,np,[]);
        J(n+(1:n))=J(n+(1:n))-dfdxp*q0;
        J(2*n+(1:n))=J(2*n+(1:n))-dfdxp*q1+tau(k+1)*dfdxp*q0;
    end
    [istau,numtau]=ismember(np,ind_tau);
    if istau
        J(2*n+(1:n))=J(2*n+(1:n))+df_dxk(:,:,numtau+1)*q0;
    end
else
    error('sys_deri_BT:call','sys_deri_BT called for 2nd derivative');
end
end
