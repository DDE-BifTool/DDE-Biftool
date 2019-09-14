function y=sys_rhs_TorusBif(x,p,ip,funcs,sys_deri,lhs_matrix)
%% r.h.s. of extended DDE for torus and period doubling bifurcation
% argument sys_deri can be either numeric (say, 1e-3, then a
% finite-difference approximation is used for the linearized system), or a
% function providing the partial derivatives of the DDE's sys_rhs
%
% $Id: sys_rhs_TorusBif.m 369 2019-08-27 00:07:02Z jansieber $
%
dim=ip.dim;
omega=p(ip.omega);
period=p(ip.period);
tau=[0,p(ip.orig_tau)];
x0=x(1:dim,:,:);
u=x(dim+1:2*dim,:,:);
v=x(2*dim+1:end,:,:);
user_p=p(1:ip.nuserpar);
y0=funcs.sys_rhs(x0,user_p);
v1=reshape(v(:,1,:),dim,[]);
u1=reshape(u(:,1,:),dim,[]);
yu=pi*lhs_matrix*omega/period*v1;
yv=-pi*lhs_matrix*omega/period*u1;
if isnumeric(sys_deri) 
    %% no user-provided derivative (sys_deri is size of deviation)
    df=@(x0,dev,ind)app_dir_deriv(@(x)funcs.sys_rhs(x,user_p),x0,dev,ind,sys_deri);
else
    %% sys_deri is user-provided function
    df=@(x0,dev,ind)VAopX(funcs.sys_deri(x0,user_p,ind-1,[],[]),dev,'*');
end
for i=1:size(x0,2)
    %% add partial derivatives of all delayed terms (incl delay zero)
    c=cos(pi*omega/period*tau(i));
    s=sin(pi*omega./period*tau(i));
    yu=yu+reshape(df(x0, c*u(:,i,:)+s*v(:,i,:),i),dim,[]);
    yv=yv+reshape(df(x0,-s*u(:,i,:)+c*v(:,i,:),i),dim,[]);
end
y=cat(1,y0,yu,yv);
end
