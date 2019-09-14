function y=rot_rhs(xx,p,A,expA,user_rhs,user_tau,lhs_matrix)
%% right-hand side in rotating coordinates
%
% $Id: rot_rhs.m 369 2019-08-27 00:07:02Z jansieber $
%
dim=size(xx,1);
nvec=size(xx,3);
omega=p(1,end);
userpar=p(1,1:end-1);
tau_ind=user_tau();
tau=[0,p(1,tau_ind)]'; % delays including 0
rs=@(x)reshape(x,dim,nvec);
for i=size(xx,2):-1:1
    xxrot(:,i,:)=expA(-omega.*tau(i))*rs(xx(:,i,:));
end
y0=reshape(user_rhs(xxrot,userpar),[dim,1,nvec]);
y=y0-reshape(lhs_matrix*A*omega*rs(xxrot(:,1,:)),[dim,1,nvec]);
end
