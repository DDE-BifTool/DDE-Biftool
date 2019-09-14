function y=rot_dirderi(xx,p,dxx,dp,A,expA,user_dirderi,user_tau,lhs_mat)
%% directional derivative or r.h.s. in rotating coordinates
[dim,d,nvec]=size(xx);
omega=p(1,end);
userpar=p(1,1:end-1);
tau_ind=user_tau();
tau=[0,p(1,tau_ind)]';   % delays including 0
dtau=[0,dp(1,tau_ind)]'; % derivatives wrt. delays including 0
dom=dp(1,end);
duserpar=dp(1,1:end-1);
rs=@(x)reshape(x,dim,nvec);
rs1=@(x)reshape(x,dim,1,nvec);
for i=d:-1:1
    dev=A*rs(xx(:,i,:))*(tau(i)*dom+dtau(i));
    xxrot(:,i,:)=expA(-omega*tau(i))*rs(  xx(:,i,:));
   dxxrot(:,i,:)=expA(-omega*tau(i))*(rs(dxx(:,i,:))-dev);
end
y0=rs1(user_dirderi(xxrot,userpar,dxxrot,duserpar));
y=y0-rs1(lhs_mat*A*(dom*rs(xxrot(:,1,:))+omega*rs(dxxrot(:,1,:))));
end