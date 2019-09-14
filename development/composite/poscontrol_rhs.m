function f=poscontrol_rhs(xx,par,ind)

% p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11
% tau s0 c K
s0=reshape(par(1,ind.s0,:),1,[]);
K=reshape(par(1,ind.K,:),1,[]);
c=reshape(par(1,ind.c,:),1,[]);
ntaup1=size(xx,2);
x=reshape(xx(1,:,:),ntaup1,[]);
s=reshape(xx(2,:,:),ntaup1,[]);
kc2=K.*c/2;
f(1,:)=-kc2.*(s(2,:)-s0);
num=-kc2.*(s(4,:)+s(2,:)-2*s0)-c.*s(1,:)+x(3,:)+x(1,:);
den=c+kc2.*(s(4,:)-s0);
f(2,:)=num./den;
end



