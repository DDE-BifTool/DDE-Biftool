function f=sd_rhs(xx,par)
%% Right-hand side for sd_demo (vectorized in x)
% p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11

f=NaN(size(xx,1),size(xx,3));
f(1,:)=(1./(par(1)+xx(2,1,:))).*(1-par(2)*xx(1,1,:).*xx(1,4,:).*xx(3,4,:)+par(3)*xx(1,2,:).*xx(2,3,:));
f(2,:)=par(4)*xx(1,1,:)./(par(1)+xx(2,1,:))+par(5)*tanh(xx(2,6,:))-1;
f(3,:)=par(6)*(xx(2,1,:)-xx(3,1,:))-par(7)*(xx(1,7,:)-xx(2,5,:)).*exp(-par(8)*xx(4,1,:));
f(4,:)=xx(1,5,:).*exp(-par(1)*xx(4,1,:))-0.1;
f(5,:)=3*(xx(1,3,:)-xx(5,1,:))-par(9); 
end
