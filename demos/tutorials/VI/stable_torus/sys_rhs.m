function f = sys_rhs(xx,par)
%%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: sys_rhs.m 134 2016-09-12 11:10:44Z mmbosschaert $
%%
f = -par(4).*xx(1,1,:)-par(1).*xx(1,2,:)-par(2).*xx(1,3,:)-par(1).*par(3).*xx(1,1,:).*(par(1).*xx(1,4,:)+par(2).*xx(1,5,:)+par(4).*xx(1,2,:))-par(2).*par(3).*xx(1,1,:).*(par(1).*xx(1,5,:)+par(2).*xx(1,6,:)+par(4).*xx(1,3,:))-par(1).^2.*par(3).^2.*xx(1,1,:).*xx(1,2,:).*(par(1).*xx(1,7,:)+par(2).*xx(1,8,:)+par(4).*xx(1,4,:))-par(1).*par(2).*par(3).^2.*xx(1,1,:).*xx(1,2,:).*(par(1).*xx(1,8,:)+par(2).*xx(1,9,:)+par(4).*xx(1,5,:))-par(2).*par(1).*par(3).^2.*xx(1,1,:).*xx(1,3,:).*(par(1).*xx(1,8,:)+par(2).*xx(1,9,:)+par(4).*xx(1,5,:))-par(2).^2.*par(3).^2.*xx(1,1,:).*xx(1,3,:).*(par(1).*xx(1,9,:)+par(2).*xx(1,10,:)+par(4).*xx(1,6,:))-1/2.*par(3).^2.*xx(1,1,:).^2.*par(1).*(par(4).^2.*xx(1,2,:)+2.*par(4).*(par(1).*xx(1,4,:)+par(2).*xx(1,5,:))+par(1).^2.*xx(1,7,:)+2.*par(1).*par(2).*xx(1,8,:)+par(2).^2.*xx(1,9,:))-1/2.*par(3).^2.*xx(1,1,:).^2.*par(2).*(par(4).^2.*xx(1,3,:)+2.*par(4).*(par(1).*xx(1,5,:)+par(2).*xx(1,6,:))+par(1).^2.*xx(1,8,:)+2.*par(1).*par(2).*xx(1,9,:)+par(2).^2.*xx(1,10,:));

end
