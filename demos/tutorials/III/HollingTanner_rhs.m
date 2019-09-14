%% Right-hand side of Holling=Tanner model with delay
function f = HollingTanner_rhs(xx,par)
%%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: HollingTanner_rhs.m 134 2016-09-12 11:10:44Z mmbosschaert $
%%

f(1,:) = (xx(1,1,:)+par(4)).*(1-xx(1,1,:)-par(4))-...
    xx(1,1,:).*xx(2,1,:)./(par(3)*xx(2,1,:)+xx(1,1,:))-par(5);
f(2,:) = par(6)*xx(2,1,:).*(par(1)-xx(2,2,:)./xx(1,2,:));
end
