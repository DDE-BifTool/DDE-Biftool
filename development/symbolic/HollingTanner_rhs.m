%% Right-hand side of Holling=Tanner model with delay
function f = HollingTanner_rhs(xx,par)
%
% $Id: HollingTanner_rhs.m 168 2017-03-03 22:04:32Z jansieber $
%
f(1,:) = (xx(1,1,:)+par(4)).*(1-xx(1,1,:)-par(4))-...
    xx(1,1,:).*xx(2,1,:)./(par(3)*xx(2,1,:)+xx(1,1,:))-par(5);
f(2,:) = par(6)*xx(2,1,:).*(par(1)-xx(2,2,:)./xx(1,2,:));
end
