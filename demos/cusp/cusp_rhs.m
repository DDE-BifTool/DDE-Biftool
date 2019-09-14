function f = cusp_rhs(xx,par)
%% right-hand side of cusp demo (see cusp_demo.m for equations)
%
% $Id: cusp_rhs.m 176 2017-03-13 00:25:33Z jansieber $
%
%par = [q11,q12,q21,e1,e2]

alpha = 1./(1+exp(-4*xx(1,2,:)))-1/2;

f(1,:) = -xx(1,1,:)+par(1)*alpha-par(2)*xx(2,2,:)+par(4);
f(2,:) = -xx(2,1,:)+par(3)*alpha+par(5);
end
