function f = sys_rhs(xx,par)

f(1) =par(1) + xx(2,1) - par(5)*xx(3,1) - par(3)*xx(1,1)^3 + par(4)*xx(1,2)^2;
f(2) =par(5) - xx(2,1) - par(6)*xx(1,1)^2;
f(3) =-par(8)*(xx(3,1) + par(2)*(par(7) - xx(1,1)));

end
