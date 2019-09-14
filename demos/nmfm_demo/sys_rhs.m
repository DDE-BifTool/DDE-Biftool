function f = sys_rhs(xx,par)

f(1,1) = -xx(1,1)-par(1)*(tanh(par(2)*xx(1,2)-1)+tanh(1))*cosh(1)^2+par(3)*(tanh(par(4)*xx(2,3)-1)+tanh(1))*cosh(1)^2;
f(2,1) = -xx(2,1)-par(1)*(tanh(par(2)*xx(2,2)-1)+tanh(1))*cosh(1)^2+par(3)*(tanh(par(4)*xx(1,3)-1)+tanh(1))*cosh(1)^2;
end
