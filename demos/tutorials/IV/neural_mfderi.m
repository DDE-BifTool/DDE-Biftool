function y = neural_mfderi(xx,par,varargin)
%%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: neural_mfderi.m 134 2016-09-12 11:10:44Z mmbosschaert $
%%

if nargin == 2
	error('SYS_MFDERI: no arguments.');
elseif nargin > 7
	error('SYS_MFDERI: too many arguments.');
end

y = 0;

numarg = nargin - 2;

switch numarg
	case 1
		u1 = varargin{1};
		y = [-u1(1,1)-par(1)*par(2)*(1-tanh(par(2)*xx(1,2)-1)^2)*cosh(1)^2*u1(1,2)+par(3)*par(4)*(1-tanh(par(4)*xx(2,3)-1)^2)*cosh(1)^2*u1(2,3);
			-u1(2,1)-par(1)*par(2)*(1-tanh(par(2)*xx(2,2)-1)^2)*cosh(1)^2*u1(2,2)+par(3)*par(4)*(1-tanh(par(4)*xx(1,3)-1)^2)*cosh(1)^2*u1(1,3)];
	case 2
		u1 = varargin{1}; u2 = varargin{2};
		y = [2*par(1)*par(2)^2*tanh(par(2)*xx(1,2)-1)*(1-tanh(par(2)*xx(1,2)-1)^2)*cosh(1)^2*u1(1,2)*u2(1,2)-2*par(3)*par(4)^2*tanh(par(4)*xx(2,3)-1)*(1-tanh(par(4)*xx(2,3)-1)^2)*cosh(1)^2*u1(2,3)*u2(2,3);
			2*par(1)*par(2)^2*tanh(par(2)*xx(2,2)-1)*(1-tanh(par(2)*xx(2,2)-1)^2)*cosh(1)^2*u1(2,2)*u2(2,2)-2*par(3)*par(4)^2*tanh(par(4)*xx(1,3)-1)*(1-tanh(par(4)*xx(1,3)-1)^2)*cosh(1)^2*u1(1,3)*u2(1,3)];
	case 3
		u1 = varargin{1}; u2 = varargin{2}; u3 = varargin{3};
		y = [(2*par(1)*par(2)^3*(1-tanh(par(2)*xx(1,2)-1)^2)^2*cosh(1)^2-4*par(1)*par(2)^3*tanh(par(2)*xx(1,2)-1)^2*(1-tanh(par(2)*xx(1,2)-1)^2)*cosh(1)^2)*u1(1,2)*u2(1,2)*u3(1,2)+(-2*par(3)*par(4)^3*(1-tanh(par(4)*xx(2,3)-1)^2)^2*cosh(1)^2+4*par(3)*par(4)^3*tanh(par(4)*xx(2,3)-1)^2*(1-tanh(par(4)*xx(2,3)-1)^2)*cosh(1)^2)*u1(2,3)*u2(2,3)*u3(2,3);
			(2*par(1)*par(2)^3*(1-tanh(par(2)*xx(2,2)-1)^2)^2*cosh(1)^2-4*par(1)*par(2)^3*tanh(par(2)*xx(2,2)-1)^2*(1-tanh(par(2)*xx(2,2)-1)^2)*cosh(1)^2)*u1(2,2)*u2(2,2)*u3(2,2)+(-2*par(3)*par(4)^3*(1-tanh(par(4)*xx(1,3)-1)^2)^2*cosh(1)^2+4*par(3)*par(4)^3*tanh(par(4)*xx(1,3)-1)^2*(1-tanh(par(4)*xx(1,3)-1)^2)*cosh(1)^2)*u1(1,3)*u2(1,3)*u3(1,3)];
	case 4
		u1 = varargin{1}; u2 = varargin{2}; u3 = varargin{3}; u4 = varargin{4};
		y = [(-16*par(1)*par(2)^4*(1-tanh(par(2)*xx(1,2)-1)^2)^2*cosh(1)^2*tanh(par(2)*xx(1,2)-1)+8*par(1)*par(2)^4*tanh(par(2)*xx(1,2)-1)^3*(1-tanh(par(2)*xx(1,2)-1)^2)*cosh(1)^2)*u1(1,2)*u2(1,2)*u3(1,2)*u4(1,2)+(16*par(3)*par(4)^4*(1-tanh(par(4)*xx(2,3)-1)^2)^2*cosh(1)^2*tanh(par(4)*xx(2,3)-1)-8*par(3)*par(4)^4*tanh(par(4)*xx(2,3)-1)^3*(1-tanh(par(4)*xx(2,3)-1)^2)*cosh(1)^2)*u1(2,3)*u2(2,3)*u3(2,3)*u4(2,3);
			(-16*par(1)*par(2)^4*(1-tanh(par(2)*xx(2,2)-1)^2)^2*cosh(1)^2*tanh(par(2)*xx(2,2)-1)+8*par(1)*par(2)^4*tanh(par(2)*xx(2,2)-1)^3*(1-tanh(par(2)*xx(2,2)-1)^2)*cosh(1)^2)*u1(2,2)*u2(2,2)*u3(2,2)*u4(2,2)+(16*par(3)*par(4)^4*(1-tanh(par(4)*xx(1,3)-1)^2)^2*cosh(1)^2*tanh(par(4)*xx(1,3)-1)-8*par(3)*par(4)^4*tanh(par(4)*xx(1,3)-1)^3*(1-tanh(par(4)*xx(1,3)-1)^2)*cosh(1)^2)*u1(1,3)*u2(1,3)*u3(1,3)*u4(1,3)];
	case 5
		u1 = varargin{1}; u2 = varargin{2}; u3 = varargin{3}; u4 = varargin{4}; u5 = varargin{5};
		y = [(88*par(1)*par(2)^5*(1-tanh(par(2)*xx(1,2)-1)^2)^2*cosh(1)^2*tanh(par(2)*xx(1,2)-1)^2-16*par(1)*par(2)^5*(1-tanh(par(2)*xx(1,2)-1)^2)^3*cosh(1)^2-16*par(1)*par(2)^5*tanh(par(2)*xx(1,2)-1)^4*(1-tanh(par(2)*xx(1,2)-1)^2)*cosh(1)^2)*u1(1,2)*u2(1,2)*u3(1,2)*u4(1,2)*u5(1,2)+(-88*par(3)*par(4)^5*(1-tanh(par(4)*xx(2,3)-1)^2)^2*cosh(1)^2*tanh(par(4)*xx(2,3)-1)^2+16*par(3)*par(4)^5*(1-tanh(par(4)*xx(2,3)-1)^2)^3*cosh(1)^2+16*par(3)*par(4)^5*tanh(par(4)*xx(2,3)-1)^4*(1-tanh(par(4)*xx(2,3)-1)^2)*cosh(1)^2)*u1(2,3)*u2(2,3)*u3(2,3)*u4(2,3)*u5(2,3);
			(88*par(1)*par(2)^5*(1-tanh(par(2)*xx(2,2)-1)^2)^2*cosh(1)^2*tanh(par(2)*xx(2,2)-1)^2-16*par(1)*par(2)^5*(1-tanh(par(2)*xx(2,2)-1)^2)^3*cosh(1)^2-16*par(1)*par(2)^5*tanh(par(2)*xx(2,2)-1)^4*(1-tanh(par(2)*xx(2,2)-1)^2)*cosh(1)^2)*u1(2,2)*u2(2,2)*u3(2,2)*u4(2,2)*u5(2,2)+(-88*par(3)*par(4)^5*(1-tanh(par(4)*xx(1,3)-1)^2)^2*cosh(1)^2*tanh(par(4)*xx(1,3)-1)^2+16*par(3)*par(4)^5*(1-tanh(par(4)*xx(1,3)-1)^2)^3*cosh(1)^2+16*par(3)*par(4)^5*tanh(par(4)*xx(1,3)-1)^4*(1-tanh(par(4)*xx(1,3)-1)^2)*cosh(1)^2)*u1(1,3)*u2(1,3)*u3(1,3)*u4(1,3)*u5(1,3)];
end
