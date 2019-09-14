function y = sys_mfderi(xx,par,varargin)

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
		y = [-par(1)*u1(1,1)+par(3)*(1-tanh(xx(2,3))^2)*u1(2,3)+par(2)*(1-tanh(xx(1,4))^2)*u1(1,4);
			-par(1)*u1(2,1)+par(4)*(1-tanh(xx(1,2))^2)*u1(1,2)+par(2)*(1-tanh(xx(2,4))^2)*u1(2,4)];
	case 2
		u1 = varargin{1}; u2 = varargin{2};
		y = [-2*par(3)*tanh(xx(2,3))*(1-tanh(xx(2,3))^2)*u1(2,3)*u2(2,3)-2*par(2)*tanh(xx(1,4))*(1-tanh(xx(1,4))^2)*u1(1,4)*u2(1,4);
			-2*par(4)*tanh(xx(1,2))*(1-tanh(xx(1,2))^2)*u1(1,2)*u2(1,2)-2*par(2)*tanh(xx(2,4))*(1-tanh(xx(2,4))^2)*u1(2,4)*u2(2,4)];
	case 3
		u1 = varargin{1}; u2 = varargin{2}; u3 = varargin{3};
		y = [(-2*par(3)*(1-tanh(xx(2,3))^2)^2+4*par(3)*tanh(xx(2,3))^2*(1-tanh(xx(2,3))^2))*u1(2,3)*u2(2,3)*u3(2,3)+(-2*par(2)*(1-tanh(xx(1,4))^2)^2+4*par(2)*tanh(xx(1,4))^2*(1-tanh(xx(1,4))^2))*u1(1,4)*u2(1,4)*u3(1,4);
			(-2*par(4)*(1-tanh(xx(1,2))^2)^2+4*par(4)*tanh(xx(1,2))^2*(1-tanh(xx(1,2))^2))*u1(1,2)*u2(1,2)*u3(1,2)+(-2*par(2)*(1-tanh(xx(2,4))^2)^2+4*par(2)*tanh(xx(2,4))^2*(1-tanh(xx(2,4))^2))*u1(2,4)*u2(2,4)*u3(2,4)];
	case 4
		u1 = varargin{1}; u2 = varargin{2}; u3 = varargin{3}; u4 = varargin{4};
		y = [(16*par(3)*(1-tanh(xx(2,3))^2)^2*tanh(xx(2,3))-8*par(3)*tanh(xx(2,3))^3*(1-tanh(xx(2,3))^2))*u1(2,3)*u2(2,3)*u3(2,3)*u4(2,3)+(16*par(2)*(1-tanh(xx(1,4))^2)^2*tanh(xx(1,4))-8*par(2)*tanh(xx(1,4))^3*(1-tanh(xx(1,4))^2))*u1(1,4)*u2(1,4)*u3(1,4)*u4(1,4);
			(16*par(4)*(1-tanh(xx(1,2))^2)^2*tanh(xx(1,2))-8*par(4)*tanh(xx(1,2))^3*(1-tanh(xx(1,2))^2))*u1(1,2)*u2(1,2)*u3(1,2)*u4(1,2)+(16*par(2)*(1-tanh(xx(2,4))^2)^2*tanh(xx(2,4))-8*par(2)*tanh(xx(2,4))^3*(1-tanh(xx(2,4))^2))*u1(2,4)*u2(2,4)*u3(2,4)*u4(2,4)];
	case 5
		u1 = varargin{1}; u2 = varargin{2}; u3 = varargin{3}; u4 = varargin{4}; u5 = varargin{5};
		y = [(-88*par(3)*(1-tanh(xx(2,3))^2)^2*tanh(xx(2,3))^2+16*par(3)*(1-tanh(xx(2,3))^2)^3+16*par(3)*tanh(xx(2,3))^4*(1-tanh(xx(2,3))^2))*u1(2,3)*u2(2,3)*u3(2,3)*u4(2,3)*u5(2,3)+(-88*par(2)*(1-tanh(xx(1,4))^2)^2*tanh(xx(1,4))^2+16*par(2)*(1-tanh(xx(1,4))^2)^3+16*par(2)*tanh(xx(1,4))^4*(1-tanh(xx(1,4))^2))*u1(1,4)*u2(1,4)*u3(1,4)*u4(1,4)*u5(1,4);
			(-88*par(4)*(1-tanh(xx(1,2))^2)^2*tanh(xx(1,2))^2+16*par(4)*(1-tanh(xx(1,2))^2)^3+16*par(4)*tanh(xx(1,2))^4*(1-tanh(xx(1,2))^2))*u1(1,2)*u2(1,2)*u3(1,2)*u4(1,2)*u5(1,2)+(-88*par(2)*(1-tanh(xx(2,4))^2)^2*tanh(xx(2,4))^2+16*par(2)*(1-tanh(xx(2,4))^2)^3+16*par(2)*tanh(xx(2,4))^4*(1-tanh(xx(2,4))^2))*u1(2,4)*u2(2,4)*u3(2,4)*u4(2,4)*u5(2,4)];
end