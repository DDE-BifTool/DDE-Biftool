function y = cusp_mfderi(xx,par,varargin)
%% higher-order derivatives in xx for cusp demo (see cusp_demo.m for equations)
%
% $Id: cusp_mfderi.m 121 2015-09-02 22:02:24Z jansieber $
%
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
		y = [-u1(1,1)+4*par(1)/(1+exp(-4*xx(1,2)))^2*exp(-4*xx(1,2))*u1(1,2)-par(2)*u1(2,2);
			-u1(2,1)+4*par(3)/(1+exp(-4*xx(1,2)))^2*exp(-4*xx(1,2))*u1(1,2)];
	case 2
		u1 = varargin{1}; u2 = varargin{2};
		y = [(32*par(1)/(1+exp(-4*xx(1,2)))^3*exp(-4*xx(1,2))^2-16*par(1)/(1+exp(-4*xx(1,2)))^2*exp(-4*xx(1,2)))*u1(1,2)*u2(1,2);
			(32*par(3)/(1+exp(-4*xx(1,2)))^3*exp(-4*xx(1,2))^2-16*par(3)/(1+exp(-4*xx(1,2)))^2*exp(-4*xx(1,2)))*u1(1,2)*u2(1,2)];
	case 3
		u1 = varargin{1}; u2 = varargin{2}; u3 = varargin{3};
		y = [(384*par(1)/(1+exp(-4*xx(1,2)))^4*exp(-4*xx(1,2))^3-384*par(1)/(1+exp(-4*xx(1,2)))^3*exp(-4*xx(1,2))^2+64*par(1)/(1+exp(-4*xx(1,2)))^2*exp(-4*xx(1,2)))*u1(1,2)*u2(1,2)*u3(1,2);
			(384*par(3)/(1+exp(-4*xx(1,2)))^4*exp(-4*xx(1,2))^3-384*par(3)/(1+exp(-4*xx(1,2)))^3*exp(-4*xx(1,2))^2+64*par(3)/(1+exp(-4*xx(1,2)))^2*exp(-4*xx(1,2)))*u1(1,2)*u2(1,2)*u3(1,2)];
	case 4
		u1 = varargin{1}; u2 = varargin{2}; u3 = varargin{3}; u4 = varargin{4};
		y = [(6144*par(1)/(1+exp(-4*xx(1,2)))^5*exp(-4*xx(1,2))^4-9216*par(1)/(1+exp(-4*xx(1,2)))^4*exp(-4*xx(1,2))^3+3584*par(1)/(1+exp(-4*xx(1,2)))^3*exp(-4*xx(1,2))^2-256*par(1)/(1+exp(-4*xx(1,2)))^2*exp(-4*xx(1,2)))*u1(1,2)*u2(1,2)*u3(1,2)*u4(1,2);
			(6144*par(3)/(1+exp(-4*xx(1,2)))^5*exp(-4*xx(1,2))^4-9216*par(3)/(1+exp(-4*xx(1,2)))^4*exp(-4*xx(1,2))^3+3584*par(3)/(1+exp(-4*xx(1,2)))^3*exp(-4*xx(1,2))^2-256*par(3)/(1+exp(-4*xx(1,2)))^2*exp(-4*xx(1,2)))*u1(1,2)*u2(1,2)*u3(1,2)*u4(1,2)];
	case 5
		u1 = varargin{1}; u2 = varargin{2}; u3 = varargin{3}; u4 = varargin{4}; u5 = varargin{5};
		y = [(122880*par(1)/(1+exp(-4*xx(1,2)))^6*exp(-4*xx(1,2))^5-245760*par(1)/(1+exp(-4*xx(1,2)))^5*exp(-4*xx(1,2))^4+153600*par(1)/(1+exp(-4*xx(1,2)))^4*exp(-4*xx(1,2))^3-30720*par(1)/(1+exp(-4*xx(1,2)))^3*exp(-4*xx(1,2))^2+1024*par(1)/(1+exp(-4*xx(1,2)))^2*exp(-4*xx(1,2)))*u1(1,2)*u2(1,2)*u3(1,2)*u4(1,2)*u5(1,2);
			(122880*par(3)/(1+exp(-4*xx(1,2)))^6*exp(-4*xx(1,2))^5-245760*par(3)/(1+exp(-4*xx(1,2)))^5*exp(-4*xx(1,2))^4+153600*par(3)/(1+exp(-4*xx(1,2)))^4*exp(-4*xx(1,2))^3-30720*par(3)/(1+exp(-4*xx(1,2)))^3*exp(-4*xx(1,2))^2+1024*par(3)/(1+exp(-4*xx(1,2)))^2*exp(-4*xx(1,2)))*u1(1,2)*u2(1,2)*u3(1,2)*u4(1,2)*u5(1,2)];
end