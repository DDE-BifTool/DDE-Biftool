function y = sys_mfderi(xx,par,varargin)

if nargin == 2
	error('SYS_MFDERI: no arguments.');
elseif nargin > 7
	error('SYS_MFDERI: too many arguments.');
end

y=0;

numarg = nargin - 2;

switch numarg
	case 2
		u1 = varargin{1}; u2 = varargin{2};
		y=sys_B(xx,par,u1(:),u2(:));
	case 3
		u1 = varargin{1}; u2 = varargin{2}; u3 = varargin{3};
		y=sys_C(xx,par,u1(:),u2(:),u3(:));
	case 4
		u1 = varargin{1}; u2 = varargin{2}; u3 = varargin{3}; u4 = varargin{4};
		y=sys_D(xx,par,u1(:),u2(:),u3(:),u4(:));
	case 5
		u1 = varargin{1}; u2 = varargin{2}; u3 = varargin{3}; u4 = varargin{4}; u5 = varargin{5};
		y=sys_E(xx,par,u1(:),u2(:),u3(:),u4(:),u5(:));

end