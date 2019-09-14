function J = cusp_deri(xx,par,nx,np,v)
%% first and 2nd order derivatives in xx and par for cusp demo (see cusp_demo.m for equations)
%
% $Id: cusp_deri.m 121 2015-09-02 22:02:24Z jansieber $
%
%%
J = [];

if length(nx) == 1 && isempty(np) && isempty(v)
	switch nx
	case 0
		J = [-1, 0; 0, -1];
	case 1
		J = [4*par(1)/(1+exp(-4*xx(1,2)))^2*exp(-4*xx(1,2)), -par(2); 4*par(3)/(1+exp(-4*xx(1,2)))^2*exp(-4*xx(1,2)), 0];
	end
elseif isempty(nx) && length(np) == 1 && isempty(v)
	switch np
	case 1
		J = [1/(1+exp(-4*xx(1,2))); 0];
	case 2
		J = [-xx(2,2); 0];
	case 3
		J = [0; 1/(1+exp(-4*xx(1,2)))];
	case 4
		J = [1; 0];
	case 5
		J = [0; 1];
	case 6
		J = [0; 0];
	end
elseif length(nx) == 1 && length(np) == 1 && isempty(v)
	switch nx
	case 0
		switch np
		case 1
			J = [0, 0; 0, 0];
		case 2
			J = [0, 0; 0, 0];
		case 3
			J = [0, 0; 0, 0];
		case 4
			J = [0, 0; 0, 0];
		case 5
			J = [0, 0; 0, 0];
		case 6
			J = [0, 0; 0, 0];
		end
	case 1
		switch np
		case 1
			J = [4/(1+exp(-4*xx(1,2)))^2*exp(-4*xx(1,2)), 0; 0, 0];
		case 2
			J = [0, -1; 0, 0];
		case 3
			J = [0, 0; 4/(1+exp(-4*xx(1,2)))^2*exp(-4*xx(1,2)), 0];
		case 4
			J = [0, 0; 0, 0];
		case 5
			J = [0, 0; 0, 0];
		case 6
			J = [0, 0; 0, 0];
		end
	end
elseif length(nx) == 2 && isempty(np) && ~isempty(v)
	nx1 = nx(1); nx2 = nx(2);
	switch nx1
	case 0
		switch nx2
		case 0
			J = [0, 0; 0, 0];
		case 1
			J = [0, 0; 0, 0];
		end
	case 1
		switch nx2
		case 0
			J = [0, 0; 0, 0];
		case 1
			J = [32*par(1)/(1+exp(-4*xx(1,2)))^3*exp(-4*xx(1,2))^2*v(1)-16*par(1)/(1+exp(-4*xx(1,2)))^2*exp(-4*xx(1,2))*v(1), 0; 32*par(3)/(1+exp(-4*xx(1,2)))^3*exp(-4*xx(1,2))^2*v(1)-16*par(3)/(1+exp(-4*xx(1,2)))^2*exp(-4*xx(1,2))*v(1), 0];
		end
	end
end
if isempty(J)
	display([nx np size(v)]);
	error('SYS_DERI: requested derivative could not be computed!');
end