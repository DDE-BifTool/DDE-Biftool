function J = neural_deri(xx,par,nx,np,v)
%%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: neural_deri.m 134 2016-09-12 11:10:44Z mmbosschaert $
%%

J = [];

if length(nx) == 1 && isempty(np) && isempty(v)
	switch nx
	case 0
		J = [-1, 0; 0, -1];
	case 1
		J = [-par(1)*par(2)*(1-tanh(par(2)*xx(1,2)-1)^2)*cosh(1)^2, 0; 0, -par(1)*par(2)*(1-tanh(par(2)*xx(2,2)-1)^2)*cosh(1)^2];
	case 2
		J = [0, par(3)*par(4)*(1-tanh(par(4)*xx(2,3)-1)^2)*cosh(1)^2; par(3)*par(4)*(1-tanh(par(4)*xx(1,3)-1)^2)*cosh(1)^2, 0];
	end
elseif isempty(nx) && length(np) == 1 && isempty(v)
	switch np
	case 1
		J = [-(tanh(par(2)*xx(1,2)-1)+tanh(1))*cosh(1)^2; -(tanh(par(2)*xx(2,2)-1)+tanh(1))*cosh(1)^2];
	case 2
		J = [-par(1)*xx(1,2)*(1-tanh(par(2)*xx(1,2)-1)^2)*cosh(1)^2; -par(1)*xx(2,2)*(1-tanh(par(2)*xx(2,2)-1)^2)*cosh(1)^2];
	case 3
		J = [(tanh(par(4)*xx(2,3)-1)+tanh(1))*cosh(1)^2; (tanh(par(4)*xx(1,3)-1)+tanh(1))*cosh(1)^2];
	case 4
		J = [par(3)*xx(2,3)*(1-tanh(par(4)*xx(2,3)-1)^2)*cosh(1)^2; par(3)*xx(1,3)*(1-tanh(par(4)*xx(1,3)-1)^2)*cosh(1)^2];
	case 5
		J = [0; 0];
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
			J = [-par(2)*(1-tanh(par(2)*xx(1,2)-1)^2)*cosh(1)^2, 0; 0, -par(2)*(1-tanh(par(2)*xx(2,2)-1)^2)*cosh(1)^2];
		case 2
			J = [-par(1)*(1-tanh(par(2)*xx(1,2)-1)^2)*cosh(1)^2+2*par(1)*xx(1,2)*tanh(par(2)*xx(1,2)-1)*par(2)*(1-tanh(par(2)*xx(1,2)-1)^2)*cosh(1)^2, 0; 0, -par(1)*(1-tanh(par(2)*xx(2,2)-1)^2)*cosh(1)^2+2*par(1)*xx(2,2)*tanh(par(2)*xx(2,2)-1)*par(2)*(1-tanh(par(2)*xx(2,2)-1)^2)*cosh(1)^2];
		case 3
			J = [0, 0; 0, 0];
		case 4
			J = [0, 0; 0, 0];
		case 5
			J = [0, 0; 0, 0];
		case 6
			J = [0, 0; 0, 0];
		end
	case 2
		switch np
		case 1
			J = [0, 0; 0, 0];
		case 2
			J = [0, 0; 0, 0];
		case 3
			J = [0, par(4)*(1-tanh(par(4)*xx(2,3)-1)^2)*cosh(1)^2; par(4)*(1-tanh(par(4)*xx(1,3)-1)^2)*cosh(1)^2, 0];
		case 4
			J = [0, par(3)*(1-tanh(par(4)*xx(2,3)-1)^2)*cosh(1)^2-2*par(3)*xx(2,3)*tanh(par(4)*xx(2,3)-1)*par(4)*(1-tanh(par(4)*xx(2,3)-1)^2)*cosh(1)^2; par(3)*(1-tanh(par(4)*xx(1,3)-1)^2)*cosh(1)^2-2*par(3)*xx(1,3)*tanh(par(4)*xx(1,3)-1)*par(4)*(1-tanh(par(4)*xx(1,3)-1)^2)*cosh(1)^2, 0];
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
		case 2
			J = [0, 0; 0, 0];
		end
	case 1
		switch nx2
		case 0
			J = [0, 0; 0, 0];
		case 1
			J = [2*par(1)*par(2)^2*tanh(par(2)*xx(1,2)-1)*(1-tanh(par(2)*xx(1,2)-1)^2)*cosh(1)^2*v(1), 0; 0, 2*par(1)*par(2)^2*tanh(par(2)*xx(2,2)-1)*(1-tanh(par(2)*xx(2,2)-1)^2)*cosh(1)^2*v(2)];
		case 2
			J = [0, 0; 0, 0];
		end
	case 2
		switch nx2
		case 0
			J = [0, 0; 0, 0];
		case 1
			J = [0, 0; 0, 0];
		case 2
			J = [0, -2*par(3)*par(4)^2*tanh(par(4)*xx(2,3)-1)*(1-tanh(par(4)*xx(2,3)-1)^2)*cosh(1)^2*v(2); -2*par(3)*par(4)^2*tanh(par(4)*xx(1,3)-1)*(1-tanh(par(4)*xx(1,3)-1)^2)*cosh(1)^2*v(1), 0];
		end
	end
end
if isempty(J)
	display([nx np size(v)]);
	error('SYS_DERI: requested derivative could not be computed!');
end
