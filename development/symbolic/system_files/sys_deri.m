function J = sys_deri(xx,par,nx,np,v)

J = [];

if length(nx) == 1 && isempty(np) && isempty(v)
	switch nx
	case 0
		J = [-3*par(3)*xx(1,1)^2,1,-par(5);-2*par(6)*xx(1,1),-1,0;par(2)*par(8),0,-par(8)];
	case 1
		J = [2*par(4)*xx(1,2),0,0;0,0,0;0,0,0];
	end
elseif isempty(nx) && length(np) == 1 && isempty(v)
	switch np
	case 1
		J = [1;0;0];
	case 2
		J = [0;0;-par(8)*(par(7) - xx(1,1))];
	end
elseif length(nx) == 1 && length(np) == 1 && isempty(v)
	switch nx
	case 0
		switch np
		case 1
			J = [0,0,0;0,0,0;0,0,0];
		case 2
			J = [0,0,0;0,0,0;par(8),0,0];
		end
	case 1
		switch np
		case 1
			J = [0,0,0;0,0,0;0,0,0];
		case 2
			J = [0,0,0;0,0,0;0,0,0];
		end
	end
elseif length(nx) == 2 && isempty(np) && ~isempty(v)
	nx1 = nx(1); nx2 = nx(2);
	switch nx1
	case 0
		switch nx2
		case 0
			J = [-6*par(3)*v(1)*xx(1,1),0,0;-2*par(6)*v(1),0,0;0,0,0];
		case 1
			J = [0,0,0;0,0,0;0,0,0];
		end
	case 1
		switch nx2
		case 0
			J = [0,0,0;0,0,0;0,0,0];
		case 1
			J = [2*par(4)*v(1),0,0;0,0,0;0,0,0];
		end
	end
end
if isempty(J)
	display([nx np size(v)]);
	error('SYS_DERI: requested derivative could not be computed!');
end