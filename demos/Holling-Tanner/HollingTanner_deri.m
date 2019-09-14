function J = HollingTanner_deri(xx,par,nx,np,v)
%% Derivatives wrt states and parameters
%
% $Id: HollingTanner_deri.m 109 2015-08-31 23:45:11Z jansieber $
%

J = [];

if length(nx) == 1 && isempty(np) && isempty(v)
	switch nx
	case 0
		J = [1-2*xx(1,1)-2*par(4)-xx(2,1)/(par(3)*xx(2,1)+xx(1,1))+xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^2, -xx(1,1)/(par(3)*xx(2,1)+xx(1,1))+xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^2*par(3); 0, par(6)*(par(1)-xx(2,2)/xx(1,2))];
	case 1
		J = [0, 0; par(6)*xx(2,1)*xx(2,2)/xx(1,2)^2, -par(6)*xx(2,1)/xx(1,2)];
	end
elseif isempty(nx) && length(np) == 1 && isempty(v)
	switch np
	case 1
		J = [0; par(6)*xx(2,1)];
	case 2
		J = [0; 0];
	case 3
		J = [xx(1,1)*xx(2,1)^2/(par(3)*xx(2,1)+xx(1,1))^2; 0];
	case 4
		J = [1-2*xx(1,1)-2*par(4); 0];
	case 5
		J = [-1; 0];
	case 6
		J = [0; xx(2,1)*(par(1)-xx(2,2)/xx(1,2))];
	end
elseif length(nx) == 1 && length(np) == 1 && isempty(v)
	switch nx
	case 0
		switch np
		case 1
			J = [0, 0; 0, par(6)];
		case 2
			J = [0, 0; 0, 0];
		case 3
			J = [xx(2,1)^2/(par(3)*xx(2,1)+xx(1,1))^2-2*xx(1,1)*xx(2,1)^2/(par(3)*xx(2,1)+xx(1,1))^3, 2*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^2-2*xx(1,1)*xx(2,1)^2/(par(3)*xx(2,1)+xx(1,1))^3*par(3); 0, 0];
		case 4
			J = [-2, 0; 0, 0];
		case 5
			J = [0, 0; 0, 0];
		case 6
			J = [0, 0; 0, par(1)-xx(2,2)/xx(1,2)];
		end
	case 1
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
			J = [0, 0; xx(2,1)*xx(2,2)/xx(1,2)^2, -xx(2,1)/xx(1,2)];
		end
	end
elseif length(nx) == 2 && isempty(np) && ~isempty(v)
	nx1 = nx(1); nx2 = nx(2);
	switch nx1
	case 0
		switch nx2
		case 0
			J = [(-2+2*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^2-2*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^3)*v(1)+(-1/(par(3)*xx(2,1)+xx(1,1))+xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^2*par(3)+xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^2-2*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^3*par(3))*v(2), (-1/(par(3)*xx(2,1)+xx(1,1))+xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^2*par(3)+xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^2-2*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^3*par(3))*v(1)+(2*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^2*par(3)-2*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^3*par(3)^2)*v(2); 0, 0];
		case 1
			J = [0, 0; par(6)*xx(2,2)/xx(1,2)^2*v(2), -par(6)/xx(1,2)*v(2)];
		end
	case 1
		switch nx2
		case 0
			J = [0, 0; 0, par(6)*xx(2,2)/xx(1,2)^2*v(1)-par(6)/xx(1,2)*v(2)];
		case 1
			J = [0, 0; -2*par(6)*xx(2,1)*xx(2,2)/xx(1,2)^3*v(1)+par(6)*xx(2,1)/xx(1,2)^2*v(2), par(6)*xx(2,1)/xx(1,2)^2*v(1)];
		end
	end
end
if isempty(J)
	display([nx np size(v)]);
	error('SYS_DERI: requested derivative could not be computed!');
end