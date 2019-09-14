function J = sys_deri_vec(xx,par,nx,np,v)

J = [];

I=ones(size(xx(1,1,:)));

if length(nx) == 1 && isempty(np) && isempty(v)
	switch nx
	case 0
		J = [(-3.*par(3).*xx(1,1,:).^2.*I),(1.*I),(-par(5).*I);(-2.*par(6).*xx(1,1,:).*I),(-1.*I),(0.*I);(par(2).*par(8).*I),(0.*I),(-par(8).*I)];
	case 1
		J = [(2.*par(4).*xx(1,2,:).*I),(0.*I),(0.*I);(0.*I),(0.*I),(0.*I);(0.*I),(0.*I),(0.*I)];
	end
elseif isempty(nx) && length(np) == 1 && isempty(v)
	switch np
	case 1
		J = [(1.*I);(0.*I);(0.*I)];
	case 2
		J = [(0.*I);(0.*I);(-par(8).*(par(7) - xx(1,1,:)).*I)];
	end
elseif length(nx) == 1 && length(np) == 1 && isempty(v)
	switch nx
	case 0
		switch np
		case 1
			J = [(0.*I),(0.*I),(0.*I);(0.*I),(0.*I),(0.*I);(0.*I),(0.*I),(0.*I)];
		case 2
			J = [(0.*I),(0.*I),(0.*I);(0.*I),(0.*I),(0.*I);(par(8).*I),(0.*I),(0.*I)];
		end
	case 1
		switch np
		case 1
			J = [(0.*I),(0.*I),(0.*I);(0.*I),(0.*I),(0.*I);(0.*I),(0.*I),(0.*I)];
		case 2
			J = [(0.*I),(0.*I),(0.*I);(0.*I),(0.*I),(0.*I);(0.*I),(0.*I),(0.*I)];
		end
	end
elseif length(nx) == 2 && isempty(np) && ~isempty(v)
	nx1 = nx(1); nx2 = nx(2);
	switch nx1
	case 0
		switch nx2
		case 0
			J = [(-6.*par(3).*v(1,:,:).*xx(1,1,:).*I),(0.*I),(0.*I);(-2.*par(6).*v(1,:,:).*I),(0.*I),(0.*I);(0.*I),(0.*I),(0.*I)];
		case 1
			J = [(0.*I),(0.*I),(0.*I);(0.*I),(0.*I),(0.*I);(0.*I),(0.*I),(0.*I)];
		end
	case 1
		switch nx2
		case 0
			J = [(0.*I),(0.*I),(0.*I);(0.*I),(0.*I),(0.*I);(0.*I),(0.*I),(0.*I)];
		case 1
			J = [(2.*par(4).*v(1,:,:).*I),(0.*I),(0.*I);(0.*I),(0.*I),(0.*I);(0.*I),(0.*I),(0.*I)];
		end
	end
end
if isempty(J)
	display([nx np size(v)]);
	error('SYS_DERI: requested derivative could not be computed!');
end