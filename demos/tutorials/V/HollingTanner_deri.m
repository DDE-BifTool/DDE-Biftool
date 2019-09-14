function J = HollingTanner_deri(xx,par,nx,np,v)
%%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: HollingTanner_deri.m 134 2016-09-12 11:10:44Z mmbosschaert $
%%

J = [];

if length(nx) == 1 && isempty(np) && isempty(v)
	switch nx
	case 0
		J = [1-2.*xx(1,1,:)-2.*par(4)-xx(2,1,:)./(par(3).*xx(2,1,:)+xx(1,1,:))+xx(1,1,:).*xx(2,1,:)./(par(3).*xx(2,1,:)+xx(1,1,:)).^2, -xx(1,1,:)./(par(3).*xx(2,1,:)+xx(1,1,:))+xx(1,1,:).*xx(2,1,:)./(par(3).*xx(2,1,:)+xx(1,1,:)).^2.*par(3); zeros(size(xx(1,1,:))), par(6).*(par(1)-xx(2,2,:)./xx(1,2,:))];
	case 1
		J = [zeros(size(xx(1,1,:))), zeros(size(xx(1,1,:))); par(6).*xx(2,1,:).*xx(2,2,:)./xx(1,2,:).^2, -par(6).*xx(2,1,:)./xx(1,2,:)];
	end
elseif isempty(nx) && length(np) == 1 && isempty(v)
	switch np
	case 1
		J = [zeros(size(xx(1,1,:))); par(6).*xx(2,1,:)];
	case 2
		J = [zeros(size(xx(1,1,:))); zeros(size(xx(1,1,:)))];
	case 3
		J = [xx(1,1,:).*xx(2,1,:).^2./(par(3).*xx(2,1,:)+xx(1,1,:)).^2; zeros(size(xx(1,1,:)))];
	case 4
		J = [1-2.*xx(1,1,:)-2.*par(4); zeros(size(xx(1,1,:)))];
	case 5
		J = [-ones(size(xx(1,1,:))); zeros(size(xx(1,1,:)))];
	case 6
		J = [zeros(size(xx(1,1,:))); xx(2,1,:).*(par(1)-xx(2,2,:)./xx(1,2,:))];
	end
elseif length(nx) == 1 && length(np) == 1 && isempty(v)
	switch nx
	case 0
		switch np
		case 1
			J = [zeros(size(xx(1,1,:))), zeros(size(xx(1,1,:))); zeros(size(xx(1,1,:))), par(6)*ones(size(xx(1,1,:)))];
		case 2
			J = [zeros(size(xx(1,1,:))), zeros(size(xx(1,1,:))); zeros(size(xx(1,1,:))), zeros(size(xx(1,1,:)))];
		case 3
			J = [xx(2,1,:).^2./(par(3).*xx(2,1,:)+xx(1,1,:)).^2-2.*xx(1,1,:).*xx(2,1,:).^2./(par(3).*xx(2,1,:)+xx(1,1,:)).^3, 2.*xx(1,1,:).*xx(2,1,:)./(par(3).*xx(2,1,:)+xx(1,1,:)).^2-2.*xx(1,1,:).*xx(2,1,:).^2./(par(3).*xx(2,1,:)+xx(1,1,:)).^3.*par(3); zeros(size(xx(1,1,:))), zeros(size(xx(1,1,:)))];
		case 4
			J = [-2, zeros(size(xx(1,1,:))); zeros(size(xx(1,1,:))), zeros(size(xx(1,1,:)))];
		case 5
			J = [zeros(size(xx(1,1,:))), zeros(size(xx(1,1,:))); zeros(size(xx(1,1,:))), zeros(size(xx(1,1,:)))];
		case 6
			J = [zeros(size(xx(1,1,:))), zeros(size(xx(1,1,:))); zeros(size(xx(1,1,:))), par(1)-xx(2,2,:)./xx(1,2,:)];
		end
	case 1
		switch np
		case 1
			J = [zeros(size(xx(1,1,:))), zeros(size(xx(1,1,:))); zeros(size(xx(1,1,:))), zeros(size(xx(1,1,:)))];
		case 2
			J = [zeros(size(xx(1,1,:))), zeros(size(xx(1,1,:))); zeros(size(xx(1,1,:))), zeros(size(xx(1,1,:)))];
		case 3
			J = [zeros(size(xx(1,1,:))), zeros(size(xx(1,1,:))); zeros(size(xx(1,1,:))), zeros(size(xx(1,1,:)))];
		case 4
			J = [zeros(size(xx(1,1,:))), zeros(size(xx(1,1,:))); zeros(size(xx(1,1,:))), zeros(size(xx(1,1,:)))];
		case 5
			J = [zeros(size(xx(1,1,:))), zeros(size(xx(1,1,:))); zeros(size(xx(1,1,:))), zeros(size(xx(1,1,:)))];
		case 6
			J = [zeros(size(xx(1,1,:))), zeros(size(xx(1,1,:))); xx(2,1,:).*xx(2,2,:)./xx(1,2,:).^2, -xx(2,1,:)./xx(1,2,:)];
		end
	end
elseif length(nx) == 2 && isempty(np) && ~isempty(v)
	nx1 = nx(1); nx2 = nx(2);
	switch nx1
	case 0
		switch nx2
		case 0
			J = [(-2+2.*xx(2,1,:)./(par(3).*xx(2,1,:)+xx(1,1,:)).^2-2.*xx(1,1,:).*xx(2,1,:)./(par(3).*xx(2,1,:)+xx(1,1,:)).^3).*v(1,:,:)+(-1./(par(3).*xx(2,1,:)+xx(1,1,:))+xx(2,1,:)./(par(3).*xx(2,1,:)+xx(1,1,:)).^2.*par(3)+xx(1,1,:)./(par(3).*xx(2,1,:)+xx(1,1,:)).^2-2.*xx(1,1,:).*xx(2,1,:)./(par(3).*xx(2,1,:)+xx(1,1,:)).^3.*par(3)).*v(2,:,:), (-1./(par(3).*xx(2,1,:)+xx(1,1,:))+xx(2,1,:)./(par(3).*xx(2,1,:)+xx(1,1,:)).^2.*par(3)+xx(1,1,:)./(par(3).*xx(2,1,:)+xx(1,1,:)).^2-2.*xx(1,1,:).*xx(2,1,:)./(par(3).*xx(2,1,:)+xx(1,1,:)).^3.*par(3)).*v(1,:,:)+(2.*xx(1,1,:)./(par(3).*xx(2,1,:)+xx(1,1,:)).^2.*par(3)-2.*xx(1,1,:).*xx(2,1,:)./(par(3).*xx(2,1,:)+xx(1,1,:)).^3.*par(3).^2).*v(2,:,:); zeros(size(xx(1,1,:))), zeros(size(xx(1,1,:)))];
		case 1
			J = [zeros(size(xx(1,1,:))), zeros(size(xx(1,1,:))); par(6).*xx(2,2,:)./xx(1,2,:).^2.*v(2,:,:), -par(6)./xx(1,2,:).*v(2,:,:)];
		end
	case 1
		switch nx2
		case 0
			J = [zeros(size(xx(1,1,:))), zeros(size(xx(1,1,:))); zeros(size(xx(1,1,:))), par(6).*xx(2,2,:)./xx(1,2,:).^2.*v(1,:,:)-par(6)./xx(1,2,:).*v(2,:,:)];
		case 1
			J = [zeros(size(xx(1,1,:))), zeros(size(xx(1,1,:))); -2.*par(6).*xx(2,1,:).*xx(2,2,:)./xx(1,2,:).^3.*v(1,:,:)+par(6).*xx(2,1,:)./xx(1,2,:).^2.*v(2,:,:), par(6).*xx(2,1,:)./xx(1,2,:).^2.*v(1,:,:)];
		end
	end
end
if isempty(J)
	display([nx np size(v)]);
	error('SYS_DERI: requested derivative could not be computed!');
end
