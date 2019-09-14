function J = sys_deri(xx,par,nx,np,v)
%%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: sys_deri.m 134 2016-09-12 11:10:44Z mmbosschaert $
%%
J = [];

if length(nx) == 1 && isempty(np) && isempty(v)
	switch nx
	case 0
		J = [-par(4)-par(1).*par(3).*(par(1).*xx(1,4,:)+par(2).*xx(1,5,:)+par(4).*xx(1,2,:))-par(2).*par(3).*(par(1).*xx(1,5,:)+par(2).*xx(1,6,:)+par(4).*xx(1,3,:))-par(1).^2.*par(3).^2.*xx(1,2,:).*(par(1).*xx(1,7,:)+par(2).*xx(1,8,:)+par(4).*xx(1,4,:))-par(1).*par(2).*par(3).^2.*xx(1,2,:).*(par(1).*xx(1,8,:)+par(2).*xx(1,9,:)+par(4).*xx(1,5,:))-par(2).*par(1).*par(3).^2.*xx(1,3,:).*(par(1).*xx(1,8,:)+par(2).*xx(1,9,:)+par(4).*xx(1,5,:))-par(2).^2.*par(3).^2.*xx(1,3,:).*(par(1).*xx(1,9,:)+par(2).*xx(1,10,:)+par(4).*xx(1,6,:))-par(3).^2.*xx(1,1,:).*par(1).*(par(4).^2.*xx(1,2,:)+2.*par(4).*(par(1).*xx(1,4,:)+par(2).*xx(1,5,:))+par(1).^2.*xx(1,7,:)+2.*par(1).*par(2).*xx(1,8,:)+par(2).^2.*xx(1,9,:))-par(3).^2.*xx(1,1,:).*par(2).*(par(4).^2.*xx(1,3,:)+2.*par(4).*(par(1).*xx(1,5,:)+par(2).*xx(1,6,:))+par(1).^2.*xx(1,8,:)+2.*par(1).*par(2).*xx(1,9,:)+par(2).^2.*xx(1,10,:))];
	case 1
		J = [-par(1)-par(1).*par(3).*xx(1,1,:).*par(4)-par(1).^2.*par(3).^2.*xx(1,1,:).*(par(1).*xx(1,7,:)+par(2).*xx(1,8,:)+par(4).*xx(1,4,:))-par(1).*par(2).*par(3).^2.*xx(1,1,:).*(par(1).*xx(1,8,:)+par(2).*xx(1,9,:)+par(4).*xx(1,5,:))-1/2.*par(3).^2.*xx(1,1,:).^2.*par(1).*par(4).^2];
	case 2
		J = [-par(2)-par(2).*par(3).*xx(1,1,:).*par(4)-par(1).*par(2).*par(3).^2.*xx(1,1,:).*(par(1).*xx(1,8,:)+par(2).*xx(1,9,:)+par(4).*xx(1,5,:))-par(2).^2.*par(3).^2.*xx(1,1,:).*(par(1).*xx(1,9,:)+par(2).*xx(1,10,:)+par(4).*xx(1,6,:))-1/2.*par(3).^2.*xx(1,1,:).^2.*par(2).*par(4).^2];
	case 3
		J = [-par(1).^2.*par(3).^2.*par(4).*xx(1,1,:).^2-par(1).^2.*par(3).^2.*par(4).*xx(1,1,:).*xx(1,2,:)-par(1).^2.*par(3).*xx(1,1,:)];
	case 4
		J = [-2.*par(1).*par(2).*par(3).^2.*par(4).*xx(1,1,:).^2-par(1).*par(2).*par(3).^2.*par(4).*xx(1,1,:).*xx(1,2,:)-par(1).*par(2).*par(3).^2.*par(4).*xx(1,1,:).*xx(1,3,:)-2.*par(1).*par(2).*par(3).*xx(1,1,:)];
	case 5
		J = [-par(2).^2.*par(3).^2.*par(4).*xx(1,1,:).^2-par(2).^2.*par(3).^2.*par(4).*xx(1,1,:).*xx(1,3,:)-par(2).^2.*par(3).*xx(1,1,:)];
	case 6
		J = [-par(1).^3.*par(3).^2.*xx(1,1,:).*xx(1,2,:)-1/2.*par(3).^2.*xx(1,1,:).^2.*par(1).^3];
	case 7
		J = [-2.*par(1).^2.*par(3).^2.*xx(1,1,:).*xx(1,2,:).*par(2)-par(2).*par(1).^2.*par(3).^2.*xx(1,1,:).*xx(1,3,:)-3/2.*par(3).^2.*xx(1,1,:).^2.*par(1).^2.*par(2)];
	case 8
		J = [-par(1).*par(2).^2.*par(3).^2.*xx(1,1,:).*xx(1,2,:)-2.*par(2).^2.*par(1).*par(3).^2.*xx(1,1,:).*xx(1,3,:)-3/2.*par(3).^2.*xx(1,1,:).^2.*par(1).*par(2).^2];
	case 9
		J = [-par(2).^3.*par(3).^2.*xx(1,1,:).*xx(1,3,:)-1/2.*par(3).^2.*xx(1,1,:).^2.*par(2).^3];
	end
elseif isempty(nx) && length(np) == 1 && isempty(v)
	switch np
	case 1
		J = [-xx(1,2,:)-par(3).*xx(1,1,:).*(par(1).*xx(1,4,:)+par(2).*xx(1,5,:)+par(4).*xx(1,2,:))-par(1).*par(3).*xx(1,1,:).*xx(1,4,:)-par(2).*par(3).*xx(1,1,:).*xx(1,5,:)-2.*par(1).*par(3).^2.*xx(1,1,:).*xx(1,2,:).*(par(1).*xx(1,7,:)+par(2).*xx(1,8,:)+par(4).*xx(1,4,:))-par(1).^2.*par(3).^2.*xx(1,1,:).*xx(1,2,:).*xx(1,7,:)-par(2).*par(3).^2.*xx(1,1,:).*xx(1,2,:).*(par(1).*xx(1,8,:)+par(2).*xx(1,9,:)+par(4).*xx(1,5,:))-par(1).*par(2).*par(3).^2.*xx(1,1,:).*xx(1,2,:).*xx(1,8,:)-par(2).*par(3).^2.*xx(1,1,:).*xx(1,3,:).*(par(1).*xx(1,8,:)+par(2).*xx(1,9,:)+par(4).*xx(1,5,:))-par(2).*par(1).*par(3).^2.*xx(1,1,:).*xx(1,3,:).*xx(1,8,:)-par(2).^2.*par(3).^2.*xx(1,1,:).*xx(1,3,:).*xx(1,9,:)-1/2.*par(3).^2.*xx(1,1,:).^2.*(par(4).^2.*xx(1,2,:)+2.*par(4).*(par(1).*xx(1,4,:)+par(2).*xx(1,5,:))+par(1).^2.*xx(1,7,:)+2.*par(1).*par(2).*xx(1,8,:)+par(2).^2.*xx(1,9,:))-1/2.*par(3).^2.*xx(1,1,:).^2.*par(1).*(2.*par(1).*xx(1,7,:)+2.*par(2).*xx(1,8,:)+2.*par(4).*xx(1,4,:))-1/2.*par(3).^2.*xx(1,1,:).^2.*par(2).*(2.*par(1).*xx(1,8,:)+2.*par(2).*xx(1,9,:)+2.*par(4).*xx(1,5,:))];
	case 2
		J = [-xx(1,3,:)-par(1).*par(3).*xx(1,1,:).*xx(1,5,:)-par(3).*xx(1,1,:).*(par(1).*xx(1,5,:)+par(2).*xx(1,6,:)+par(4).*xx(1,3,:))-par(2).*par(3).*xx(1,1,:).*xx(1,6,:)-par(1).^2.*par(3).^2.*xx(1,1,:).*xx(1,2,:).*xx(1,8,:)-par(1).*par(3).^2.*xx(1,1,:).*xx(1,2,:).*(par(1).*xx(1,8,:)+par(2).*xx(1,9,:)+par(4).*xx(1,5,:))-par(1).*par(2).*par(3).^2.*xx(1,1,:).*xx(1,2,:).*xx(1,9,:)-par(1).*par(3).^2.*xx(1,1,:).*xx(1,3,:).*(par(1).*xx(1,8,:)+par(2).*xx(1,9,:)+par(4).*xx(1,5,:))-par(2).*par(1).*par(3).^2.*xx(1,1,:).*xx(1,3,:).*xx(1,9,:)-2.*par(2).*par(3).^2.*xx(1,1,:).*xx(1,3,:).*(par(1).*xx(1,9,:)+par(2).*xx(1,10,:)+par(4).*xx(1,6,:))-par(2).^2.*par(3).^2.*xx(1,1,:).*xx(1,3,:).*xx(1,10,:)-1/2.*par(3).^2.*xx(1,1,:).^2.*par(1).*(2.*par(1).*xx(1,8,:)+2.*par(2).*xx(1,9,:)+2.*par(4).*xx(1,5,:))-1/2.*par(3).^2.*xx(1,1,:).^2.*(par(4).^2.*xx(1,3,:)+2.*par(4).*(par(1).*xx(1,5,:)+par(2).*xx(1,6,:))+par(1).^2.*xx(1,8,:)+2.*par(1).*par(2).*xx(1,9,:)+par(2).^2.*xx(1,10,:))-1/2.*par(3).^2.*xx(1,1,:).^2.*par(2).*(2.*par(1).*xx(1,9,:)+2.*par(2).*xx(1,10,:)+2.*par(4).*xx(1,6,:))];
	case 3
		J = [-par(1).*xx(1,1,:).*(par(1).*xx(1,4,:)+par(2).*xx(1,5,:)+par(4).*xx(1,2,:))-par(2).*xx(1,1,:).*(par(1).*xx(1,5,:)+par(2).*xx(1,6,:)+par(4).*xx(1,3,:))-2.*par(1).^2.*par(3).*xx(1,1,:).*xx(1,2,:).*(par(1).*xx(1,7,:)+par(2).*xx(1,8,:)+par(4).*xx(1,4,:))-2.*par(1).*par(2).*par(3).*xx(1,1,:).*xx(1,2,:).*(par(1).*xx(1,8,:)+par(2).*xx(1,9,:)+par(4).*xx(1,5,:))-2.*par(2).*par(1).*par(3).*xx(1,1,:).*xx(1,3,:).*(par(1).*xx(1,8,:)+par(2).*xx(1,9,:)+par(4).*xx(1,5,:))-2.*par(2).^2.*par(3).*xx(1,1,:).*xx(1,3,:).*(par(1).*xx(1,9,:)+par(2).*xx(1,10,:)+par(4).*xx(1,6,:))-par(3).*xx(1,1,:).^2.*par(1).*(par(4).^2.*xx(1,2,:)+2.*par(4).*(par(1).*xx(1,4,:)+par(2).*xx(1,5,:))+par(1).^2.*xx(1,7,:)+2.*par(1).*par(2).*xx(1,8,:)+par(2).^2.*xx(1,9,:))-par(3).*xx(1,1,:).^2.*par(2).*(par(4).^2.*xx(1,3,:)+2.*par(4).*(par(1).*xx(1,5,:)+par(2).*xx(1,6,:))+par(1).^2.*xx(1,8,:)+2.*par(1).*par(2).*xx(1,9,:)+par(2).^2.*xx(1,10,:))];
	case 4
		J = [-xx(1,1,:)-par(1).*par(3).*xx(1,1,:).*xx(1,2,:)-par(2).*par(3).*xx(1,1,:).*xx(1,3,:)-par(1).^2.*par(3).^2.*xx(1,1,:).*xx(1,2,:).*xx(1,4,:)-par(1).*par(2).*par(3).^2.*xx(1,1,:).*xx(1,2,:).*xx(1,5,:)-par(2).*par(1).*par(3).^2.*xx(1,1,:).*xx(1,3,:).*xx(1,5,:)-par(2).^2.*par(3).^2.*xx(1,1,:).*xx(1,3,:).*xx(1,6,:)-1/2.*par(3).^2.*xx(1,1,:).^2.*par(1).*(2.*par(1).*xx(1,4,:)+2.*par(2).*xx(1,5,:)+2.*par(4).*xx(1,2,:))-1/2.*par(3).^2.*xx(1,1,:).^2.*par(2).*(2.*par(1).*xx(1,5,:)+2.*par(2).*xx(1,6,:)+2.*par(4).*xx(1,3,:))];
	end
elseif length(nx) == 1 && length(np) == 1 && isempty(v)
	switch nx
	case 0
		switch np
		case 1
			J = [-par(3).*(par(1).*xx(1,4,:)+par(2).*xx(1,5,:)+par(4).*xx(1,2,:))-par(1).*par(3).*xx(1,4,:)-par(2).*par(3).*xx(1,5,:)-2.*par(1).*par(3).^2.*xx(1,2,:).*(par(1).*xx(1,7,:)+par(2).*xx(1,8,:)+par(4).*xx(1,4,:))-par(1).^2.*par(3).^2.*xx(1,2,:).*xx(1,7,:)-par(2).*par(3).^2.*xx(1,2,:).*(par(1).*xx(1,8,:)+par(2).*xx(1,9,:)+par(4).*xx(1,5,:))-par(1).*par(2).*par(3).^2.*xx(1,2,:).*xx(1,8,:)-par(2).*par(3).^2.*xx(1,3,:).*(par(1).*xx(1,8,:)+par(2).*xx(1,9,:)+par(4).*xx(1,5,:))-par(2).*par(1).*par(3).^2.*xx(1,3,:).*xx(1,8,:)-par(2).^2.*par(3).^2.*xx(1,3,:).*xx(1,9,:)-par(3).^2.*xx(1,1,:).*(par(4).^2.*xx(1,2,:)+2.*par(4).*(par(1).*xx(1,4,:)+par(2).*xx(1,5,:))+par(1).^2.*xx(1,7,:)+2.*par(1).*par(2).*xx(1,8,:)+par(2).^2.*xx(1,9,:))-par(3).^2.*xx(1,1,:).*par(1).*(2.*par(1).*xx(1,7,:)+2.*par(2).*xx(1,8,:)+2.*par(4).*xx(1,4,:))-par(3).^2.*xx(1,1,:).*par(2).*(2.*par(1).*xx(1,8,:)+2.*par(2).*xx(1,9,:)+2.*par(4).*xx(1,5,:))];
		case 2
			J = [-par(1).*par(3).*xx(1,5,:)-par(3).*(par(1).*xx(1,5,:)+par(2).*xx(1,6,:)+par(4).*xx(1,3,:))-par(2).*par(3).*xx(1,6,:)-par(1).^2.*par(3).^2.*xx(1,2,:).*xx(1,8,:)-par(1).*par(3).^2.*xx(1,2,:).*(par(1).*xx(1,8,:)+par(2).*xx(1,9,:)+par(4).*xx(1,5,:))-par(1).*par(2).*par(3).^2.*xx(1,2,:).*xx(1,9,:)-par(1).*par(3).^2.*xx(1,3,:).*(par(1).*xx(1,8,:)+par(2).*xx(1,9,:)+par(4).*xx(1,5,:))-par(2).*par(1).*par(3).^2.*xx(1,3,:).*xx(1,9,:)-2.*par(2).*par(3).^2.*xx(1,3,:).*(par(1).*xx(1,9,:)+par(2).*xx(1,10,:)+par(4).*xx(1,6,:))-par(2).^2.*par(3).^2.*xx(1,3,:).*xx(1,10,:)-par(3).^2.*xx(1,1,:).*par(1).*(2.*par(1).*xx(1,8,:)+2.*par(2).*xx(1,9,:)+2.*par(4).*xx(1,5,:))-par(3).^2.*xx(1,1,:).*(par(4).^2.*xx(1,3,:)+2.*par(4).*(par(1).*xx(1,5,:)+par(2).*xx(1,6,:))+par(1).^2.*xx(1,8,:)+2.*par(1).*par(2).*xx(1,9,:)+par(2).^2.*xx(1,10,:))-par(3).^2.*xx(1,1,:).*par(2).*(2.*par(1).*xx(1,9,:)+2.*par(2).*xx(1,10,:)+2.*par(4).*xx(1,6,:))];
		case 3
			J = [-par(1).*(par(1).*xx(1,4,:)+par(2).*xx(1,5,:)+par(4).*xx(1,2,:))-par(2).*(par(1).*xx(1,5,:)+par(2).*xx(1,6,:)+par(4).*xx(1,3,:))-2.*par(1).^2.*par(3).*xx(1,2,:).*(par(1).*xx(1,7,:)+par(2).*xx(1,8,:)+par(4).*xx(1,4,:))-2.*par(1).*par(2).*par(3).*xx(1,2,:).*(par(1).*xx(1,8,:)+par(2).*xx(1,9,:)+par(4).*xx(1,5,:))-2.*par(2).*par(1).*par(3).*xx(1,3,:).*(par(1).*xx(1,8,:)+par(2).*xx(1,9,:)+par(4).*xx(1,5,:))-2.*par(2).^2.*par(3).*xx(1,3,:).*(par(1).*xx(1,9,:)+par(2).*xx(1,10,:)+par(4).*xx(1,6,:))-2.*par(3).*xx(1,1,:).*par(1).*(par(4).^2.*xx(1,2,:)+2.*par(4).*(par(1).*xx(1,4,:)+par(2).*xx(1,5,:))+par(1).^2.*xx(1,7,:)+2.*par(1).*par(2).*xx(1,8,:)+par(2).^2.*xx(1,9,:))-2.*par(3).*xx(1,1,:).*par(2).*(par(4).^2.*xx(1,3,:)+2.*par(4).*(par(1).*xx(1,5,:)+par(2).*xx(1,6,:))+par(1).^2.*xx(1,8,:)+2.*par(1).*par(2).*xx(1,9,:)+par(2).^2.*xx(1,10,:))];
		case 4
			J = [-1-par(1).*par(3).*xx(1,2,:)-par(2).*par(3).*xx(1,3,:)-par(1).^2.*par(3).^2.*xx(1,2,:).*xx(1,4,:)-par(1).*par(2).*par(3).^2.*xx(1,2,:).*xx(1,5,:)-par(2).*par(1).*par(3).^2.*xx(1,3,:).*xx(1,5,:)-par(2).^2.*par(3).^2.*xx(1,3,:).*xx(1,6,:)-par(3).^2.*xx(1,1,:).*par(1).*(2.*par(1).*xx(1,4,:)+2.*par(2).*xx(1,5,:)+2.*par(4).*xx(1,2,:))-par(3).^2.*xx(1,1,:).*par(2).*(2.*par(1).*xx(1,5,:)+2.*par(2).*xx(1,6,:)+2.*par(4).*xx(1,3,:))];
		end
	case 1
		switch np
		case 1
			J = [-1-par(3).*xx(1,1,:).*par(4)-2.*par(1).*par(3).^2.*xx(1,1,:).*(par(1).*xx(1,7,:)+par(2).*xx(1,8,:)+par(4).*xx(1,4,:))-par(1).^2.*par(3).^2.*xx(1,1,:).*xx(1,7,:)-par(2).*par(3).^2.*xx(1,1,:).*(par(1).*xx(1,8,:)+par(2).*xx(1,9,:)+par(4).*xx(1,5,:))-par(1).*par(2).*par(3).^2.*xx(1,1,:).*xx(1,8,:)-1/2.*par(3).^2.*xx(1,1,:).^2.*par(4).^2];
		case 2
			J = [-par(1).^2.*par(3).^2.*xx(1,1,:).*xx(1,8,:)-par(1).*par(3).^2.*xx(1,1,:).*(par(1).*xx(1,8,:)+par(2).*xx(1,9,:)+par(4).*xx(1,5,:))-par(1).*par(2).*par(3).^2.*xx(1,1,:).*xx(1,9,:)];
		case 3
			J = [-par(1).*xx(1,1,:).*par(4)-2.*par(1).^2.*par(3).*xx(1,1,:).*(par(1).*xx(1,7,:)+par(2).*xx(1,8,:)+par(4).*xx(1,4,:))-2.*par(1).*par(2).*par(3).*xx(1,1,:).*(par(1).*xx(1,8,:)+par(2).*xx(1,9,:)+par(4).*xx(1,5,:))-par(3).*xx(1,1,:).^2.*par(1).*par(4).^2];
		case 4
			J = [-par(1).^2.*par(3).^2.*xx(1,1,:).*xx(1,4,:)-par(1).*par(2).*par(3).^2.*xx(1,1,:).*xx(1,5,:)-par(1).*par(3).^2.*par(4).*xx(1,1,:).^2-par(1).*par(3).*xx(1,1,:)];
		end
	case 2
		switch np
		case 1
			J = [-par(2).*par(3).^2.*xx(1,1,:).*(par(1).*xx(1,8,:)+par(2).*xx(1,9,:)+par(4).*xx(1,5,:))-par(1).*par(2).*par(3).^2.*xx(1,1,:).*xx(1,8,:)-par(2).^2.*par(3).^2.*xx(1,1,:).*xx(1,9,:)];
		case 2
			J = [-1-par(3).*xx(1,1,:).*par(4)-par(1).*par(3).^2.*xx(1,1,:).*(par(1).*xx(1,8,:)+par(2).*xx(1,9,:)+par(4).*xx(1,5,:))-par(1).*par(2).*par(3).^2.*xx(1,1,:).*xx(1,9,:)-2.*par(2).*par(3).^2.*xx(1,1,:).*(par(1).*xx(1,9,:)+par(2).*xx(1,10,:)+par(4).*xx(1,6,:))-par(2).^2.*par(3).^2.*xx(1,1,:).*xx(1,10,:)-1/2.*par(3).^2.*xx(1,1,:).^2.*par(4).^2];
		case 3
			J = [-par(2).*xx(1,1,:).*par(4)-2.*par(1).*par(2).*par(3).*xx(1,1,:).*(par(1).*xx(1,8,:)+par(2).*xx(1,9,:)+par(4).*xx(1,5,:))-2.*par(2).^2.*par(3).*xx(1,1,:).*(par(1).*xx(1,9,:)+par(2).*xx(1,10,:)+par(4).*xx(1,6,:))-par(3).*xx(1,1,:).^2.*par(2).*par(4).^2];
		case 4
			J = [-par(1).*par(2).*par(3).^2.*xx(1,1,:).*xx(1,5,:)-par(2).^2.*par(3).^2.*xx(1,1,:).*xx(1,6,:)-par(2).*par(3).^2.*par(4).*xx(1,1,:).^2-par(2).*par(3).*xx(1,1,:)];
		end
	case 3
		switch np
		case 1
			J = [-2.*par(1).*par(3).^2.*par(4).*xx(1,1,:).^2-2.*par(1).*par(3).^2.*par(4).*xx(1,1,:).*xx(1,2,:)-2.*par(1).*par(3).*xx(1,1,:)];
		case 2
			J = zeros(size(xx(1,1,:)));
		case 3
			J = [-2.*par(1).^2.*par(3).*par(4).*xx(1,1,:).^2-2.*par(1).^2.*par(3).*par(4).*xx(1,1,:).*xx(1,2,:)-par(1).^2.*xx(1,1,:)];
		case 4
			J = [-par(1).^2.*par(3).^2.*xx(1,1,:).^2-par(1).^2.*par(3).^2.*xx(1,1,:).*xx(1,2,:)];
		end
	case 4
		switch np
		case 1
			J = [-2.*par(2).*par(3).^2.*par(4).*xx(1,1,:).^2-par(2).*par(3).^2.*par(4).*xx(1,1,:).*xx(1,2,:)-par(2).*par(3).^2.*par(4).*xx(1,1,:).*xx(1,3,:)-2.*par(2).*par(3).*xx(1,1,:)];
		case 2
			J = [-2.*par(1).*par(3).^2.*par(4).*xx(1,1,:).^2-par(1).*par(3).^2.*par(4).*xx(1,1,:).*xx(1,2,:)-par(1).*par(3).^2.*par(4).*xx(1,1,:).*xx(1,3,:)-2.*par(1).*par(3).*xx(1,1,:)];
		case 3
			J = [-4.*par(1).*par(2).*par(3).*par(4).*xx(1,1,:).^2-2.*par(1).*par(2).*par(3).*par(4).*xx(1,1,:).*xx(1,2,:)-2.*par(1).*par(2).*par(3).*par(4).*xx(1,1,:).*xx(1,3,:)-2.*par(1).*par(2).*xx(1,1,:)];
		case 4
			J = [-2.*par(1).*par(2).*par(3).^2.*xx(1,1,:).^2-par(1).*par(2).*par(3).^2.*xx(1,1,:).*xx(1,2,:)-par(1).*par(2).*par(3).^2.*xx(1,1,:).*xx(1,3,:)];
		end
	case 5
		switch np
		case 1
			J = zeros(size(xx(1,1,:)));
		case 2
			J = [-2.*par(2).*par(3).^2.*par(4).*xx(1,1,:).^2-2.*par(2).*par(3).^2.*par(4).*xx(1,1,:).*xx(1,3,:)-2.*par(2).*par(3).*xx(1,1,:)];
		case 3
			J = [-2.*par(2).^2.*par(3).*par(4).*xx(1,1,:).^2-2.*par(2).^2.*par(3).*par(4).*xx(1,1,:).*xx(1,3,:)-par(2).^2.*xx(1,1,:)];
		case 4
			J = [-par(2).^2.*par(3).^2.*xx(1,1,:).^2-par(2).^2.*par(3).^2.*xx(1,1,:).*xx(1,3,:)];
		end
	case 6
		switch np
		case 1
			J = [-3.*par(1).^2.*par(3).^2.*xx(1,1,:).*xx(1,2,:)-3/2.*par(3).^2.*xx(1,1,:).^2.*par(1).^2];
		case 2
			J = zeros(size(xx(1,1,:)));
		case 3
			J = [-par(1).^3.*par(3).*xx(1,1,:).^2-2.*par(1).^3.*par(3).*xx(1,1,:).*xx(1,2,:)];
		case 4
			J = zeros(size(xx(1,1,:)));
		end
	case 7
		switch np
		case 1
			J = [-3.*par(1).*par(2).*par(3).^2.*xx(1,1,:).^2-4.*par(1).*par(2).*par(3).^2.*xx(1,1,:).*xx(1,2,:)-2.*par(1).*par(2).*par(3).^2.*xx(1,1,:).*xx(1,3,:)];
		case 2
			J = [-2.*par(1).^2.*par(3).^2.*xx(1,1,:).*xx(1,2,:)-par(1).^2.*par(3).^2.*xx(1,1,:).*xx(1,3,:)-3/2.*par(3).^2.*xx(1,1,:).^2.*par(1).^2];
		case 3
			J = [-3.*par(1).^2.*par(2).*par(3).*xx(1,1,:).^2-4.*par(1).^2.*par(2).*par(3).*xx(1,1,:).*xx(1,2,:)-2.*par(1).^2.*par(2).*par(3).*xx(1,1,:).*xx(1,3,:)];
		case 4
			J = zeros(size(xx(1,1,:)));
		end
	case 8
		switch np
		case 1
			J = [-par(2).^2.*par(3).^2.*xx(1,1,:).*xx(1,2,:)-2.*par(2).^2.*par(3).^2.*xx(1,1,:).*xx(1,3,:)-3/2.*par(3).^2.*xx(1,1,:).^2.*par(2).^2];
		case 2
			J = [-3.*par(1).*par(2).*par(3).^2.*xx(1,1,:).^2-2.*par(1).*par(2).*par(3).^2.*xx(1,1,:).*xx(1,2,:)-4.*par(1).*par(2).*par(3).^2.*xx(1,1,:).*xx(1,3,:)];
		case 3
			J = [-3.*par(1).*par(2).^2.*par(3).*xx(1,1,:).^2-2.*par(1).*par(2).^2.*par(3).*xx(1,1,:).*xx(1,2,:)-4.*par(1).*par(2).^2.*par(3).*xx(1,1,:).*xx(1,3,:)];
		case 4
			J = zeros(size(xx(1,1,:)));
		end
	case 9
		switch np
		case 1
			J = zeros(size(xx(1,1,:)));
		case 2
			J = [-3.*par(2).^2.*par(3).^2.*xx(1,1,:).*xx(1,3,:)-3/2.*par(3).^2.*xx(1,1,:).^2.*par(2).^2];
		case 3
			J = [-par(2).^3.*par(3).*xx(1,1,:).^2-2.*par(2).^3.*par(3).*xx(1,1,:).*xx(1,3,:)];
		case 4
			J = zeros(size(xx(1,1,:)));
		end
	end
elseif length(nx) == 2 && isempty(np) && ~isempty(v)
	nx1 = nx(1); nx2 = nx(2);
	switch nx1
	case 0
		switch nx2
		case 0
			J = [(-par(3).^2.*par(1).*(par(4).^2.*xx(1,2,:)+2.*par(4).*(par(1).*xx(1,4,:)+par(2).*xx(1,5,:))+par(1).^2.*xx(1,7,:)+2.*par(1).*par(2).*xx(1,8,:)+par(2).^2.*xx(1,9,:))-par(3).^2.*par(2).*(par(4).^2.*xx(1,3,:)+2.*par(4).*(par(1).*xx(1,5,:)+par(2).*xx(1,6,:))+par(1).^2.*xx(1,8,:)+2.*par(1).*par(2).*xx(1,9,:)+par(2).^2.*xx(1,10,:))).*v(1,:,:)];
		case 1
			J = [(-par(1).*par(3).*par(4)-par(1).^2.*par(3).^2.*(par(1).*xx(1,7,:)+par(2).*xx(1,8,:)+par(4).*xx(1,4,:))-par(1).*par(2).*par(3).^2.*(par(1).*xx(1,8,:)+par(2).*xx(1,9,:)+par(4).*xx(1,5,:))-par(3).^2.*xx(1,1,:).*par(1).*par(4).^2).*v(1,:,:)];
		case 2
			J = [(-par(2).*par(3).*par(4)-par(1).*par(2).*par(3).^2.*(par(1).*xx(1,8,:)+par(2).*xx(1,9,:)+par(4).*xx(1,5,:))-par(2).^2.*par(3).^2.*(par(1).*xx(1,9,:)+par(2).*xx(1,10,:)+par(4).*xx(1,6,:))-par(3).^2.*xx(1,1,:).*par(2).*par(4).^2).*v(1,:,:)];
		case 3
			J = [(-2.*par(1).^2.*par(3).^2.*par(4).*xx(1,1,:)-par(1).^2.*par(3).^2.*par(4).*xx(1,2,:)-par(1).^2.*par(3)).*v(1,:,:)];
		case 4
			J = [(-4.*par(1).*par(2).*par(3).^2.*par(4).*xx(1,1,:)-par(1).*par(2).*par(3).^2.*par(4).*xx(1,2,:)-par(1).*par(2).*par(3).^2.*par(4).*xx(1,3,:)-2.*par(1).*par(2).*par(3)).*v(1,:,:)];
		case 5
			J = [(-2.*par(2).^2.*par(3).^2.*par(4).*xx(1,1,:)-par(2).^2.*par(3).^2.*par(4).*xx(1,3,:)-par(2).^2.*par(3)).*v(1,:,:)];
		case 6
			J = [(-par(1).^3.*par(3).^2.*xx(1,1,:)-par(1).^3.*par(3).^2.*xx(1,2,:)).*v(1,:,:)];
		case 7
			J = [(-3.*par(1).^2.*par(2).*par(3).^2.*xx(1,1,:)-2.*par(1).^2.*par(2).*par(3).^2.*xx(1,2,:)-par(1).^2.*par(2).*par(3).^2.*xx(1,3,:)).*v(1,:,:)];
		case 8
			J = [(-3.*par(1).*par(2).^2.*par(3).^2.*xx(1,1,:)-par(1).*par(2).^2.*par(3).^2.*xx(1,2,:)-2.*par(1).*par(2).^2.*par(3).^2.*xx(1,3,:)).*v(1,:,:)];
		case 9
			J = [(-par(2).^3.*par(3).^2.*xx(1,1,:)-par(2).^3.*par(3).^2.*xx(1,3,:)).*v(1,:,:)];
		end
	case 1
		switch nx2
		case 0
			J = [(-par(1).*par(3).*par(4)-par(1).^2.*par(3).^2.*(par(1).*xx(1,7,:)+par(2).*xx(1,8,:)+par(4).*xx(1,4,:))-par(1).*par(2).*par(3).^2.*(par(1).*xx(1,8,:)+par(2).*xx(1,9,:)+par(4).*xx(1,5,:))-par(3).^2.*xx(1,1,:).*par(1).*par(4).^2).*v(1,:,:)];
		case 1
			J = zeros(size(xx(1,1,:)));
		case 2
			J = zeros(size(xx(1,1,:)));
		case 3
			J = [-par(3).^2.*xx(1,1,:).*par(1).^2.*par(4).*v(1,:,:)];
		case 4
			J = [-par(3).^2.*xx(1,1,:).*par(1).*par(4).*par(2).*v(1,:,:)];
		case 5
			J = zeros(size(xx(1,1,:)));
		case 6
			J = [-par(3).^2.*xx(1,1,:).*par(1).^3.*v(1,:,:)];
		case 7
			J = [-2.*par(3).^2.*xx(1,1,:).*par(1).^2.*par(2).*v(1,:,:)];
		case 8
			J = [-par(3).^2.*xx(1,1,:).*par(1).*par(2).^2.*v(1,:,:)];
		case 9
			J = zeros(size(xx(1,1,:)));
		end
	case 2
		switch nx2
		case 0
			J = [(-par(2).*par(3).*par(4)-par(1).*par(2).*par(3).^2.*(par(1).*xx(1,8,:)+par(2).*xx(1,9,:)+par(4).*xx(1,5,:))-par(2).^2.*par(3).^2.*(par(1).*xx(1,9,:)+par(2).*xx(1,10,:)+par(4).*xx(1,6,:))-par(3).^2.*xx(1,1,:).*par(2).*par(4).^2).*v(1,:,:)];
		case 1
			J = zeros(size(xx(1,1,:)));
		case 2
			J = zeros(size(xx(1,1,:)));
		case 3
			J = zeros(size(xx(1,1,:)));
		case 4
			J = [-par(3).^2.*xx(1,1,:).*par(1).*par(4).*par(2).*v(1,:,:)];
		case 5
			J = [-par(3).^2.*xx(1,1,:).*par(2).^2.*par(4).*v(1,:,:)];
		case 6
			J = zeros(size(xx(1,1,:)));
		case 7
			J = [-par(3).^2.*xx(1,1,:).*par(1).^2.*par(2).*v(1,:,:)];
		case 8
			J = [-2.*par(3).^2.*xx(1,1,:).*par(1).*par(2).^2.*v(1,:,:)];
		case 9
			J = [-par(3).^2.*xx(1,1,:).*par(2).^3.*v(1,:,:)];
		end
	case 3
		switch nx2
		case 0
			J = [(-2.*par(1).^2.*par(3).^2.*par(4).*xx(1,1,:)-par(1).^2.*par(3).^2.*par(4).*xx(1,2,:)-par(1).^2.*par(3)).*v(1,:,:)];
		case 1
			J = [-par(3).^2.*xx(1,1,:).*par(1).^2.*par(4).*v(1,:,:)];
		case 2
			J = zeros(size(xx(1,1,:)));
		case 3
			J = zeros(size(xx(1,1,:)));
		case 4
			J = zeros(size(xx(1,1,:)));
		case 5
			J = zeros(size(xx(1,1,:)));
		case 6
			J = zeros(size(xx(1,1,:)));
		case 7
			J = zeros(size(xx(1,1,:)));
		case 8
			J = zeros(size(xx(1,1,:)));
		case 9
			J = zeros(size(xx(1,1,:)));
		end
	case 4
		switch nx2
		case 0
			J = [(-4.*par(1).*par(2).*par(3).^2.*par(4).*xx(1,1,:)-par(1).*par(2).*par(3).^2.*par(4).*xx(1,2,:)-par(1).*par(2).*par(3).^2.*par(4).*xx(1,3,:)-2.*par(1).*par(2).*par(3)).*v(1,:,:)];
		case 1
			J = [-par(3).^2.*xx(1,1,:).*par(1).*par(4).*par(2).*v(1,:,:)];
		case 2
			J = [-par(3).^2.*xx(1,1,:).*par(1).*par(4).*par(2).*v(1,:,:)];
		case 3
			J = zeros(size(xx(1,1,:)));
		case 4
			J = zeros(size(xx(1,1,:)));
		case 5
			J = zeros(size(xx(1,1,:)));
		case 6
			J = zeros(size(xx(1,1,:)));
		case 7
			J = zeros(size(xx(1,1,:)));
		case 8
			J = zeros(size(xx(1,1,:)));
		case 9
			J = zeros(size(xx(1,1,:)));
		end
	case 5
		switch nx2
		case 0
			J = [(-2.*par(2).^2.*par(3).^2.*par(4).*xx(1,1,:)-par(2).^2.*par(3).^2.*par(4).*xx(1,3,:)-par(2).^2.*par(3)).*v(1,:,:)];
		case 1
			J = zeros(size(xx(1,1,:)));
		case 2
			J = [-par(3).^2.*xx(1,1,:).*par(2).^2.*par(4).*v(1,:,:)];
		case 3
			J = zeros(size(xx(1,1,:)));
		case 4
			J = zeros(size(xx(1,1,:)));
		case 5
			J = zeros(size(xx(1,1,:)));
		case 6
			J = zeros(size(xx(1,1,:)));
		case 7
			J = zeros(size(xx(1,1,:)));
		case 8
			J = zeros(size(xx(1,1,:)));
		case 9
			J = zeros(size(xx(1,1,:)));
		end
	case 6
		switch nx2
		case 0
			J = [(-par(1).^3.*par(3).^2.*xx(1,1,:)-par(1).^3.*par(3).^2.*xx(1,2,:)).*v(1,:,:)];
		case 1
			J = [-par(3).^2.*xx(1,1,:).*par(1).^3.*v(1,:,:)];
		case 2
			J = zeros(size(xx(1,1,:)));
		case 3
			J = zeros(size(xx(1,1,:)));
		case 4
			J = zeros(size(xx(1,1,:)));
		case 5
			J = zeros(size(xx(1,1,:)));
		case 6
			J = zeros(size(xx(1,1,:)));
		case 7
			J = zeros(size(xx(1,1,:)));
		case 8
			J = zeros(size(xx(1,1,:)));
		case 9
			J = zeros(size(xx(1,1,:)));
		end
	case 7
		switch nx2
		case 0
			J = [(-3.*par(1).^2.*par(2).*par(3).^2.*xx(1,1,:)-2.*par(1).^2.*par(2).*par(3).^2.*xx(1,2,:)-par(1).^2.*par(2).*par(3).^2.*xx(1,3,:)).*v(1,:,:)];
		case 1
			J = [-2.*par(3).^2.*xx(1,1,:).*par(1).^2.*par(2).*v(1,:,:)];
		case 2
			J = [-par(3).^2.*xx(1,1,:).*par(1).^2.*par(2).*v(1,:,:)];
		case 3
			J = zeros(size(xx(1,1,:)));
		case 4
			J = zeros(size(xx(1,1,:)));
		case 5
			J = zeros(size(xx(1,1,:)));
		case 6
			J = zeros(size(xx(1,1,:)));
		case 7
			J = zeros(size(xx(1,1,:)));
		case 8
			J = zeros(size(xx(1,1,:)));
		case 9
			J = zeros(size(xx(1,1,:)));
		end
	case 8
		switch nx2
		case 0
			J = [(-3.*par(1).*par(2).^2.*par(3).^2.*xx(1,1,:)-par(1).*par(2).^2.*par(3).^2.*xx(1,2,:)-2.*par(1).*par(2).^2.*par(3).^2.*xx(1,3,:)).*v(1,:,:)];
		case 1
			J = [-par(3).^2.*xx(1,1,:).*par(1).*par(2).^2.*v(1,:,:)];
		case 2
			J = [-2.*par(3).^2.*xx(1,1,:).*par(1).*par(2).^2.*v(1,:,:)];
		case 3
			J = zeros(size(xx(1,1,:)));
		case 4
			J = zeros(size(xx(1,1,:)));
		case 5
			J = zeros(size(xx(1,1,:)));
		case 6
			J = zeros(size(xx(1,1,:)));
		case 7
			J = zeros(size(xx(1,1,:)));
		case 8
			J = zeros(size(xx(1,1,:)));
		case 9
			J = zeros(size(xx(1,1,:)));
		end
	case 9
		switch nx2
		case 0
			J = [(-par(2).^3.*par(3).^2.*xx(1,1,:)-par(2).^3.*par(3).^2.*xx(1,3,:)).*v(1,:,:)];
		case 1
			J = zeros(size(xx(1,1,:)));
		case 2
			J = [-par(3).^2.*xx(1,1,:).*par(2).^3.*v(1,:,:)];
		case 3
			J = zeros(size(xx(1,1,:)));
		case 4
			J = zeros(size(xx(1,1,:)));
		case 5
			J = zeros(size(xx(1,1,:)));
		case 6
			J = zeros(size(xx(1,1,:)));
		case 7
			J = zeros(size(xx(1,1,:)));
		case 8
			J = zeros(size(xx(1,1,:)));
		case 9
			J = zeros(size(xx(1,1,:)));
		end
	end
end
if isempty(J)
	display([nx np size(v)]);
	error('SYS_DERI: requested derivative could not be computed!');
end
