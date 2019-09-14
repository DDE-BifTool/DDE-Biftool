function dtau=sd_dtau(ind,xx,par,nx,np)
% p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11
dtau=[];
if length(nx)==1 && isempty(np)
    %% first order derivatives wrt state variables:
    if nx==0 % derivative wrt x(t)
        if ind==3
            dtau(1:5)=0;
            dtau(2)=par(5)*par(10)*xx(2,2);
        elseif ind==4
            dtau(1)=xx(2,3)/(1+xx(1,1)*xx(2,3))^2;
            dtau(2:5)=0;
        elseif ind==5
            dtau(1:5)=0;
            dtau(4)=1;
        elseif ind==6
            dtau(5)=1;
        else
            dtau(1:5)=0;
        end;
    elseif nx==1 % derivative wrt x(t-tau1)
        if ind==3
            dtau(1:5)=0;
            dtau(2)=par(5)*par(10)*xx(2,1);
        else
            dtau(1:5)=0;
        end;
    elseif nx==2 % derivative wrt x(t-tau2)
        if ind==4
            dtau(1:5)=0;
            dtau(2)=xx(1,1)/(1+xx(1,1)*xx(2,3))^2;
        else
            dtau(1:5)=0;
        end;
    else
        dtau(1:5)=0;
    end;
elseif isempty(nx) && length(np)==1,
    %% First order derivatives wrt parameters
    if ind==1 && np==10
        dtau=1;
    elseif ind==2 && np==11
        dtau=1;
    elseif ind==3 && np==5
        dtau=par(10)*xx(2,1)*xx(2,2);
    elseif ind==3 && np==10
        dtau=par(5)*xx(2,1)*xx(2,2);
    else
        dtau=0;
    end;
elseif length(nx)==2 && isempty(np),
    %% Second order derivatives wrt state variables
    dtau=zeros(5);
    if ind==3
        if (nx(1)==0 && nx(2)==1) || (nx(1)==1 && nx(2)==0)
            dtau(2,2)=par(5)*par(10);
        end
    elseif ind==4
        if nx(1)==0 && nx(2)==0
            dtau(1,1)=-2*xx(2,3)*xx(2,3)/(1+xx(1,1)*xx(2,2))^3;
        elseif nx(1)==0 && nx(2)==2
            dtau(1,2)=(1-xx(1,1)*xx(2,3))/(1+xx(1,1)*xx(2,2))^3;
        elseif nx(1)==2 && nx(2)==0
            dtau(2,1)=(1-xx(1,1)*xx(2,3))/(1+xx(1,1)*xx(2,2))^3;
        elseif nx(1)==2 && nx(2)==2
            dtau(2,2)=-2*xx(1,1)*xx(1,1)/(1+xx(1,1)*xx(2,2))^3;
        end
    end
elseif length(nx)==1 && length(np)==1,
    %% Mixed state parameter derivatives
    dtau(1:5)=0;
    if ind==3
        if nx==0 && np==5
            dtau(2)=par(10)*xx(2,2);
        elseif nx==0 && np==10
            dtau(2)=par(5)*xx(2,2);
        elseif nx==1 && np==5
            dtau(2)=par(10)*xx(2,1);
        elseif nx==1 && np==10
            dtau(2)=par(5)*xx(2,1);
        end
    end
end
%% Otherwise raise error
% Raise error if the requested derivative does not exist
if isempty(dtau)
    error(['SYS_DTAU, delay %d: requested derivative ',...
        'nx=%d, np=%d does not exist!'],ind, nx, np);
end
end

