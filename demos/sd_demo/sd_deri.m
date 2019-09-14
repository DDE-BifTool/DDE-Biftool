%% partial derivatives of sd_demo right-hand side (vectorized in x)
%
% $Id: sd_deri.m 20 2014-04-11 19:27:33Z jan.sieber $
%
%%
function J=sd_deri(xx,par,nx,np,v)
%%
%#ok<*NASGU>
% p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11
n=size(xx,1);
nvec=size(xx,3);
J=[];
x1=xx(1,1,:);
x2=xx(2,1,:);
x3=xx(3,1,:);
x4=xx(4,1,:);
x5=xx(5,1,:); %#ok<*NASGU>

x1_tau1=xx(1,2,:);
x2_tau1=xx(2,2,:);

x1_tau2=xx(1,3,:);
x2_tau2=xx(2,3,:);

x1_tau3=xx(1,4,:);
x3_tau3=xx(3,4,:);

x1_tau4=xx(1,5,:);
x2_tau4=xx(2,5,:);

x2_tau5=xx(2,6,:);

x1_tau6=xx(1,7,:);

p1=par(1);
p2=par(2);
p3=par(3);
p4=par(4);
p5=par(5);
p6=par(6);
p7=par(7);
p8=par(8);
p9=par(9);

if length(nx)==1 && isempty(np) && isempty(v)
    %% first order derivatives wrt state variables
    J=zeros(n,n,nvec);
    if nx==0 % derivative wrt x(t)
        J(1,1,:)=-p2*x1_tau3.*x3_tau3./(p1+x2);
        J(1,2,:)=-(1-p2*x1.*x1_tau3.*x3_tau3+p3*x1_tau1.*x2_tau2)./((p1+x2).^2);
        J(2,1,:)=p4./(p1+x2);
        J(2,2,:)=-p4*x1./((p1+x2).^2);
        J(3,2,:)=p6;
        J(3,3,:)=-p6;
        J(3,4,:)=-p7*(x1_tau6-x2_tau4).*exp(-p8*x4)*(-p8);
        J(4,4,:)=x1_tau4.*exp(-p1*x4)*(-p1);
        J(5,5,:)=-3;
    elseif nx==1 % derivative wrt x(t-tau1)
        J(1,1,:)=p3*x2_tau2./(p1+x2);
        J(5,5,:)=0;
    elseif nx==2 % derivative wrt x(t-tau2)
        J(1,2,:)=p3*x1_tau1./(p1+x2);
        J(5,1,:)=3;
        J(5,5,:)=0;
    elseif nx==3 % derivative wrt x(t-tau3)
        J(1,1,:)=-p2*x1.*x3_tau3./(p1+x2);
        J(1,3,:)=-p2*x1.*x1_tau3./(p1+x2);
        J(5,5,:)=0;
    elseif nx==4 % derivative wrt x(t-tau4)
        J(3,2,:)=p7*exp(-p8*x4);
        J(4,1,:)=exp(-p1*x4);
        J(5,5,:)=0;
    elseif nx==5 % derivative wrt x(t-tau5)
        J(2,2,:)=p5*(1-tanh(x2_tau5).^2);
        J(5,5,:)=0;
    elseif nx==6 % derivative wrt x(t-tau6)
        J(3,1,:)=-p7*exp(-p8*x4);
        J(5,5,:)=0;
    end
elseif isempty(nx) && length(np)==1 && isempty(v)
    %% first order derivatives wrt parameters
    J=zeros(n,1,nvec);
    if np==1 % derivative wrt p1
        J(1,1,:)=-(1-p2*x1.*x1_tau3.*x3_tau3+p3*x1_tau1.*x2_tau2)./((p1+x2).^2);
        J(2,1,:)=-(p4*x1)./((p1+x2).^2);
        J(4,1,:)=x1_tau4.*exp(-p1*x4).*(-x4);
        J(5,1,:)=0;
    elseif np==2 % derivative wrt p2
        J(1,1,:)=-x1.*x1_tau3.*x3_tau3./(p1+x2);
        J(5,1,:)=0;
    elseif np==3 % derivative wrt p3
        J(1,1,:)=x1_tau1.*x2_tau2./(p1+x2);
        J(5,1,:)=0;
    elseif np==4 % derivative wrt p4
        J(2,1,:)=x1./(p1+x2);
        J(5,1,:)=0;
    elseif np==5 % derivative wrt p5
        J(2,1,:)=tanh(x2_tau5);
        J(5,1,:)=0;
    elseif np==6 % derivative wrt p6
        J(3,1,:)=x2-x3;
        J(5,1,:)=0;
    elseif np==7 % derivative wrt p7
        J(3,1,:)=-(x1_tau6-x2_tau4).*exp(-p8*x4);
        J(5,1,:)=0;
    elseif np==8 % derivative wrt p8
        J(3,1,:)=-p7*(x1_tau6-x2_tau4).*exp(-p8*x4).*(-x4);
        J(5,1,:)=0;
    elseif np==9 % derivative wrt p9
        J(5,1,:)=-1;
    elseif np==10 || np==11 % derivative wrt tau1, tau2
        J(5,1,:)=0;
    end
    %% second order derivatives wrt state variables
elseif length(nx)==2 && isempty(np) && ~isempty(v),
    J=zeros(n,n,nvec);
    v=reshape(v,[n,1,nvec]);
    if nx(1)==0 % first derivative wrt x(t)
        if nx(2)==0
            J(1,1,:)=p2*x1_tau3.*x3_tau3.*v(2,1,:)./((p1+x2).^2);
            J(1,2,:)=p2*x1_tau3.*x3_tau3*v(1,1,:)./((p1+x2).^2)+...
                2*(1-p2*x1.*x1_tau3.*x3_tau3+p3*x1_tau1.*x2_tau2).*v(2,1,:)/((p1+x2)^3);
            J(2,1,:)=-p4*v(2,1,:)./((p1+x2).^2);
            J(2,2,:)=-p4*v(1,1,:)./((p1+x2).^2)+2*p4*x1.*v(2,1,:)./((p1+x2).^3);
            J(3,4,:)=-p7*(x1_tau6-x2_tau4).*exp(-p8*x4)*(p8)^2.*v(4,1,:);
            J(4,4,:)=x1_tau4.*exp(-p1*x4)*(p1)^2*v(4,1,:);
            J(5,5,:)=0;
        elseif nx(2)==1
            J(1,1,:)=-p3*x2_tau2.*v(2,1,:)./((p1+x2).^2);
            J(5,5,:)=0;
        elseif nx(2)==2
            J(1,2,:)=-p3*x1_tau1.*v(2,1,:)./((p1+x2).^2);
            J(5,5,:)=0;
        elseif nx(2)==3
            J(1,1,:)=-p2*x3_tau3.*v(1,1,:)./(p1+x2)+p2*x1.*x3_tau3.*v(2,1,:)./((p1+x2).^2);
            J(1,3,:)=-p2*x1_tau3.*v(1,1,:)./(p1+x2)+p2*x1.*x1_tau3.*v(2,1,:)./((p1+x2).^2);
            J(5,5,:)=0;
        elseif nx(2)==4
            J(3,2,:)=p7*exp(-p8*x4)*(-p8).*v(4,1,:);
            J(4,1,:)=exp(-p1*x4)*(-p1).*v(4,1,:);
            J(5,5,:)=0;
        elseif nx(2)==6
            J(3,1,:)=-p7*exp(-p8*x4)*(-p8)*v(4,1,:);
            J(5,5,:)=0;
        else
            J(5,5,:)=0;
        end
    elseif nx(1)==1 % first derivative wrt x(t-tau1)
        if nx(2)==0
            J(1,2,:)=-p3*x2_tau2*v(1,1,:)./((p1+x2).^2);
            J(5,5,:)=0;
        elseif nx(2)==2
            J(1,2,:)=p3*v(1,1,:)./(p1+x2);
            J(5,5,:)=0;
        else
            J(5,5,:)=0;
        end
    elseif nx(1)==2 % first derivative wrt x(t-tau2)
        if nx(2)==0
            J(1,2,:)=-p3*x1_tau1.*v(2,1,:)./((p1+x2).^2);
            J(5,5,:)=0;
        elseif nx(2)==1
            J(1,1,:)=p3*v(2,1,:)./(p1+x2);
            J(5,5,:)=0;
        else
            J(5,5,:)=0;
        end
    elseif nx(1)==3 % first derivative wrt x(t-tau3)
        if nx(2)==0
            J(1,1,:)=(-p2*x3_tau3.*v(1,1,:)-p2*x1_tau3.*v(3,1,:))./(p1+x2);
            J(1,2,:)=p2*x1.*x3_tau3.*v(1,1,:)./((p1+x2).^2)+p2*x1.*x1_tau3.*v(3,1,:)./((p1+x2).^2);
            J(5,5,:)=0;
        elseif nx(2)==3
            J(1,3,:)=-p2*x1.*v(1,1,:)./(p1+x2);
            J(1,1,:)=-p2*x1.*v(3,1,:)./(p1+x2);
            J(5,5,:)=0;
        else
            J(5,5,:)=0;
        end
    elseif nx(1)==4 % first derivative wrt x(t-tau4)
        if nx(2)==0
            J(3,4,:)=p7*exp(-p8*x4)*(-p8).*v(2,1,:);
            J(4,4,:)=exp(-p1*x4)*(-p1).*v(1,1,:);
            J(5,5,:)=0;
        else
            J(5,5,:)=0;
        end
    elseif nx(1)==5 % first derivative wrt x(t-tau5)
        if nx(2)==5
            th=tanh(x2_tau5);
            J(2,2,:)=-2*p5*th.*(1-th.*th).*v(2,1,:);
            J(5,5,:)=0;
        else
            J(5,5,:)=0;
        end
    elseif nx(1)==6 % first derivative wrt x(t-tau6)
        if nx(2)==0
            J(3,4,:)=-p7*exp(-p8*x4)*(-p8).*v(1,1,:);
            J(5,5,:)=0;
        else
            J(5,5,:)=0;
        end
    end
elseif length(nx)==1 && length(np)==1 && isempty(v)
    %% mixed state parameter derivatives:
    J=zeros(n,n,nvec);
    if nx==0 % derivative wrt x(t)
        if np==1 % derivative wrt p1
            J(1,1,:)=p2*x1_tau3.*x3_tau3./((p1+x2).^2);
            J(1,2,:)=2*(1-p2*x1.*x1_tau3.*x3_tau3+p3*x1_tau1.*x2_tau2)./((p1+x2).^3);
            J(2,1,:)=-p4./((p1+x2).^2);
            J(2,2,:)=2*p4*x1./((p1+x2).^3);
            J(4,4,:)=x1_tau4.*exp(-p1*x4).*(-x4-1);
            J(5,5,:)=0;
        elseif np==2 % derivative wrt p2
            J(1,1,:)=-x1_tau3.*x3_tau3./(p1+x2);
            J(1,2,:)=x1.*x1_tau3.*x3_tau3./((p1+x2).^2);
            J(5,5,:)=0;
        elseif np==3 % derivative wrt p3
            J(1,2,:)=-x1_tau1.*x2_tau2./((p1+x2).^2);
            J(5,5,:)=0;
        elseif np==4 % derivative wrt p4
            J(2,1,:)=1./(p1+x2);
            J(2,2,:)=-x1./((p1+x2).^2);
            J(5,5,:)=0;
        elseif np==6 % derivative wrt p6
            J(3,2,:)=1;
            J(3,3,:)=-1;
            J(5,5,:)=0;
        elseif np==7 % derivative wrt p7
            J(3,4,:)=-(x1_tau6-x2_tau4).*exp(-p8*x4)*(-p8);
            J(5,5,:)=0;
        elseif np==8 % derivative wrt p8
            J(3,4,:)=-p7*(x1_tau6-x2_tau4).*exp(-p8*x4).*(-x4-1);
            J(5,5,:)=0;
        end
    elseif nx==1 % derivative wrt x(t-tau1)
        if np==1 % derivative wrt p1
            J(1,1,:)=-p3*x2_tau2./((p1+x2).^2);
            J(5,5,:)=0;
        elseif np==3 % derivative wrt p3
            J(1,1,:)=x2_tau2./(p1+x2);
            J(5,5,:)=0;
        end
    elseif nx==2 % derivative wrt x(t-tau2)
        if np==1 % derivative wrt p1
            J(1,2,:)=-p3*x1_tau1./((p1+x2).^2);
            J(5,5,:)=0;
        elseif np==3 % derivative wrt p3
            J(1,2,:)=x1_tau1./(p1+x2);
            J(5,5,:)=0;
        end
    elseif nx==3 % derivative wrt x(t-tau3)
        if np==1 % derivative wrt p1
            J(1,1,:)=p2*x1.*x3_tau3./((p1+x2).^2);
            J(1,3,:)=p2*x1.*x1_tau3./((p1+x2).^2);
            J(5,5,:)=0;
        elseif np==2 % derivative wrt p2
            J(1,1,:)=-x1.*x3_tau3./(p1+x2);
            J(1,3,:)=-x1.*x1_tau3./(p1+x2);
            J(5,5,:)=0;
        end
    elseif nx==4 % derivative wrt x(t-tau4)
        if np==1 % derivative wrt p1
            J(4,1,:)=exp(-p1*x4).*(-x4);
            J(5,5,:)=0;
        elseif np==7 % derivative wrt p7
            J(3,2,:)=exp(-p8*x4);
            J(5,5,:)=0;
        elseif np==8 % derivative wrt p8
            J(3,2,:)=p7*exp(-p8*x4).*(-x4);
            J(5,5,:)=0;
        end
    elseif nx==5 % derivative wrt x(t-tau5)
        if np==5 % derivative wrt p5
            J(2,2,:)=1-tanh(x2_tau5).^2;
            J(5,5,:)=0;
        end
    elseif nx==6 % derivative wrt x(t-tau6)
        if np==7 % derivative wrt p7
            J(3,1,:)=-exp(-p8*x4);
            J(5,5,:)=0;
        elseif np==8 % derivative wrt p8
            J(3,1,:)=-p7*exp(-p8*x4).*(-x4);
            J(5,5,:)=0;
        end
    end
end

if isempty(J)
    error('SYS_DERI:order',...
        'SYS_DERI: requested derivative nx=%, np=%d, size(v)=(%d,%d) could not be computed!',...
        nx,np,size(v));
end
end
