function J=sys_deri_TorusBif(xx,par,nx,np,v,hjac,...
    ip,funcs,trfuncs)
%% partial derivatives of r.h.s of extended DDE for torus or period doubling bifurcation
%
% $Id: sys_deri_TorusBif.m 369 2019-08-27 00:07:02Z jansieber $
%
user_lhs_matrix=funcs.lhs_matrix(ip.dim);
ind_tau=funcs.sys_tau();
if length(nx)==1 && isempty(np) && isempty(v)
    %% first order derivatives of the state
    J=sys_drhs_dx_TorusBif(nx,xx,par(1:ip.nuserpar),par(ip.omega),...
        par(ip.period),[0,par(ind_tau)],ip.dim,funcs.sys_deri,user_lhs_matrix);
elseif isempty(nx) && length(np)==1 && isempty(v)
    %% first-order parameter derivatives
    J=sys_drhs_dp_TorusBif(np,xx,par,ip.omega,ip.period,ind_tau,...
        ip.dim,funcs.sys_deri,user_lhs_matrix);
else
    %% shouldn't be needed
    J=df_deriv(trfuncs,xx,par,nx,np,v,hjac);
end
end

function J=sys_drhs_dx_TorusBif(ind,x,p,omega,period,tau,dim,sys_deri,lhs_matrix)
%% derivative of rhs of extended DDE for torus bif wrt to x(:,i)
nvec=size(x,3);
x0=x(1:dim,:,:);
u=x(dim+1:2*dim,:,:);
v=x(2*dim+1:end,:,:);
Jxx=sys_deri(x0,p,ind,[],[]);
Jxu=zeros(dim,dim,nvec);
Jxv=zeros(dim,dim,nvec);
Jux=zeros(dim,dim,nvec);
Juu=zeros(dim,dim,nvec);
Juv=zeros(dim,dim,nvec);
Jvx=zeros(dim,dim,nvec);
Jvu=zeros(dim,dim,nvec);
Jvv=zeros(dim,dim,nvec);
PiOmT=pi*omega/period;
pid=PiOmT*lhs_matrix*eye(dim);
pid=pid(:,:,ones(1,nvec));
if ind==0
    Juv=Juv+pid;
    Jvu=Jvu-pid;
end
for i=1:size(x0,2)
    c=cos(pi*omega/period*tau(i));
    s=sin(pi*omega/period*tau(i));
    Jux=Jux+sys_deri(x0,p,[i-1,ind],[], c*u(:,i,:)+s*v(:,i,:));
    Jvx=Jvx+sys_deri(x0,p,[i-1,ind],[],-s*u(:,i,:)+c*v(:,i,:));
    if i-1==ind
        difp=sys_deri(x0,p,ind,[],[]);
        Juu=Juu+difp*c;
        Juv=Juv+difp*s;
        Jvu=Jvu-difp*s;
        Jvv=Jvv+difp*c;
    end
end
J=cat(1,...
    cat(2,Jxx,Jxu,Jxv),...
    cat(2,Jux,Juu,Juv),...
    cat(2,Jvx,Jvu,Jvv));
end
function J=sys_drhs_dp_TorusBif(ind,x,par,ind_omega,ind_period,ind_tau,...
        dim,sys_deri,lhs_matrix)
%% derivative of rhs of extended DDE for torus bif wrt to parameters
nvec=size(x,3);
x0=x(1:dim,:,:);
u=x(dim+1:2*dim,:,:);
v=x(2*dim+1:end,:,:);
omega=par(ind_omega);
period=par(ind_period);
tau=[0,par(ind_tau)];
p=par(1:ind_omega-1);
PiOmT=pi*omega/period;
if ind<ind_omega
    %% derivative wrt system parameter
    Jx=sys_deri(x0,p,[],ind,[]);
    Ju=zeros(dim,1,nvec);
    Jv=zeros(dim,1,nvec);
    for i=1:size(x0,2)
        c=cos(PiOmT*tau(i));
        s=sin(PiOmT*tau(i));
        difxp=sys_deri(x0,p,i-1,ind,[]);
        Ju=Ju+VAopX(difxp, c*u(:,i,:)+s*v(:,i,:),'*');
        Jv=Jv+VAopX(difxp,-s*u(:,i,:)+c*v(:,i,:),'*');
        if i>1 && ind==ind_tau(i-1)
            difx=sys_deri(x0,p,i-1,[],[]);
            Ju=Ju+VAopX(difx*PiOmT,-s*u(:,i,:)+c*v(:,i,:),'*');
            Jv=Jv+VAopX(difx*PiOmT,-c*u(:,i,:)-s*v(:,i,:),'*');
        end
    end
    J=cat(1,Jx,Ju,Jv);
elseif ind==ind_omega || ind==ind_period
    Jx=zeros(dim,1,nvec);
    sq=@(x)reshape(x,size(x,1),[]);
    ex=@(x)reshape(x,size(x,1),1,[]);
    Ju=pi*ex(lhs_matrix*sq(v(:,1,:)));
    Jv=-pi*ex(lhs_matrix*sq(u(:,1,:)));
    for i=1:size(x0,2)
        c=cos(PiOmT*tau(i));
        s=sin(PiOmT*tau(i));
        difx=sys_deri(x0,p,i-1,[],[]);
        Ju=Ju+VAopX(difx*pi*tau(i),-s*u(:,i,:)+c*v(:,i,:),'*');
        Jv=Jv+VAopX(difx*pi*tau(i),-c*u(:,i,:)-s*v(:,i,:),'*');
    end
    J=cat(1,Jx,Ju,Jv);
    if ind==ind_omega
        J=J/period;
    else
        J=-J*omega/period^2;
    end
end
end
