function y=sys_rhs_SD_POfold(x,p,ip,funcs,sys_deri,sys_dtau)
%% rhs for fold of periodic orbits in SD-DDEs
% xtau_ind(i,:) contains which columns of x(1:dim,:) correspond to
% x(tau(i)+tau(1:end))
% argument sys_deri can be either numeric (say, 1e-3, then a
% finite-difference approximation is used for the linearized system), or a
% function providing the partial derivatives of the DDE's sys_rhs
%
% $Id: sys_rhs_SD_POfold.m 369 2019-08-27 00:07:02Z jansieber $
%
ntau=ip.orig_ntau;
dim=ip.dim;
beta=p(ip.beta);
period=p(ip.period);
nvec=size(x,3);
xall=x(1:dim,:,:);
x0=xall(:,1:ntau+1,:);
xp=xall(:,ntau+1+(1:ntau),:);
v=x(dim+1:end,1:ntau+1,:);
user_p=p(1:ip.nuserpar);
y0=funcs.sys_rhs(x0,user_p);
if isnumeric(sys_deri) 
    %% no user-provided derivative (sys_deri is size of deviation)
    df=@(dev,ind)app_dir_deriv(@(x)funcs.sys_rhs(x,user_p),x0,dev,ind,sys_deri);
else
    %% sys_deri is user-provided function
    df=@(dev,ind)VAopX(sys_deri(x0,user_p,ind-1,[],[]),dev,'*');
end
if isnumeric(sys_dtau) 
    %% no user-provided derivative (sys_dtau is size of deviation)
    dtau=@(dev,itau,ind)app_dir_deriv(@(x)funcs.sys_tau(itau,x,user_p),...
        x0(:,1:itau,:),dev,ind,sys_dtau);
else
    %% sys_dtau is user-provided function
    dtau=@(dev,itau,ind)VAopX(sys_dtau(itau,x0(:,1:itau,:),user_p,ind-1,[]),dev,'*');
end
%% accumulate xpj=x'(-tau_j), xxj=dxj/dx*v and xTj=dxj/dperiod
on=ones(dim,1);
xx=NaN(dim,ntau+1,nvec);
xT=xx;
xx(:,1,:)=v(:,1,:);
tau=NaN(1,ntau+1,nvec);
xT(:,1,:)=0;
tau(1,1,:)=0; % count of tau's includes tau0=0
for j=2:ntau+1
    tau(1,j,:)=funcs.sys_tau(j-1,x0(:,1:j-1,:),user_p);
    sumdtau_xx=zeros(1,1,nvec);
    sumdtau_xT=zeros(1,1,nvec);
    for k=1:j-1
        sumdtau_xx=sumdtau_xx+reshape(dtau(xx(:,k,:),j-1,k),size(sumdtau_xx));
        sumdtau_xT=sumdtau_xT+reshape(dtau(xT(:,k,:),j-1,k),size(sumdtau_xT));
    end
    xx(:,j,:)=v(:,j,:)-xp(:,j-1,:).*sumdtau_xx(on,1,:);
    xT(:,j,:)=xp(:,j-1,:).*(tau(on,j,:)/period-sumdtau_xT(on,1,:));
end

%% accumulate rhs
y1=beta/period*y0;
for j=1:ntau+1
    dev=xx(:,j,:)+beta*xT(:,j,:);
    y1=y1+reshape(df(dev,j),dim,nvec);
end
y=cat(1,y0,y1);
end
