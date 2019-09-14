function J=sys_deri_psolfold(xx,par,nx,np,v,hjac,...
    ind_beta,ind_period,dim,ntau,funcs,ext_rhs,use_df_deriv)
%% partial derivatives of r.h.s of extended DDE for fold of periodic orbits w constant delays
%
% $Id: sys_deri_psolfold.m 346 2019-05-13 05:41:50Z jansieber $
%%
ind_tau=funcs.sys_tau();
nvec=size(xx,3);
pfuncs=struct('sys_rhs',ext_rhs);
if length(nx)==1 && isempty(np) && isempty(v)
    %% first order derivatives of the state are prepared:
    if ~use_df_deriv
        J=sys_drhs_dx_psolfold(nx,xx,par(1:ind_beta-1),par(ind_beta),...
            par(ind_period),[0,par(ind_tau)],dim,ntau,funcs.sys_deri);
    else
        J=df_deriv(pfuncs,xx,par,nx,np,v,hjac);
    end
elseif isempty(nx) && length(np)==1 && isempty(v)
    %% first-order parameter derivatives
    if np<ind_beta
        %% derivatives wrt problem parameters
        if ~use_df_deriv
            J=sys_drhs_dp_POfold(np,xx,par(1:ind_beta-1),par(ind_beta),...
                par(ind_period),ind_tau,funcs.sys_rhs,dim,xtau_ind,funcs.sys_deri);
        else
            J=df_deriv(pfuncs,xx,par,nx,np,v,hjac);
        end
    elseif  np==ind_beta
        %% derivative wrt beta
        par(ind_beta)=1;
        fp=pfuncs.sys_rhs([xx(1:dim,:,:);zeros(dim,size(xx,2),nvec)],par);
        J=[zeros(dim,1,nvec);reshape(fp(dim+1:end,:),[dim,1,nvec])];
    elseif np==ind_period
        %% derivative wrt copy of period
        par(ind_period)=-par(ind_period)^2;
        fp=pfuncs.sys_rhs([xx(1:dim,:,:);zeros(dim,size(xx,2),nvec)],par);
        J=[zeros(dim,1,nvec);reshape(fp(dim+1:end,:),[dim,1,nvec])];
    elseif np>ind_period
        J=zeros(2*dim,1,nvec);
    else
        J=df_deriv(pfuncs,xx,par,nx,np,v,hjac);
    end
else
    %% shouldn't be needed
    J=df_deriv(pfuncs,xx,par,nx,np,v,hjac);
end
end

function J=sys_drhs_dx_psolfold(ind,x,p,beta,period,tau,dim,ntau,sys_deri)
%% derivative of rhs of extended DDE for fold of periodic orbits wrt to x(:,i)
nvec=size(x,3);
xall=x(1:dim,:,:);
x0=xall(:,1:ntau+1,:);
v=x(dim+1:end,1:ntau+1,:);
xp=xall(1:dim,ntau+1+(1:ntau),:);
if ind==0
    Jxx=sys_deri(x0,p,ind,[],[]);
    Jxv=zeros(dim,dim,nvec);
    Jvv=Jxx;
    Jvx=Jxx*beta/period+sys_deri(x0,p,[0,ind],[],v(:,1,:));
    for i=2:size(x0,2)
        fdev=xp(:,i-1,:)*tau(i)*beta/period;
        deviation=v(:,i,:)+reshape(fdev,dim,1,[]);
        Jvx=Jvx+sys_deri(x0,p,[i-1,ind],[],deviation);
    end
    if ind>0
        x_ind=xall(:,xtau_ind(ind+1,:),:);
        Jvx=Jvx+VAopX(Jxx,sys_deri(x_ind,p,0,[],[])*tau(ind+1)*beta/period,'*');
    end
else
    Jxx=zeros(dim,dim,nvec);
    Jxv=zeros(dim,dim,nvec);
    Jvv=zeros(dim,dim,nvec);
    x_ir=xall(:,xtau_ind(ir,:),:);
    x_ic=xall(:,xtau_ind(ic,:),:);
    Jtau=VAopX(tau(ir)*sys_deri(x0,p,ir-1,[],[]),sys_deri(x_ir,p,ic-1,[],[]),'*');
    if ir~=ic
        Jtau=Jtau+VAopX(tau(ic)*sys_deri(x0,p,ic-1,[],[]),sys_deri(x_ic,p,ir-1,[],[]),'*');
    end
    Jvx=beta/period*Jtau;
end
J=cat(1,...
    cat(2,Jxx,Jxv),...
    cat(2,Jvx,Jvv));
end
function J=sys_drhs_dp_POfold(ind,x,p,beta,period,ind_tau,sys_rhs,dim,xtau_ind,sys_deri)
%% derivative of rhs of extended DDE for fold of periodic orbits wrt to x(:,i)
xall=x(1:dim,:,:);
x0=xall(:,xtau_ind(1,:),:);
v=x(dim+1:end,xtau_ind(1,:),:);
tau=[0,p(ind_tau)];
Jxp=sys_deri(x0,p,[],ind,[]);
Jvp=Jxp*beta/period+VAopX(sys_deri(x0,p,0,ind,[]),v(:,1,:),'*');
for i=2:size(x0,2)
    x_ind=xall(:,xtau_ind(i,:),:);
    fdev=sys_rhs(x_ind,p)*tau(i)*beta/period;
    deviation=v(:,i,:)+reshape(fdev,dim,1,[]);
    difp=sys_deri(x0,p,i-1,[],[]);
    Jvp=Jvp+VAopX(sys_deri(x0,p,i-1,ind,[]),deviation,'*')+...
        VAopX(difp,sys_deri(x_ind,p,[],ind,[])*beta*tau(i)/period,'*');
    if ind==ind_tau(i-1)
        dfdev=VAopX(difp*beta/period,sys_rhs(x_ind,p),'*');
        Jvp=Jvp+reshape(dfdev,dim,1,[]);
    end
end
J=cat(1,Jxp,Jvp);
end
