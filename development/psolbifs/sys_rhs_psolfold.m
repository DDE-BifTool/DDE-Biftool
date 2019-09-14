function y=sys_rhs_psolfold(x,p,phat,dx,dp,dphat,...
    funcs,free_par_orig,free_par,dtaudp)
%% rhs of extended DDE for fold of periodic orbits (const delays)
% first x(:,2:ntau+1,:) contain x(t-tau(j)), ntau+(1:ntau)
% contain x'(t-tau(j)), j=1:ntau
% argument sys_deri can be either numeric (say, 1e-3, then a
% finite-difference approximation is used for the linearized system), or a
% function providing the partial derivatives of the DDE's sys_rhs
%
% $Id: sys_rhs_psolfold.m 362 2019-07-14 15:49:40Z jansieber $
%
dim=size(x,1)/2;
nvec=size(x,3);
nfp=length(free_par);
nfpo=length(free_par_orig);
ov=ones(1,nvec);
on=ones(dim,1);
if ~funcs.tp_del
    itau=funcs.sys_tau();
    ntau=length(itau)+1;
    tau=[0,p(itau)]';
    dtaudphat=dtaudp(:,free_par_orig); % deriv of delay wrt orig. var. pars
else
    ntau=funcs.sys_ntau()+1;
end
x0=x(1:dim,1:ntau,:);            % x(t-tau)
v=x(dim+1:end,1:ntau,:);        % v(t-tau)
xp=x(1:dim,ntau+(1:ntau),:); % x'(t-tau)
beta=phat(1);                   % var. par for period   
ph=reshape(phat(2:end),[],1);   % var. parameters of orig. problem
ry=@(y)reshape(y,dim,nvec);
rx=@(x)reshape(x,dim,[],nvec);
%% compute deviations entering the variational eq
if ~funcs.tp_del % constant delays
    fac=beta*tau-dtaudphat*ph;      % factor for x' in var. eq. (d x 1)
    facv=repmat(fac(:)',dim,1,nvec);
    dev=reshape(v+xp.*facv,[],nvec); %
else % state-dependent delays
    dev=zeros(dim,ntau,nvec);
    dev(:,1,:)=v(:,1,:);
    for k=2:ntau
        tau=funcs.sys_tau(k-1,x0,p);
        [Jtx,Jtp]=deri1(funcs.sys_dtau,x0,p,free_par_orig,k-1);
        Jtph=Jtp(:,free_par_orig,:);
        facp=beta*tau-VAopX(Jtph,phat(:,ov),'*');
        dev(:,k,:)=v(:,k,:)+xp(:,k,:).*rx(facp(on,:));
        facx=VAopX(Jtx(:,1:dim*(k-1),:),reshape(dev(:,1:k-1,:),dim*(k-1),nvec),'*');
        dev(:,k,:)=dev(:,k,:)-xp(:,k,:).*rx(facx(on,:));
    end
    dev=reshape(dev,[],nvec);
end
%% partial derivatives wrt p and delayed x's
[Jx,Jp]=deri1(funcs.sys_deri,x0,p,free_par); % Jacobians
Jph=Jp(:,free_par_orig,:);
if isempty(dx) 
    %% sys_rhs - partial derivative wrt period
    %% partial derivative wrt x,parameters
    y0=funcs.sys_rhs(x0,p);
    y1=ry(beta*xp(:,1,:));
    y1=y1+VAopX(Jx,dev,'*')+VAopX(Jph,ph(:,ov),'*');
else
    %% sys_deri (valid only for constant delays)
    dx0=dx(1:dim,1:ntau,:);            % dx(t-tau)
    dv=dx(dim+1:end,1:ntau,:);        % dv(t-tau)
    dxp=dx(1:dim,ntau+(1:ntau),:); % x'(t-tau)
    dp=reshape(dp,[],1);
    dbeta=dphat(1);
    dph=reshape(dphat(2:end),[],1);   % dev of var. parameters of orig. problem
    rJ=@(J)reshape(J,dim,[],nvec);
    y0=VAopX(Jx,reshape(dx0,[],nvec),'*')+VAopX(Jp(:,free_par,:),dp(free_par,ov),'*')+...
        VAopX(Jph,dph(:,ov),'*');
    [Jxx,Jxp,Jpp]=deri2(funcs.sys_deri,x0,p,free_par,free_par_orig,dx0);
    Jpp=Jpp(:,free_par,free_par_orig,:);
    Jxph=Jxp(:,:,free_par_orig,:);
    Jxp=Jxp(:,:,free_par,:);
    y1=ry(dbeta*xp(:,1,:)+beta*dxp(:,1,:));
    dfac=beta*dtaudp*dp+tau*dbeta-dtaudphat*dph;      % factor for x' in var. eq. (d x 1)
    dfacv=repmat(dfac(:)',dim,1,nvec);
    ydev=reshape(dv+dxp.*facv+xp.*dfacv,[],nvec);
    y1=y1+VAopX(Jx,ydev,'*')+VAopX(Jph,dph(:,ov),'*');
    for k=1:ntau
        y1=y1+VAopX(rJ(Jxx(:,:,k,:)),dev,'*');
    end
    for i=1:nfp
        y1=y1+VAopX(rJ(Jxp(:,:,i,:)),dev,'*')*dp(free_par(i));
    end
    for j=1:nfpo
        y1=y1+VAopX(rJ(Jxph(:,:,j,:)),dev,'*')*phat(j);
        y1=y1+VAopX(rJ(Jpp(:,:,j,:)),dp(free_par),'*')*phat(j);
    end
end
y=cat(1,y0,y1);
end
%%
function [Jx,Jp]=deri1(df_inp,x,p,free_par,varargin)
istau=~isempty(varargin);
df=@(ix,ip,dev)df_inp(varargin{:},x,p,ix,ip,dev);
[dim,ntau,nvec]=size(x);
ntau_c=ntau;
jdim1=dim;
if istau
    ntau_c=varargin{1}-1;
    jdim1=1;
end
Jx=zeros(jdim1,dim,ntau,nvec);
for i=1:ntau_c
    Jx(:,:,i,:)=reshape(df(i-1,[],[]),jdim1,dim,1,nvec);
end
Jx=reshape(Jx,jdim1,dim*ntau,nvec);
Jp=zeros(jdim1,size(p,2),nvec);
for j=1:length(free_par)
    Jp(:,free_par(j),:)=reshape(df([],free_par(j),[]),jdim1,1,nvec);
end
end
%%
function [Jxx,Jxp,Jpp]=deri2(df_inp,x,p,free_par,free_par_orig,xdev)
df=@(ix,ip,dev)df_inp(x,p,ix,ip,dev);
[dim,ntau,nvec]=size(x);
nfp=length(free_par);
nfpo=length(free_par_orig);
np=size(p,2);
Jxx=NaN(dim,dim,ntau,ntau,nvec);
Jxp=zeros(dim,dim,ntau,np,nvec);
Jpp=zeros(dim,np,np,nvec);
for i=1:ntau
    for k=1:ntau
        Jxx(:,:,i,k,:)=reshape(df([i,k]-1,[],xdev(:,k,:)),dim,dim,1,1,nvec);
    end
    for j=1:nfp
        Jxp(:,:,i,free_par(j),:)=reshape(df(i-1,free_par(j),[]),dim,dim,1,1,nvec);
    end
end
Jxx=reshape(Jxx,dim,dim*ntau,ntau,nvec);
Jxp=reshape(Jxp,dim,dim*ntau,np,nvec);
for l=1:nfp
    for j=1:nfpo
        Jpp(:,free_par(l),free_par(j),:)=reshape(df(free_par(l),free_par(j),[]),dim,1,1,nvec);
    end
end
end
