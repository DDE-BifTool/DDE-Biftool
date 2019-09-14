function dtau=df_derit(funcs,delay_nr,xx,par,nx,np,hjac)
%% approximate derivatives of delays wrt state and parameter with finite difference
%
% function dtau=df_derit(funcs,delay_nr,xx,par,nx,np,hjac)
% INPUT:
%   funcs  problem functions (only sys_tau field needed)
%   delay_nr delay number
%	xx state variable and delayed state variables columnwise
%	par list of parameter values
%	nx empty or list of requested state-derivatives (numbers of delay or zero)
%	np empty or list of requested parameter-derivatives
%   hjac (optional) relative finite differencing step
% OUTPUT:
%	dtau result of derivatives on delay function
% COMMENT:
%	the numerical derivatives are evaluated using forward differences

% (c) DDE-BIFTOOL v. 2.00, 30/11/2001
%
% $Id: df_derit.m 296 2018-09-24 21:36:56Z jansieber $
%
%%
% first order derivative discretisation parameters:
sys_tau=funcs.sys_tau;
if nargin==6
    hjac=1e-6;
end

abs_eps_x1=hjac;
abs_eps_x2=hjac;
abs_eps_p1=hjac;
abs_eps_p2=hjac;
rel_eps_x1=hjac;
rel_eps_x2=hjac;
rel_eps_p1=hjac;
rel_eps_p2=hjac;

xsize=size(xx);
n=xsize(1);
ndelays=xsize(2);
vecdim=xsize(3:end);
nvec=prod(vecdim);
xx=reshape(xx,[n,ndelays,nvec]);
dtau=[];

% first order state derivatives:
if length(nx)==1 && isempty(np)
    tau=sys_tau(delay_nr,xx,par);
    dtau=NaN(1,n,nvec);
    for j=1:n
        xx_eps=xx;
        eps=abs_eps_x1+rel_eps_x1*abs(xx(j,nx+1,:));
        xx_eps(j,nx+1,:)=xx(j,nx+1,:)+eps;
        dtau(1,j,:)=(sys_tau(delay_nr,xx_eps,par)-tau)./eps;
    end
    % first order parameter derivatives:
elseif isempty(nx) && length(np)==1,
    tau=sys_tau(delay_nr,xx,par);
    par_eps=par;
    eps=abs_eps_p1+rel_eps_p1*abs(par(np));
    par_eps(np)=par(np)+eps;
    dtau=(sys_tau(delay_nr,xx,par_eps)-tau)/eps;
    % second order state derivatives:
elseif length(nx)==2 && isempty(np)
    dtau_init=(df_derit(funcs,delay_nr,xx,par,nx(1),[],hjac));
    dtau=NaN(n,n,nvec);
    on=ones(n,1);
    for j=1:n
        xx_eps=xx;
        eps=abs_eps_x2+rel_eps_x2*abs(xx(j,nx(2)+1,:));
        xx_eps(j,nx(2)+1)=xx_eps(j,nx(2)+1,:)+eps;
        dtau(j,:,:)=((df_derit(funcs,delay_nr,xx_eps,par,nx(1),[],hjac))-dtau_init)./eps(1,on,:);
    end
    dtau=permute(dtau,[2,1,3]);
    % mixed state parameter derivatives:
elseif length(nx)==1 && ~isempty(np)
    dtau=df_derit(funcs,delay_nr,xx,par,nx(1),[],hjac);
    par_eps=par;
    eps=abs_eps_p2+rel_eps_p2*abs(par(np));
    par_eps(np)=par(np)+eps;
    dtau=(df_derit(funcs,delay_nr,xx,par_eps,nx(1),[],hjac)-dtau)/eps;
end;

if isempty(dtau)
    error('SYS_DTAU: requested derivative of %d for nx=%d np=%d does not exist!',...
        delay_nr,nx,np);
end
dtau=reshape(dtau,[size(dtau,1),size(dtau,2),vecdim]);
end
