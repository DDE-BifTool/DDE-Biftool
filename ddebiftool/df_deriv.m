function J=df_deriv(funcs,xx,par,nx,np,v,hjac)
%% approximate derivatives of r.h.s wrt state and parameter with finite differences
%
% function J=df_deriv(funcs,xx,par,nx,np,v [,hjac])
% INPUT:
%   funcs  problem functions (only sys_rhs field needed)
%	xx state variable and delayed state variables columnwise
%	par list of parameter values
%	nx empty or list of requested state-derivatives (numbers of delay or zero)
%	np empty or list of requested parameter-derivatives
%	v matrix to multiply result with for 2nd xx derivative
%   hjac (optional) length of deviation (absolute and relative)
% OUTPUT:
%	J result of derivatives on righthandside multiplied with v
% COMMENT:
%	the numerical derivatives are evaluated using forward differences

% (c) DDE-BIFTOOL v. 1.00, 11/03/2000
%
% $Id: df_deriv.m 296 2018-09-24 21:36:56Z jansieber $
%
%%
% first order derivative discretisation parameters:

sys_rhs=funcs.sys_rhs;
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
J=[];

% first order derivatives of the state:
if length(nx)==1 && isempty(np) && isempty(v),
    J=zeros(n,n,nvec);
    on=ones(n,1);
    f0=sys_rhs(xx,par);
    for j=1:n
        xx_eps=xx;
        eps=abs_eps_x1+rel_eps_x1*abs(xx(j,nx+1,:));
        xx_eps(j,nx+1,:)=xx(j,nx+1,:)+eps;
        J(:,j,:)=reshape(sys_rhs(xx_eps,par)-f0,[n,1,nvec])./eps(on,1,:);
    end;
    % first order parameter derivatives:
elseif isempty(nx) && length(np)==1 && isempty(v),
    par_eps=par;
    eps=abs_eps_p1+rel_eps_p1*abs(par(np));
    par_eps(np)=par(np)+eps;
    J=(sys_rhs(xx,par_eps)-sys_rhs(xx,par))./eps;
    J=reshape(J,n,1,nvec);
    % second order state derivatives:
elseif length(nx)==2 && isempty(np) && ~isempty(v),
    J=zeros(n,n,nvec);
    on=ones(n,1);
    v=reshape(v,[n,1,nvec]);
    if nvec>1
        mult=@(a,b)VAopX(a,b,'*');
    else
        mult=@(a,b)a*b;
    end
    J0=mult(df_deriv(funcs,xx,par,nx(1),[],[],hjac),v);
    for j=1:n
        J(:,j,:)=J0;
        xx_eps=xx;
        eps=abs_eps_x2+rel_eps_x2*abs(xx(j,nx(2)+1,:));
        xx_eps(j,nx(2)+1,:)=xx(j,nx(2)+1,:)+eps;
        Je=mult(df_deriv(funcs,xx_eps,par,nx(1),[],[],hjac),v);
        J(:,j,:)=(Je-J0)./eps(on,1,:);
    end
% mixed state parameter derivatives:
elseif length(nx)==1 && length(np)==1 && isempty(v),
    J=df_deriv(funcs,xx,par,nx(1),[],[],hjac);
    par_eps=par;
    eps=abs_eps_p2+rel_eps_p2*abs(par(np));
    par_eps(np)=par(np)+eps;
    J=(df_deriv(funcs,xx,par_eps,nx(1),[],[],hjac)-J)/eps;
end;

if isempty(J)
    sv=num2str(size(v));
    error(['SYS_DERI: requested derivative nx=%d, np=%d, size(v)=(',sv,...
        ')  does not exist!',nx,np,sv]);
end;
J=reshape(J,[size(J,1),size(J,2),vecdim]);
end

