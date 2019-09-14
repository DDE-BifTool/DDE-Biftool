function [y,cfpol]=nmfm_dirmf_combined(funcs,x0,p0,dev,order,varargin)
%% compute derivative of order for composite functional in single direction
% at equilibrium
%
% [d^k/dh^k] F(x+hv,p+hq)=[d^k/dh^k] f(y(h),p+hq)
%
% assume y=sum yj h^j/j!, j=0..k, determine y first:
%
% yj=j [d^(j-1)/dh^(j-1)] v(-tau(y(h),p+hq))
%
% y0=x(-tau)=pt.x (equilibrium)
% y1=v(-tau)
%
% $Id: nmfm_dirmf_combined.m 309 2018-10-28 19:02:42Z jansieber $
%%
default={'print',0,'cfpol',[]};
options=dde_set_options(default,varargin,'pass_on');
% coeffs for chain rule for up to maxorder
if isempty(options.cfpol)
    options.cf=dde_chainrule_combinatorics(order);
    for k=length(options.cf):-1:1
        cfpol{k}=dde_pol_from_chain(options.cf(k));
    end
else
    cfpol=options.cfpol;
end
ntau=funcs.sys_ntau();
ntaup1=ntau+1;
ndim=size(x0,1);
nvec=length(dev);
npvec=size(p0,3);
assert(size(x0,2)==nvec);
assert(size(p0,1)==1 && (npvec==1 || npvec==nvec));
%% compute delays at equilibrium
x0=repmat(reshape(x0,ndim,1,nvec),[1,ntaup1,1]);
taus=zeros(1,ntaup1,nvec);
np=size(p0,2);
for i=2:ntaup1
    taus(1,i,:)=funcs.sys_tau(i-1,x0(:,1:i-1,:),p0);
end
% determine deviations v of solutions and their derivatives in time, and
% deviation q of parameters. v(j,k,l) is (l-1)st deriv of v_j at tau_k
% (starting to count tau at tau_1=0).
nz=ndim*ntaup1;
v=zeros(ndim,ntaup1,nvec,order+1);
q=NaN(np,nvec);
z=zeros(nz,nvec,order+1);
z(1:nz,:,1)=reshape(x0,nz,nvec);
z(nz+(1:np),:,1)=reshape(repmat(p0,[1,1,nvec/npvec]),np,nvec);
for l=1:nvec
    for i=order+1:-1:1
        v0=real(nmfm_dev_call(dev(l),-taus(1,:,l),i-1));
        v(:,:,l,i)=v0(1:ndim,:);
    end
    q(:,l)=v0(ndim+(1:np),1);
end
z(1:nz,:,2)=reshape(v(:,:,:,1),nz,nvec);
z(nz+(1:np),:,2)=reshape(q,np,nvec);
dtauvals=cat(3,reshape(taus,ntaup1,nvec,1), zeros(ntaup1,nvec,order-1));
f_dtau=@(order,y0,dy)dtaufunc(order,y0,dy,funcs.sys_dirdtau,ndim,ntaup1,nvec,npvec);
for i=2:order
    %% compute d^i/dh^i [tau(y(h),p+hq]
    dtauvals(:,:,i)=nmfm_chainrule(f_dtau,z(:,:,1:i),i-1,'cf',cfpol{i-1});
    %% compute d^(i+1)/dh^(i+1) y(h)
    dvdti=@(order,tau0,dtau)dvfunc(order,dtau,v);
    z(1:nz,:,i+1)=nmfm_chainrule(dvdti,dtauvals(:,:,1:i),i-1,'cf',cfpol{i-1});
    z(1:nz,:,i+1)=i*z(1:nz,:,i+1);
end
f_drhs=@(order,y0,dy)drhsfunc(order,y0,dy,funcs.sys_dirderi,ndim,ntaup1,nvec,npvec);
y=nmfm_chainrule(f_drhs,z,order,'cf',cfpol{order});
end
%% derivative(s) of v(-tau)
function zval=dvfunc(order,dtau,v)
ndim=size(v,1);
sgn=mod(order-1,2)*2-1;
dtau=reshape(dtau,[1,size(dtau)]);
zval=sgn*v(:,:,:,order+1).*(dtau(ones(1,ndim),:,:)).^order;
zval=reshape(zval,[],size(v,3));
end
%% derivative(s) of sys_tau
function dtauval=dtaufunc(order,y0,dy,dtau,ndim,ntaup1,nvec,nvecp)
x0=reshape(y0(1:ndim*ntaup1,:),ndim,ntaup1,nvec);
p0=reshape(y0(ndim*ntaup1+1:end,1:nvecp),1,[],nvecp);
v =reshape(dy(1:ndim*ntaup1,:),size(x0));
q =reshape(dy(ndim*ntaup1+1:end,1:nvecp),size(p0));
dtauval=zeros(ntaup1,nvec);
for i=1:ntaup1-1
    dtauval(i+1,:)=dtau{order}(i,x0(:,1:i,:),p0,v(:,1:i,:),q);
end
end
%% derivative(s) of sys_rhs
function drhsval=drhsfunc(order,y0,dy,rhsdirderi,ndim,ntaup1,nvec,nvecp)
x0=reshape(y0(1:ndim*ntaup1,:),ndim,ntaup1,nvec);
p0=reshape(y0(ndim*ntaup1+1:end,:),1,[],nvecp);
v =reshape(dy(1:ndim*ntaup1,:),size(x0));
q =reshape(dy(ndim*ntaup1+1:end,:),size(p0));
drhsval=rhsdirderi{order}(x0,p0,v,q);
end
