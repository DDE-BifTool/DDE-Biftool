%% Eigenvalues of constant coefficient linear DDE using Chebyshev interpolation
%
% $Id: dde_stst_eig_cheb.m 369 2019-08-27 00:07:02Z jansieber $
%% 
function stability=dde_stst_eig_cheb(A,tau,varargin)
if length(tau)==size(A,3)-1
    taus=[0;tau(:)];
elseif length(tau)==size(A,3)
    taus=tau(:);
else
    error('dde_stst_eig_cheb:delay',...
        'dde_stst_eig_cheb: number of derivatives=%d, but number of delays=%d',...
        size(A,3),length(tau));
end
%% cut delays to non-negative values for now
taus=max(0,taus);
%% set options
[dim,dum,ndel]=size(A); %#ok<ASGLU>
default={'ncheb',[],'minimal_real_part',-1,...
    'root_accuracy',1e-6,'maxsize',200,'inisize',max(2,ceil(40/dim)),'nchebfac',2,...
    'max_number_of_eigenvalues',20,'min_number_of_eigenvalues',[],...
    'scale_w',true,'closest',[],'tauscal_limit',1e-5,...
    'imagthreshold',1e-10,'discard_accuracy_factor',1e5,'lhs_matrix',eye(size(A,1)),...
    'VWout',false};
options=dde_set_options(default,varargin,'pass_on','method');
if isempty(options.inisize)
    options.inisize=max(2,ceil(40/dim));
end
if isempty(options.minimal_real_part)
    options.minimal_real_part=-1;
end
taumax=max(taus(:));
stability=dde_stst_stability();
if options.VWout
    stability.mesh=[];
end
%% estimate scaling of delay problem
c_est=0;
for i=1:ndel
    c_est=c_est+norm(A(:,:,i),'inf')*abs(taus(i));
end
%% special case that all delays are 0 or no delays present
if c_est==0
    itau0=find(taus==0);
    A0=A(:,:,itau0(1));
    for i=2:length(itau0)
        A0=A0+A(:,:,itau0(i));
    end
    [stability.l0,stability.v,stability.w,v,W]=dde_eigsort(A0,options,dim,1); %#ok<ASGLU>
    stability.l0=stability.l0(:).';
    stability.l1=stability.l0(:);
    [res,stability.w]=check_residual(stability.l0,stability.v,stability.w,W,A,taus,options);
    stability.err=max(abs(res),[],1);
    stability=discard_roots(stability,options);
    if options.VWout
        stability.mesh=0;
    end
    return
elseif taumax<options.tauscal_limit
    %% scale only up to tauscal_limit and append  A(ndel+1)=0
    taus=[taus;options.tauscal_limit];
    A=cat(3,A,zeros(dim,dim));
    taumax=options.tauscal_limit;
    scale=2/taumax;
    tauscal=taus*scale;
    Ascal=A/scale;
else
    %% rescale delays and derivatives such that max delay=2
    scale=2/taumax;
    tauscal=taus*scale;
    Ascal=A/scale;
end
%% choose initial number of chebyshebev points
if isempty(options.ncheb)
    ncheb=max(4,norm(A(:,:,1)));
    if isfinite(options.minimal_real_part)
        for i=2:length(taus)
            ncheb=ncheb+norm(A(:,:,i))*exp(abs(options.minimal_real_part)*taus(i));
        end
        ncheb=round(ncheb);
    end
else
    ncheb=options.ncheb;
end
if isempty(options.min_number_of_eigenvalues)
    options.min_number_of_eigenvalues=dim;
end
ncheb=max(ceil(min(ncheb,options.inisize/dim)),ceil(options.min_number_of_eigenvalues/dim));
%% compute approx eigenvalues and gradually double ncheb if error too large
err=Inf;
lambda=NaN(0,1);
while err>options.root_accuracy && ncheb<=options.maxsize/dim
    [D,Etau,x]=dde_dchebary_loc(ncheb,tauscal);
    D=kron(D,eye(dim));
    D(1:dim,:)=0;
    for i=1:ndel
        D(1:dim,:)=D(1:dim,:)+Ascal(:,:,i)*kron(Etau(i,:),eye(dim));
    end
    [lambda,v,w,V,W]=dde_eigsort(D,options,dim,scale);
    %% check residuals
    [res,w,W]=check_residual(lambda,v,w,W,A,taus,options);
    if ~isempty(res)
        err=max(abs(res(isfinite(res(:)))));
    else
        err=0;
    end
    ncheb=round(ncheb*options.nchebfac);
end
stability=dde_stst_stability('l0',lambda(:).','l1',lambda(:),...
    'err',max(abs(res),[],1),'v',v,'w',w);
if options.VWout
    nlam=numel(lambda);
    stability.v=reshape(V,dim,[],nlam);
    stability.w=reshape(W,dim,[],nlam);
    stability.mesh=x/scale;
end
stability=discard_roots(stability,options);
end
%% Matrices corresponding to interpolant at tau and derivative of interpolant 
function [D,Etau,x,w]=dde_dchebary_loc(n,tau)
tau=tau(:);
taumax=max(abs(tau));
% Nodes & weights for Chebyshev interpolation
[x,w]=dde_coll_chebxwt(n,1:n+1,[-1,1]);
w=w(end:-1:1)*(2/taumax)^n;
x=(x(end:-1:1)-1)/2*taumax;
on=ones(1,n+1);
ot=ones(1,length(tau));
%% interpolant for tau
w_by_tau_m_x=w(ot,:)./(-tau(:,on)-x(ot,:));
denw=sum(w_by_tau_m_x,2);
isinf=~isfinite(w_by_tau_m_x(:));
Etau=w_by_tau_m_x./denw(:,on);
Etau(isinf)=1;
%% derivative
wrepn=w(on,:);
xreph=x(on,:);
denom=(xreph-xreph');
denom(1:n+2:end)=Inf;
D=-wrepn./wrepn'./denom;
Drsum=sum(D,2);
D(1:n+2:end)=-Drsum;
end
%% get, sort and scale eigenvalues
function [lambda,v,w,V,W]=dde_eigsort(D,opts,dim,scale)
mp=matrix_or_pair(D,opts.lhs_matrix);
[V,S,W]=eig(mp{:});
lambda=diag(S)*scale;
lambda(~isfinite(lambda))=-Inf;
[dum,iy]=sort(imag(lambda),'descend'); %#ok<ASGLU>
lambda=lambda(iy);
V=V(:,iy);
W=W(:,iy);
if isempty(opts.closest)
    [rlam,ix]=sort(real(lambda),'descend');
    sel=rlam>=opts.minimal_real_part;
else
    [rlam,ix]=sort(abs(lambda-opts.closest),'ascend');
    sel=true(size(rlam));
end    
lambda=lambda(ix);
sel(1:min(length(lambda),opts.min_number_of_eigenvalues))=true;
%% don't discard complex conjugate of last selected eigenvalue (if fewer than max)
if ~sel(end)
    lastsel=find(sel,1,'last');
    if abs(real(lambda(lastsel))-real(lambda(lastsel+1)))<opts.imagthreshold
        sel(lastsel+1)=true;
    end
end
imax=min(opts.max_number_of_eigenvalues,length(sel));
sel(imax+1:end)=false;
lambda=lambda(sel);
V=V(:,ix);
v=V(1:dim,sel);
vn_inv=diag(1./sqrt(sum(v.*conj(v),1)));
v=v*vn_inv;
W=W(:,ix);
w=W(1:dim,sel);
if length(lambda)<opts.min_number_of_eigenvalues
    lambda(end+1:opts.min_number_of_eigenvalues)=-Inf;
    v(:,end+1:opts.min_number_of_eigenvalues)=NaN;
    w(:,end+1:opts.min_number_of_eigenvalues)=NaN;
end
V=V(:,sel)*vn_inv;
W=W(:,sel);
end
%% check residual for eigenvalues
function [res,w,W]=check_residual(lambda,v,w,W,A,taus,options)
nl=length(lambda(~isnan(lambda)));
dim=size(A,1);
res=zeros(dim,nl);
Delta_f=@(lam,args)dde_stst_ch_matrix(A,taus,lam,'lhs_matrix',options.lhs_matrix,args{:});
for i=nl:-1:1
    Delta=Delta_f(lambda(i),{});
    Delta_pv=Delta_f(lambda(i),{'deri',1})*v(:,i);
    if options.scale_w
        wscal=(w(:,i)'*Delta_pv)';
        w(:,i)=w(:,i)/wscal;
        W(:,i)=W(:,i)/wscal;
    end
    res(:,i)=Delta*v(:,i);
end
end
%%
function stability=discard_roots(stability,options)
keep=stability.err<=options.root_accuracy*options.discard_accuracy_factor;
if any(~keep)
    stability.discarded.l0=stability.l0(~keep);
    stability.discarded.err=stability.err(~keep);
end
stability.l0=stability.l0(keep);
stability.l1=stability.l1(keep);
stability.err=stability.err(keep);
if options.VWout
    stability.v=stability.v(:,:,keep);
    stability.w=stability.w(:,:,keep);
else
    stability.v=stability.v(:,keep);
    stability.w=stability.w(:,keep);
end
end
%% create extended lhs matrix
function L=matrix_or_pair(D,lhsmat)
L={D};
if norm(lhsmat-eye(size(lhsmat)),'inf')==0
    return
end
L={D,blkdiag(lhsmat,eye(size(D,1)-size(lhsmat,1)))};
end
