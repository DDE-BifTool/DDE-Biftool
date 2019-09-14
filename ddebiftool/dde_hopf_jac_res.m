function [J,res]=dde_hopf_jac_res(funcs,pt,free_par,method,varargin)
%% Jacobian for Hopf problem
% function [J,res]=dde_hopf_jac_res(funcs,pt,free_par,method,pref,varargin)
% INPUT:
%   funcs problem function
%	x current Hopf solution guess in R^n
%	omega current Hopf frequency guess in R
%	v current eigenvector guess in C^n
%	par current parameters values in R^p
%	free_par free parameter numbers in N^d 
%	c normalization vector of v in C^(1 x n)
% OUTPUT: 
%	J jacobian in R^(3n+2+s x 3n+1+p)
%	res residual in R^(3n+2+s x 1)

% (c) DDE-BIFTOOL v. 2.00, 30/11/2001
%
% $Id: dde_hopf_jac_res.m 369 2019-08-27 00:07:02Z jansieber $
%
%%
default={'pref',repmat(pt,0,1)};
options=dde_set_options(default,varargin,'pass_on');
pref=options.pref;
[ind,len]=dde_ind_from_point(pt,free_par);
[A,taus]=dde_stst_linearize_rhs(funcs,pt);
taus=taus(:).';
m=size(A,3);
om=ones(m,1);
xx=pt.x(:,om);
n=size(xx,1);
n_fp=length(free_par);
z=1i*pt.omega;
eztau=exp(-z*taus);
lhs_mat=funcs.lhs_matrix(n);
D=dde_stst_ch_matrix(A,taus,z,'lhs_matrix',lhs_mat);
res=[...
    funcs.sys_rhs(xx,pt.parameter);...
    real(D*pt.v);
    imag(D*pt.v)];
%% Jacobian residual
J=zeros(length(res),len);
nf=size(A,1);
J(1:nf,ind.x)=sum(A,3);
J(1:nf,ind.parameter)=dde_stst_linearize_rhs(funcs,pt,'free_par',free_par);
dD=dde_stst_ch_matrix(A,taus,z,'deri',1,'lhs_matrix',lhs_mat);
dDvom=1i*dD*pt.v;
df=@(i,j,k)funcs.sys_deri(xx,pt.parameter,i,j,k);
dtau=@(i,j,k)funcs.sys_dtau(i,xx(:,1:i),pt.parameter,j,k);
dDdxv=zeros(n);
for j=0:m-1
    for k=0:m-1
        dDdxv=dDdxv-df([k j],[],real(pt.v))*eztau(k+1);
        dDdxv=dDdxv-1i*df([k j],[],imag(pt.v))*eztau(k+1);
        if funcs.tp_del~=0 && k>0 && j<k
            dDdxv=dDdxv+z*A(:,:,k+1)*eztau(k+1)*(pt.v*dtau(k,j,[]));
        end
    end
end
% 2nd deriv wrt parameters
if ~funcs.tp_del
    n_tau=funcs.sys_tau();
end
dDdpv=zeros(n,n_fp);
for k=1:n_fp
    for j=0:m-1
        dDdp=df(j,free_par(k),[]);
        dDdpv(:,k)=dDdpv(:,k)-dDdp*pt.v*eztau(j+1);
        if funcs.tp_del==0 && j>0
            if j>=1 && n_tau(j)==free_par(k)
                dDdpv(:,k)=dDdpv(:,k)+z*A(:,:,j+1)*pt.v*eztau(j+1);
            end
        elseif funcs.tp_del~=0 && j>0
            d=dtau(j,[],free_par(k)); % d=d(tau(j))/d(free_par(k))
            dDdpv(:,k)=dDdpv(:,k)+z*A(:,:,j+1)*pt.v*eztau(j+1)*d;
        end
    end
end
%% assemble jacobian
ic=               [ind.x(:);   ind.v.re(:); ind.v.im(:); ind.omega;  ind.parameter(:)];
J(nf+(1:nf),  ic)=[real(dDdxv), real(D),    -imag(D),    real(dDvom), real(dDdpv)];
J(2*nf+(1:nf),ic)=[imag(dDdxv), imag(D),     real(D),    imag(dDvom), imag(dDdpv)];
%% append normalization condition
resv=pt.v'*pt.v-1;
Jcv=zeros(1,len);
Jcv(1,[ind.v.re(:); ind.v.im(:)])=[2*real(pt.v)',  2*imag(pt.v)'];
%% if pref is non-empty minimize phase difference
if ~isempty(pref) &&  ~strcmp(method.preprocess,'dde_jac2square_preprocess')
    resph=imag(pref.v'*pt.v);
    Jph=zeros(1,len);
    Jph(1,[ind.v.re(:),ind.v.im(:)])=[-imag(pref.v)',real(pref.v)'];
else
    resph=zeros(0,1);
    Jph=zeros(0,len);
end
res=[res;resv;resph];
J=[J;Jcv;Jph];
end
