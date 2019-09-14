function [vals,coll]=dde_coll_sysvals(funcs,coll,par,free_par,deriv_order,T,output)
%% call right-hand sides and derivatives and derivatives of delays in all points listed in xx
% extracted by Jan Sieber from DDE-Biftool's psol_jac to enable
% vectorisation
%
% $Id: dde_coll_sysvals.m 336 2019-05-09 17:04:13Z jansieber $
%
%% number of delays (incl 0)
d=length(deriv_order);
[n,neqs]=size(coll.profile{1,1});
if funcs.p_vectorized
    oneqs=ones(1,neqs);
    par_neqs=par(1,:,oneqs);
else
    par_neqs=par;
end
%% rescale profiles and Jacobians by period
%% Rescale derivatives by period
max_deriv=size(coll.profile,2)-1;
for k=0:max_deriv
    scal=T^k;
    for t_i=1:d
        coll.jac{t_i,k+1}=coll.jac{t_i,k+1}/scal;
        coll.profile{t_i,k+1}=coll.profile{t_i,k+1}/scal;
    end
end
%% obtain all values of sys_rhs, sys_deriv and sys_dxtau
xfvec=zeros(n,d,neqs);
for t_i=d:-1:1
    xfvec(:,t_i,:)=reshape(coll.profile{t_i,deriv_order(t_i)+1},n,1,neqs);
end
%% compute r.h.s.
vals.f=funcs.sys_rhs(xfvec,par_neqs);
vals.nf=size(vals.f,1);
vals.f=vals.f(:);
vals.tau=reshape(coll.tau,neqs*d,1);
if strcmp(output,'res')
    return
end
%% compute derivatives (assuming vectorization)
bdiag=@(n1,M)sparse_blkdiag(reshape(M,n1,n,[]));
%% derivatives wrt x of rhs and tau
for t_i=d:-1:1
    vals.dfdx{t_i}=bdiag(vals.nf,funcs.sys_deri(xfvec,par_neqs,t_i-1,[],[]));
end
vals.dfdx=cat(2,vals.dfdx{:});
%% derivatives wrt parameters of rhs
np=length(free_par);
vals.dfdp=zeros(vals.nf,neqs,np);
for p_i=1:np
    vals.dfdp(:,:,p_i)=funcs.sys_deri(xfvec,par_neqs,[],free_par(p_i),[]);
end
vals.dfdp=reshape(vals.dfdp,vals.nf*neqs,np);
%% return for constant delay case
if ~funcs.tp_del
    %% which parameters are delays?
    [delay_is_par,delay_pos]=ismember([0,funcs.sys_tau()],free_par);
    delay_pos=delay_pos(delay_is_par);
    dtauc_dp=full(sparse(find(delay_is_par),delay_pos,ones(length(delay_pos),1),d,np));
    vals.dtau_dp=reshape(repmat(reshape(dtauc_dp,1,d,np),neqs,1,1),neqs*d,np);
    return
end
%% state dependent delays
% derivatives of tau wrt x
%% obtain all values of sys_rhs, sys_deriv and sys_dxtau
xtauvec=zeros(n,d,neqs);
for t_i=d:-1:1
    xtauvec(:,t_i,:)=reshape(coll.profile{t_i,1},n,1,neqs);
end
vals.dtau_dx=repmat({sparse(neqs,n*neqs)},d,d);
for t_i=d:-1:2
    for t_k=1:t_i-1
        vals.dtau_dx{t_i,t_k}=...
            bdiag(1,funcs.sys_dtau(t_i-1,xtauvec(:,1:t_i-1,:),par_neqs,t_k-1,[]));
    end
end
vals.dtau_dx=cell2mat(vals.dtau_dx);%reshape(cat(1,vals.dtau_dx{:}),neqs*d,n*neqs*d);
%% derivatives of tau wrt parameters
vals.dtau_dp=zeros(neqs,d,np);
for t_i=2:d
    for p_i=1:np
        vals.dtau_dp(:,t_i,p_i)=reshape(...
            funcs.sys_dtau(t_i-1,xtauvec(:,1:t_i-1,:),par_neqs,[],free_par(p_i)),...
            neqs,[]);
    end
end
vals.dtau_dp=reshape(vals.dtau_dp,neqs*d,np);
end
