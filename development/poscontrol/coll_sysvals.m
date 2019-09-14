function vals=coll_sysvals(funcs,y,par,free_par)
%% call right-hand sides and derivatives and derivatives of delays in all points listed in xx
% extracted by Jan Sieber from DDE-Biftool's psol_jac to enable
% vectorisation
%
% $Id$
%
if funcs.tp_del
    d=funcs.sys_ntau()+1;
else
    d=length(funcs.sys_tau())+1;
end
[n,neqs]=size(y{1,1});
oneqs=ones(1,neqs);
par_neqs=par(1,:,oneqs);
%% obtain all values of sys_rhs, sys_deriv and sys_dxtau
xxvec=NaN(n,d,neqs);
if funcs.tp_del
    vals.dtau_dx=cell(d,d-1);
end
for t_i=d:-1:1
    xxvec(:,t_i,:)=reshape(y{t_i,1},n,1,neqs);
end
%% compute derivatives (assuming vectorization)
vals.f=funcs.sys_rhs(xxvec,par_neqs);
vals.nf=size(vals.f,1);
vals.f=vals.f(:);
%% derivatives wrt x of rhs and tau
for t_i=d:-1:1
    vals.dfdx{t_i}=funcs.sys_deri(xxvec,par_neqs,t_i-1,[],[]);
    vals.dfdx{t_i}=sparse_blkdiag(reshape(vals.dfdx{t_i},vals.nf,n,[]));
    if funcs.tp_del && t_i>1
        for t_k=1:t_i-1
            vals.dtau_dx{t_i,t_k}=...
                funcs.sys_dtau(t_i-1,xxvec(:,1:t_i-1,:),par_neqs,t_k-1,[]);
            vals.dtau_dx{t_i,t_k}=sparse_blkdiag(reshape(vals.dtau_dx{t_i,t_k},1,n,[]));
        end
    end
end
%% derivatives wrt parameters of rhs and tau
vals.dfdp=zeros(vals.nf,neqs,length(free_par));
vals.dtau_dp=repmat({zeros(1,neqs,length(free_par))},1,d);
for p_i=1:length(free_par)
    vals.dfdp(:,:,p_i)=funcs.sys_deri(xxvec,par_neqs,[],free_par(p_i),[]);
end
if funcs.tp_del
    for t_i=d:-1:2
        dtau_dp=NaN(length(free_par),neqs);
        for p_i=1:length(free_par)
            dtau_dp(p_i,:)=funcs.sys_dtau(t_i-1,xxvec(:,1:t_i-1,:),par_neqs,[],free_par(p_i));
        end
        vals.dtau_dp{t_i}=dtau_dp';
    end
end
end