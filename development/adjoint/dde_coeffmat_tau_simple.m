function [W,neqs]=dde_coeffmat_tau_simple(tau,pt,varargin)
%% Returns matrix W, W' (and W'') mapping base points to collocation points 
% W{i}{j} is mapping of base points to (i-1)th derivative at time t-tau(j)
% where counting of tau starts with delay tau(1)=0.
%% optional parameters
% c pertmis t ospecify collocation points inside collocation interval,
% scaled to [0,1] (if c_is_tvals is false). If c_is_tvals is true then the
% vector c is treated as the requested time points on the entire interval.
% 'kron' switches on/off if W should be expanded to system size, 'nderivs'
% sets whether W, W' and W'' should be computed.
default={'c',[],'c_is_tvals',false,'kron',true,'nderivs',2};
options=dde_set_options(default,varargin,'pass_on');
nderivs=options.nderivs+1;
%% define problem & problem dimensions
n=size(pt.profile,1);      % dimension of x
nint=(length(pt.mesh)-1)/pt.degree; % number of collocation intervals
if ~options.c_is_tvals
    neqs=pt.degree*nint;    % number of equations from collocation points
else
    neqs=length(options.c); % only evaluate res,Jac at times c
end
tcoarse=pt.mesh(1:pt.degree:end); % boundaries of collocation intervals
h_int=diff(tcoarse); % lengths of collocation intervals
%% obtain collocation parameters and times
if ~options.c_is_tvals
    if isempty(options.c)
        c=poly_gau(pt.degree);
        c=c(:)';
    else
        c=options.c(:)';
    end
    t_c=tcoarse(ones(length(c),1),1:end-1)+...
        c(ones(length(h_int),1),:)'.*h_int(ones(length(c),1),:);
else
    t_c=options.c;
end
t_c=t_c(:)'; % evaluation points for DE residuals
%% pre-allocate needed temporary arrays
% Lagrange stamp, mapping delayed values,1st and 2nd derivs onto profile
Pb=repmat({zeros(pt.degree+1,neqs)},1,nderivs);  
W={repmat({zeros(neqs,size(pt.profile,2))},1,nderivs)};
polycoeffs={@poly_elg,@poly_del,@poly_d2l};
%% compute where delayed x are located
c_tau=t_c-tau/pt.period; % time points t-tau/T for current colloc point t
c_tau_mod=mod(c_tau,1);% time points t-tau/T for current colloc point t (mod[0,1])
% find position of t-tau/T in its colloc interval, rescaled to [0,1]
[~,ibase,c_tau_trans]=psol_eva(pt.mesh,pt.mesh,c_tau_mod,pt.degree);
% starting indices of collocation intervals for t-tau (mod[0,1])
index_b_mod=(ibase-1)*pt.degree+1;
% length of collocation interval in which t-tau lies
h_int_del=h_int(ibase);
for j=1:nderivs
    Pb{j}=polycoeffs{j}(pt.degree,c_tau_trans)./...
        h_int_del(ones(pt.degree+1,1),:).^(j-1);
end
[irow,m_ind]=ndgrid(1:neqs,0:pt.degree);
for j=1:nderivs
    Pb{j}=reshape(Pb{j}',[],1);
    icol=index_b_mod(irow)+m_ind;
    W{j}=sparse(irow(:),icol(:),Pb{j},neqs,length(pt.mesh));
    if options.kron
        W{j}=kron(W{j},eye(n));
    end
end
neqs=size(W{1},1)/n;
end

