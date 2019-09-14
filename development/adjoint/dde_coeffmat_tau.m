function [W,W_unwrap_struc,neqs]=dde_coeffmat_tau(tau,pt,varargin)
%% optional
default={'wrapJ',true,'c',[],'c_is_tvals',false,'Dtmat',eye(size(pt.profile,1)),'kron',true,...
    'nderivs',2};
options=dde_set_options(default,varargin,'pass_on');
nderivs=options.nderivs+1;
W_unwrap_struc=struct();
%% define problem & problem dimensions
n=size(pt.profile,1);      % dimension of x
%nf=size(options.Dtmat,1); % dimension of f
nint=(length(pt.mesh)-1)/pt.degree; % number of collocation intervals
if ~options.c_is_tvals
    neqs=pt.degree*nint;          % number of equations from collocation points
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
res={repmat({zeros(neqs,size(pt.profile,2))},1,nderivs)};
polycoeffs={@poly_elg,@poly_del,@poly_d2l};
%% compute where delayed x are located
c_tau=t_c-tau/period; % time points t-tau/T for current colloc point t
c_tau_mod=mod(c_tau,1);% time points t-tau/T for current colloc point t (mod[0,1])
% find position of t-tau/T in its colloc interval, rescaled to [0,1]
[~,ibase,c_tau_trans]=psol_eva(pt.mesh,pt.mesh,c_tau_mod,pt.degree);
% starting indices of collocation intervals for t-tau (mod[0,1])
index_b_mod=(ibase-1)*pt.degree+1;
% starting indices of collocation intervals for t-tau
index_b=index_b_mod+round(c_tau-c_tau_mod)*nint*pt.degree;
% length of collocation interval in which t-tau lies
h_int_del=h_int(ibase);
for j=1:nderivs
    Pb{j}=polycoeffs{j}(pt.degree,c_tau_trans)./...
        h_int_del(ones(pt.degree+1,1),:).^(j-1);
end
if ~options.wrapJ
    indmin=min(index_b);
    indshift=1-indmin;
    indmax=max(index_b)+pt.degree;
    indrg=indmin:indmax;
    t_ind=mod(indrg-1,nint*pt.degree)+1;
    t_shift=floor((indrg-1)/(nint*pt.degree));
    extmesh=pt.mesh(t_ind)+t_shift;
    index_b_shift=index_b+indshift;
    nout=2;
    ind={index_b_mod,index_b_shift};
    res(2)={repmat({zeros(d,neqs,indshift+indmax)},1,nderivs)};
    ntst={size(pt.profile,2),indshift+indmax};
else
    nout=1;
    ind={index_b_mod};
    ntst={size(pt.profile,2)};
end
[k_ind,m_ind]=ndgrid(1:neqs,0:pt.degree);
for l=1:nout
    for j=1:nderivs
        Pb{j}=reshape(Pb{j}',[],1);
        irow=(k_ind-1)*d+1;
        icol=ind{l}(1+(k_ind-1)*d)+m_ind;
        res{l}{j}=sparse(irow(:),icol(:),Pb{j},neqs,ntst{l});
        if options.kron
            res{l}{j}=kron(res{l}{j},eye(n));
        end
    end
end
if ~options.wrapJ
    xx_ext=NaN(size(pt.profile,1),indshift+indmax);
    for k=1:size(xx_ext,2)
        xx_ext(:,k)=pt.profile(:,mod(k-indshift-1,size(pt.profile,2)-1)+1);
    end
    pt_ext=pt;
    pt_ext.mesh=extmesh;
    pt_ext.profile=xx_ext;
    W_unwrap_struc=struct('W',res(2),'point',pt_ext,'ind1',indshift);
end
W=res{1};
end

