function [W,W_unwrap_struc,neqs]=dde_coeffmat(n_tau,pt,varargin)
%% optional
default={'unwrapJ',false,'c',[],'c_is_tvals',false,'Dtmat',eye(size(pt.profile,1)),'kron',true,...
    'nderivs',2};
options=dde_set_options(default,varargin,'pass_on');
nderivs=options.nderivs+1;
W_unwrap_struc=struct();
%% define problem & problem dimensions
n=size(pt.profile,1);      % dimension of x
%nf=size(options.Dtmat,1); % dimension of f
tau=pt.parameter(n_tau);  % delay values
tau=[0;tau(:)];  % incl tau=0
d=length(n_tau)+1; % number of delays
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
Pb=repmat({zeros(d,pt.degree+1,neqs)},1,nderivs);  
res={repmat({zeros(d,neqs,size(pt.profile,2))},1,nderivs)};
c_tau=zeros(d,neqs); % time points t-tau(i-1) for current colloc point t
c_tau_mod=zeros(d,neqs); % time points t-tau(i-1) for current colloc point t (mod[0,1])
index_b=zeros(d,neqs); % starting indices of collocation intervals for t-tau(i-1)
index_b_mod=zeros(d,neqs); % starting indices of collocation intervals for t-tau(i-1) (mod[0,1])
h_int_del=zeros(d,neqs); % length of collocation interval in which t-tau(i-1) lies
c_tau_trans=zeros(d,neqs); % position of t-tau(i) in its colloc interval, rescaled to [0,1]
tT=tau(:)/pt.period;
tT=repmat(tT,1,neqs);
polycoeffs={@poly_elg,@poly_del,@poly_d2l};
for t_i=1:d
    c_tau(t_i,:)=t_c-tT(t_i,:);
    c_tau_mod(t_i,:)=mod(c_tau(t_i,:),1);
    [~,ibase,c_tau_trans(t_i,:)]=psol_eva(pt.mesh,pt.mesh,c_tau_mod(t_i,:),pt.degree); 
    index_b_mod(t_i,:)=(ibase-1)*pt.degree+1;
    index_b(t_i,:)=index_b_mod(t_i,:)+round(c_tau(t_i,:)-c_tau_mod(t_i,:))*nint*pt.degree;
    h_int_del(t_i,:)=h_int(ibase);
    for j=1:nderivs
        Pb{j}(t_i,:,:)=polycoeffs{j}(pt.degree,c_tau_trans(t_i,:))./...
            h_int_del(t_i*ones(pt.degree+1,1),:).^(j-1);
    end
end
if options.unwrapJ
    indmin=min(index_b(:));
    indshift=1-indmin;
    indmax=max(index_b(:))+pt.degree;
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
[t_ind,k_ind,m_ind]=ndgrid(1:d,1:neqs,0:pt.degree);
for l=1:nout
    for j=1:nderivs
        Pb{j}=reshape(permute(Pb{j},[1,3,2]),[],1);
        irow=t_ind+(k_ind-1)*d;
        icol=ind{l}(t_ind+(k_ind-1)*d)+m_ind;
        res{l}{j}=sparse(irow(:),icol(:),Pb{j},neqs*d,ntst{l});
        if options.kron
            res{l}{j}=kron(res{l}{j},eye(n));
        end
    end
end
if options.unwrapJ
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

