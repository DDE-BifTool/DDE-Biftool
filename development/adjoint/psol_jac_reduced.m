function [J,res,tT,extmesh,Jstruc]=psol_jac_reduced(funcs,psol,free_par,varargin)
psol_prof=psol.profile;
period=psol.period;
par=psol.parameter;
%% optional
default={'wrapJ',true,'Dtmat',eye(size(psol_prof,1)),'period',true};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
%% define problem & problem dimensions
n=size(psol.profile,1);      % dimension of x
nf=size(options.Dtmat,1); % dimension of f
np=length(free_par);      % number of free parameters
%% constant delay
if funcs.tp_del
    n_tau=funcs.sys_tau(); % delay numbers
    tau=par(n_tau);  % delay values
    d=length(n_tau)+1; % number of delays
    tau=[0;tau(:)]';  % incl tau=0
else
    d=funcs.sys_ntau();
    tau=0;
end
pt=psol;
extmesh=pt.mesh;
W=cell(1,funcs.sys_ntau()+1);
Ws=W;
for t_i=1:d
    [W{t_i},Ws{t_i},neqs]=dde_coeffmat_tau(tau(:,t_i),psol,pass_on{:});
    tau(:,1)=zeros(neqs,1);
    if ~options.wrapJ && length(extmesh)<length(pt.mesh)
        pt=Ws{t_i}.pt;
        extmesh=pt.mesh;
        W{t_i}=Ws{t_i}.W;
    end
    if t_1==1
        xx=NaN(n,d,neqs);
    end
    xx(:,t_i,:)=reshape(W{t_i}{1}*pt.profile(:),n,neqs);
    if funcs.tp_del && t_i<d
        tau(:,t_i+1)=sys_tau(funcs,t_i,xx(:,1:t_i,:),pt.paramaeter);
    end
end
%% check if Jacobian should be wrapped around (mod[0,1])
if options.wrapJ
    options.period=false;
    free_par=[];
end
tT=tau/psol.period;
%% init J, res:
Jstruc.cols.x=1:n:size(pt.profile,2);
Jstruc.cols.xrg=1:Jstruc.cols.x(end)+n-1;
Jstruc.rows.de=1:size(W{1},1);
Jstruc.cols.nx=n;
Jstruc.rows.nf=nf;
Jstruc.cols.np=np;
Jstruc.free_par=free_par;
[~,par_is_delay]=ismember(free_par,n_tau);
[~,delay_is_par]=ismember([0,n_tau],free_par);
Jstruc.delay_is_par=delay_is_par;
Jstruc.par_is_delay=par_is_delay;
if options.period
    Jde_dT=zeros(nf,neqs);
    Jstruc.cols.T=Jstruc.cols.xrg(end)+1;
    Jstruc.cols.nT=1;
else
    Jstruc.cols.T=[];
    Jstruc.cols.nT=0;
end
if np>0
    Jde_dp=zeros(nf,neqs,np);
    Jstruc.cols.p=Jstruc.cols.xrg(end)+Jstruc.cols.nT+(1:np);
else
    Jstruc.cols.p=[];
    Jstruc.cols.np=0;
end
%% obtain all values of sys_rhs, sys_deriv and sys_dxtau
xx=reshape(W{1}*psol.profile(:),n,d,neqs);
dtx=reshape(W{2}*psol.profile(:),n,d,neqs);
vals=psol_sysvals(funcs,xx,par,free_par,'fdim',nf);
%% assemble Jacobian
Dtmatext=sparse([options.Dtmat,repmat(zeros(size(options.Dtmat)),1,n_tau)]);
Dtmat=repmat({Dtmatext},1,neqs);
Dtmat=blkdiag(Dtmat{:});
Jde_dx=Dtmat*W{2};
res=Jde_dx*psol.profile(:)-psol.period*vals.f(:);
dfdx=num2cell(reshape(vals.dfdx,nf,n*d,[]),[1,2]);
dfdx=sparse(blkdiag(dfdx{:}));
Jde_dx=Jde_dx-period*dfdx*W{1};
%% derivative for period if requested
if Jstruc.cols.nT>0
    %% add -f for dT in J:
    tT_ext=num2cell(repmat(sparse(tau/period),1,neqs));
    tT_mat=kron(blkdiag(tT_ext{:}),eye(n));
    Jde_dT=-vals.f(:)-dfdx*tT_mat*dtx(:);
end
%% derivative for parameters if present
if Jstruc.cols.np>0
    Jde_dp=-period*reshape(permute(vals.dfdp,[1,3,2]),neqs*nf,[]);
    for t_i=1:d
        if Jstruc.cols.np>0 && delay_is_par(t_i)>0
            tT_ext=num2cell(repmat(sparse(1,t_i,1,1,d),1,neqs));
            tT_mat=kron(blkdiag(tT_ext{:}),eye(n));
            Jde_dp(:,delay_is_par(t_i))=Jde_dp(:,delay_is_par(t_i))+dfdx*tT_mat*dtx(:);
        end
    end
end
%% reshape and combine
J=Jde_dx;
if Jstruc.cols.nT>0
    J=[J,Jde_dT(:)];
end
if Jstruc.cols.np
    J=[J,reshape(Jde_dp,[nf*neqs,np])];
end
end
%% wrapper around sys_tau
function tau=sys_tau_vec(funcs,ti,xx,par)
if ~funcs.x_vectorized
    for i=size(xx,3):-1:1
        tau(i)=funcs.sys_tau(ti,xx(:,1:ti),par);
    end
else
    tau=funcs.sys_tau(ti,xx(:,1:ti,:),par);
end
    tau=tau(:);
end