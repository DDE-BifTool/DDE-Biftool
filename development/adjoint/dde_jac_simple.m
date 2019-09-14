function [J,res,Jstruc]=dde_jac_simple(funcs,pt,free_par,varargin)
%% return Jacobian and residual
%% optional
default={'nf',size(pt.profile,1),'period',true};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
%% define problem & problem dimensions
n=size(pt.profile,1);      % dimension of x
nf=options.nf; % dimension of f
np=length(free_par);      % number of free parameters
%% constant delay
n_tau=funcs.sys_tau(); % delay numbers
tau=pt.parameter(n_tau);  % delay values
d=length(n_tau)+1; % number of delays
tau=[0;tau(:)]';  % incl tau=0
%% generate W, W' for each delay
% W{t_i}{1}) is W for delay t_i-1 (t_i=1 corresponds to tau=0)
W=cell(1,d);
for t_i=1:d
    [W{t_i},neqs]=dde_coeffmat_tau_simple(tau(:,t_i),pt,pass_on{:});
end
%% store structure of Jacobian
Jstruc.cols.x=1:n:size(pt.profile,2);
Jstruc.cols.xrg=1:Jstruc.cols.x(end)+n-1;
Jstruc.rows.de=1:size(W{1},1);
Jstruc.cols.nx=n;
Jstruc.rows.nf=nf;
Jstruc.cols.np=np;
Jstruc.free_par=free_par;
%% which parameters are delays?
[~,Jstruc.par_is_delay]=ismember(free_par,n_tau);
[~,Jstruc.delay_is_par]=ismember([0,n_tau],free_par);
if options.period
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
for t_i=d:-1:1
    xxvec(:,:,t_i)=reshape(W{t_i}{1}*pt.profile(:),n,neqs);
    dtx{t_i}=W{t_i}{2}*pt.profile(:);
end
%% compute derivatives (assuming vectorization in x)
f=funcs.sys_rhs(xxvec,pt.parameter);
f=f(:);
for t_i=d:-1:1
    dfdx{t_i}=funcs.sys_deri(permute(xxvec,[1,3,2]),pt.parameter,t_i-1,[],[]);
    dfdx{t_i}=sparse_blkdiag(reshape(dfdx{t_i},nf,n,[]));
end
for p_i=length(free_par):-1:1
    dfdp(:,:,p_i)=funcs.sys_deri(xxvec,pt.parameter,[],free_par(p_i),[]);
end
%% assemble Jacobian
res=dtx{1}-pt.period*f;
Jde_dx=W{1}{2};
for t_i=d:-1:1
    Jde_dx=Jde_dx-pt.period*dfdx{t_i}*W{t_i}{1};
end
%% derivative for period if requested
if Jstruc.cols.nT>0
    %% add -f for dT in J:
    Jde_dT=-f;
    for t_i=2:d
        Jde_dT=Jde_dT-dfdx{t_i}*dtx{t_i}*tau(t_i)/pt.period;
    end
end
%% derivative for parameters if present
if Jstruc.cols.np>0
    Jde_dp=-pt.period*reshape(dfdp,neqs*nf,[]);
    for t_i=1:d % if delay is parameter
        is_tp=Jstruc.delay_is_par(t_i);
        if is_tp>0
            Jde_dp(:,is_tp)=Jde_dp(:,is_tp)+dfdx{t_i}*dtx{t_i};
        end
    end
end
%% combine overall matrix
J=Jde_dx;
if Jstruc.cols.nT>0
    J=[J,Jde_dT(:)];
end
if Jstruc.cols.np
    J=[J,reshape(Jde_dp,[nf*neqs,np])];
end
end
