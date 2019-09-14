function sysvals=psol_sysvals(funcs,xx,par,free_par,varargin)
%% call right-hand sides and derivatives and derivatives of delays in all points listed in xx
% extracted by Jan Sieber from DDE-Biftool's psol_jac to enable vectorisation
%
% $Id: psol_sysvals.m 19 2014-04-11 14:15:36Z jan.sieber $
%
default={'fdim',size(xx,1)};
options=dde_set_options(default,varargin,'pass_on');
%% use functions from funcs
sys_rhs=funcs.sys_rhs;
sys_deri=funcs.sys_deri;
%% problem dimensions
n=size(xx,1);
d=size(xx,2);
nvec=size(xx,3);
%% pre-allocate needed temporary arrays
f=zeros(options.fdim,nvec);
dfdx=zeros(options.fdim,n,d,nvec);  %df/d[xx[k]] for k=1..d+1
dfdp=zeros(options.fdim,length(free_par),nvec);  %df/dp for k=1..d+1
%% check if vectorization has to be replaced by for loop
if nvec>1 && (~isfield(funcs,'x_vectorized')||~funcs.x_vectorized)
    slices=num2cell(1:nvec);
else
    slices={1:nvec};
end
for i=1:length(slices)
    ind=slices{i};
    %% compute f(x,x_tau):
    f(:,ind)=sys_rhs(xx(:,:,ind),par);
    %% determine df/dp:
    for p_i=1:length(free_par)
        dfdp(:,p_i,ind)=sys_deri(xx(:,:,ind),par,[],free_par(p_i),[]);
    end
    %% precompute all Jacobians for delayed values after delayed values are known
    for t_i=1:d
        dfdx(:,:,t_i,ind)=sys_deri(xx(:,:,ind),par,t_i-1,[],[]);
    end
end
sysvals=struct('f',f,'dfdx',dfdx,'dfdp',dfdp);
end
