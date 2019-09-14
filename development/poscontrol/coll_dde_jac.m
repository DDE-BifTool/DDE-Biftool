function [J,res,tT,extmesh]=coll_dde_jac(funcs,pt,free_par,varargin)
%% return Jacobian and residual of collocation problem
% function [J,res,tT,extmesh]=coll_dde_jac(funcs,pt,free_par,...)
%% INPUT:
% *  funcs problem functions
% * pt point structure (psol or hcli)
% * free_par indices of free parameters numbers in N^np
% 
%% OUTPUT:
%	J jacobian in R^(n*deg*l+n+s x n*deg*l+1+n+np)
%	res residual in R^(n*deg*l+n+s)
%   tT delays, scaled by period
%   extmesh mesh of time points, extended back to -max(tT(:))
%
% modified from psol_jac to permit arbitrary nesting of state-dependent
% delays, vectorisation and re-use for computation of Floquet multipliers,
% extra conditions for state-dependent delay (p_tau, etc), connecting
% orbits
%
%% Optional inputs:
%
% * c collocation parameters in [0,1]^m
% * If 'wrapJ' (default true) the returned jacobian is augmented with
% derivative wrt period and free_par and wrapped around periodically inside
% [0,1] if ~wrapJ the returned jacobian puts its entries into the interval
%  [-min(delay,0),max(1-[0,delays])], no augmentation is done. This
%  Jacobian can be used for computation of Floquet multipliers and modes
%
% If 'c_is_tvals' (default false) then the values in c (must be not empty) are
% taken as those values of time in [0,1] of the entire period, where
% residual and Jacobian are calculated.
%
% The argument 'Dtmat' (default=Id) is multiplied with the time
% derivative. Dtmat==zeros() permits algebraic equations.
%
% The argument 'bc' (default true) controls whether boundary conditions
% should be attached.
%
% $Id$
%%
%% optional (args are also passed on to coll_dde_map)
default={'period',true,'Dtmat',[]};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
%% get values at collocation points
[coll,extmesh]=coll_dde_map(funcs,pt,'nderivs',1,pass_on{:});
%% define problem & problem dimensions
par=pt.parameter;
[n,neqs]=size(coll.profile{1,1}); % dimension & number of collocation points
oneqs=ones(1,neqs);
d=size(coll.profile,1);           % number of delays + 1 (counting delay 0)
%% which parameters are delays?
if ~funcs.tp_del
    delay_is_par=ismember([0,funcs.sys_tau()],free_par);
end
F=coll_sysvals(funcs,coll.profile,par,free_par);
%% pre-factor of derivative if required
if isempty(options.Dtmat)
    options.Dtmat=eye(F.nf);
end
Dts=sparse_blkdiag(options.Dtmat(:,:,oneqs));
%% assemble residual
res=Dts*coll.profile{1,2}(:)-pt.period*F.f;
%% assemble Jacobian
dX=coll.jac(:,1);  % stages for derivatives wrt x
dT=cell(d,1);      % stages for derivatives wrt T
dP=cell(d,1);      % stages for derivatives wrt parameters
dP{1}=sparse(n*neqs,length(free_par));
dT{1}=sparse(n*neqs,1);
Jde_dx=Dts*coll.jac{1,2};
Jde_dp=-pt.period*reshape(F.dfdp,neqs*F.nf,length(free_par));
Jde_dT=-F.f(:);
yp=coll.profile(:,2);
for t_i=1:d
    Jde_dx=Jde_dx-pt.period*F.dfdx{t_i}*dX{t_i};
    Jde_dT=Jde_dT-pt.period*F.dfdx{t_i}*dT{t_i};
    if funcs.tp_del
        Jde_dp=Jde_dp-pt.period*F.dfdx{t_i}*dP{t_i};
    elseif delay_is_par(t_i)>0
        Jde_dp(:,delay_is_par(t_i))=...
            Jde_dp(:,delay_is_par(t_i)) + F.dfdx{t_i}*yp{t_i}(:);
    end
    if t_i==d
        break
    end
    ypti_T= sparse_blkdiag(reshape(yp{t_i+1}/pt.period,1,1,n*neqs));
    tau_by_T=repmat(coll.tau(:,t_i+1)',n,1)/pt.period;
    dT{t_i+1}=ypti_T*tau_by_T(:);
    if funcs.tp_del
        dtau_dp=reshape(repmat(reshape(F.dtau_dp{t_i+1},...
            1,neqs,length(free_par)),[n,1,1]),n*neqs,length(free_par));
        dP{t_i+1}=-ypti_T*dtau_dp;
        dX=d_accu(t_i+1,dX,F.dtau_dx,ypti_T,n);
        dT=d_accu(t_i+1,dT,F.dtau_dx,ypti_T,n);
        dP=d_accu(t_i+1,dP,F.dtau_dx,ypti_T,n);
    end
end
%% combine matrices
J=Jde_dx;
if options.period
    J=[J,Jde_dT];
end
J=[J,Jde_dp];
tT=coll.tau/pt.period;
end
%% accummulate dependence of dX, dT, dP through past delays
function D=d_accu(t_i,D,dtau_dx,ypti_T,n)
sum_tau_D=sparse(size(dtau_dx{t_i,1},1),size(D{1},2));
for t_k=1:t_i-1
    sum_tau_D=sum_tau_D+dtau_dx{t_i,t_k}*D{t_k};
end
sum_tau_D=kron(sum_tau_D,ones(n,1));
D{t_i}=D{t_i}-ypti_T*sum_tau_D;
end
