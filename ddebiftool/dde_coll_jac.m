function [J,res,tT,extmesh]=dde_coll_jac(funcs,pt,free_par,varargin)
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
% The argument
%
% $Id: dde_coll_jac.m 336 2019-05-09 17:04:13Z jansieber $
%%
%% optional (args are also passed on to coll_dde_map)
default={'Dtmat',[],'output','J','matrix','sparse'};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
%% are delayed arguments derivatives?
if isfield(funcs,'delayed_derivs')
    deriv_order=[0,funcs.delayed_derivs];
    max_deriv=max(deriv_order)+1;
else
    deriv_order=[];
    max_deriv=1;
end
%% get values at collocation points
[coll,extmesh]=dde_coll_map(funcs,pt,'nderivs',max_deriv,'free_par_ind',free_par,...
    'deriv_order',deriv_order,pass_on{:});
%% define problem & problem dimensions
par=pt.parameter;
[n,neqs]=size(coll.profile{1,1}); % dimension & number of collocation points
oneqs=ones(1,neqs);
d=size(coll.profile,1);           % number of delays + 1 (counting delay 0)
tT=coll.tau/pt.period;            % rescaled delays
if isempty(deriv_order)
    deriv_order=zeros(1,d);
end
%% evaluate  r.h.s. and delays at requested collocation points
[F,coll]=dde_coll_sysvals(funcs,coll,par,free_par,deriv_order,pt.period,options.output);
%% pre-factor of derivative if required
if isempty(options.Dtmat)
    options.Dtmat=eye(F.nf);
end
Dts=sparse_blkdiag(options.Dtmat(:,:,oneqs));
%% assemble residual & return if only residual is required
res=Dts*coll.profile{1,2}(:)-F.f;
if strcmp(options.output,'res')
    J=res;
    return
end
%% assemble Jacobian
% select interpolation matrices and their derivatives
bdiag=@(n1,vec)sparse_blkdiag(reshape(vec,n1,1,[]));
xind=@(x,ord,dim)cat(dim,x{sub2ind(size(coll.jac),1:size(coll.jac,1),ord+1)});
dXk=xind(coll.jac,deriv_order,1);
yk =        reshape(xind(coll.profile,deriv_order,  3),[],1);
yk1=bdiag(n,reshape(xind(coll.profile,deriv_order+1,3),[],1));
if funcs.tp_del
    y1=bdiag(n,reshape(cat(3,coll.profile{:,2}),[],1));
    Mat=speye(n*neqs*d)+y1*F.dtau_dx;
    [M.L,M.U,M.P,M.Q,M.R]=lu(Mat);
    dX0=linsolve(M,cat(1,coll.jac{:,1}));
    dXk=dXk-yk1*F.dtau_dx*dX0;
end
%% Derivative wrt profile
J.profile=Dts*coll.jac{1,2}-F.dfdx*dXk;
if ~strcmp(options.matrix,'sparse')
    J.profile=full(J.profile);
end
assert(size(J.profile,1)==length(res));
if strcmp(options.output,'profile')
    J=J.profile;
    return
end
%% derivative wrt parameter
dP=yk1*F.dtau_dp;
if funcs.tp_del
    dP0=linsolve(M,-y1*F.dtau_dp);
    dP=dP+yk1*F.dtau_dx*dP0;
end
J.parameter=-F.dfdp+F.dfdx*dP;
%% derivative wrt period
dT=(bdiag(n*neqs,yk)*deriv_order(:)-yk1*F.tau)/pt.period;
if funcs.tp_del
    dT0=linsolve(M,dT);
    dT=dT-yk1*F.dtau_dx*dT0;
end
J.period=-Dts*coll.profile{1,2}(:)/pt.period+F.dfdx*dT;
assert(...
    size(J.period,1)==length(res) && ...
    size(J.parameter,1)==length(res))
end
%%
function x=linsolve(M,y)
    x=M.Q*(M.U\(M.L\(M.P*(M.R\y))));
end