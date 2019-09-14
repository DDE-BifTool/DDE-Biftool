function [J,res,tT,extmesh]=psol_jac_combined(funcs,pt,free_par,varargin)
%% residual & Jacobian of collocation problem for periodic orbits
% function [J,res,tT,extmesh]=psol_jac(c,T,profile,t,deg,par,free_par,phase,varargin)
% INPUT:
%   funcs problem functions
%	c collocation parameters in [0,1]^m
%   pt: point with fields mesh,period, profile, degree, parameter
%	free_par free parameters numbers in N^np
%	phase use phase condition or not (s = 1 or 0)
%   wrapJ (optional key-value pair, default true) 
%        wrap time points periodically into [0,1]
% OUTPUT:
%	J jacobian in R^(n*deg*l+n+s x n*deg*l+1+n+np)
%	res residual in R^(n*deg*l+n+s)
%   tT delays, scaled by period
%   extmesh mesh of time points, extended back to -max(tT(:))
%
% modified to permit arbitrary nesting of state-dependent delays,
% vectorisation and optional re-use for computation of Floquet multipliers
% and extra conditions for state-dependent delay (p_tau, etc)
%
% Optional inputs:
%
% If 'wrapJ' (default true) the returned jacobian is augmented with
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

% $Id$
%
% 
%
%% optional
default={'bc',true,'rotationcheck',true,'ph',true,'extra_condition',repmat(struct(),1,0),...
    'reference',pt};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
[J,res,tT,extmesh]=coll_dde_jac(funcs,pt,free_par,pass_on{:});
%% periodicity condition:
ind=dde_ind_from_psol(pt,free_par);
if options.bc
    resbc=pt.profile(:,1)-pt.profile(:,end);
    if options.rotationcheck
        resbc=mod(resbc+pi,2*pi)-pi;
    end
    res=[res;resbc(:)];
    n=length(resbc);
    Jbc=sparse(n,size(J,2));
    Jbc(:,ind.profile(:,1))=eye(n);
    Jbc(:,ind.profile(:,end))=-eye(n);
    J=[J;Jbc];
end
%% phase condition:
if options.ph 
    p0=p_axpy(0,pt,[]);
    [resph,p_Jph_dx]=p_dot(p0,pt,'derivatives',[0,1]);
    Jph_dx=sparse(1,size(J,2));
    Jph_dx(1:ind.profile(end))=p_Jph_dx.profile(:)';
    res=[res;resph];
    J=[J;Jph_dx];
end
end