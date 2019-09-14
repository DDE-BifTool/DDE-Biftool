function [y,y2,h]=nmfm_mfderiv(funcs,point,devs,varargin)
%% compute higher-order derivative of r.h.s with deviation arguments devs
%
% y=nmfm_deriv(funs,point,devs)
%
% treats the derivative as a derivative of the abstract fuctional f(x_t) in
% a function space ($f: C([-taumax,0];R^n)->R^n$). The base point x_t is an
% equilibrium such that x_t=point.x=xeq is constant. devs is an array of
% history functions of length k and the output y in R^n is 
% $d^k f(xeq)[devs(1),...,devs(k)]$
%
% For nd constant delays the user may provide sys_mfderi returning the kth
% derivative d^k sys_rhs(xx,devs) where xx has nd+1 equal columns and devs
% is a cell array with k cells of (nd+1) x k arrays of deviations.
% Otherwise, the user may provide sys_dirderi, a directional high-order
% derivative in one direction.
%
%% Inputs
%
% * funcs: user-provided problem functions
% * point: equilibrium point (of type stst or local bifurcation)
% * devs: array of history function structures.
% optional inputs are passed on to subroutines.
%% Outputs
%
% * y (n x 1) array derivative
% * y2 (n x 1) another estimate for the derivative (may be different if
% finite difference quotients are used)
% * h interploation spacing used for finite differencing
%
% $Id$
%%
order=length(devs);
%% for constant delay and user-provided sys_mfderi, call sys_mfderi
if ~funcs.tp_del && funcs.sys_mfderi_provided
    ind_tau=funcs.sys_tau();
    xx=point.x(:,ones(1,length(ind_tau)+1));
    tau0=[0,point.parameter(ind_tau)];
    for i=order:-1:1
        devvals{i}=nmfm_dev_call(devs(i),-tau0);
    end
    y=funcs.sys_mfderi(xx,point.parameter,devvals{:});
    y2=y;
    h=0;
    return
end
%% Otherwise, apply polarization, 
% first find out if some devs are equal
groups=nmfm_dev_group(devs);
%% real coefficients and possibly complex factors for polarization identity
[pdevs,factors]=nmfm_pol_from_devs(devs,order,groups);
%% apply finite differences or user-provided sys_dirderi
[y,y2,h]=nmfm_ftau_diff(funcs,order,point,pdevs,factors,varargin{:});
end
