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
% $Id: nmfm_mfderiv.m 315 2019-01-29 19:42:21Z jansieber $
%%
order=length(devs);
%% check if parameter deviations involve delays
%if isempty(funcs.sys_mfderi)
funcs=nmfm_checktau(funcs,point,devs);
%end
%% for constant delay, no parameter deviations and user-provided sys_mfderi, 
% (backward compatibility) call sys_mfderi
if ~isempty(funcs.sys_mfderi) && ~funcs.tp_del
    [success,y]=nmfm_try_mfderi(funcs,point,devs);
    if success
        y2=y;
        h=0;
        return
    end
end
%% Otherwise, apply polarization, 
% first find out if some devs are equal
groups=nmfm_dev_group(devs);
%% real coefficients and possibly complex factors for polarization identity
[pdevs,factors]=nmfm_pol_from_devs(devs,order,groups);
%% if no derivatives are provided for rhs (or delay) use pure finite difference
if funcs.sys_dirderi_provided==0 || (funcs.tp_del && funcs.sys_dirdtau_provided==0)
    rhs_scalar=@(h)nmfm_ftau_dev(funcs,point,pdevs,factors,h);
    [y,y2,h]=nmfm_dfdx_scalar(rhs_scalar,order,varargin{:});
    return
end
%% if derivatives are provided, fill up with finite differences
if funcs.sys_dirderi_provided<order
    derivbase=@(xx,p,dxx,dp)funcs.sys_dirderi{funcs.sys_dirderi_provided}(xx,p,dxx,dp);
    for ord=funcs.sys_dirderi_provided+1:order
        funcs.sys_dirderi{ord}=...
            @(x,p,dx,dp)dde_dirderiv(derivbase,...
            x,p,dx,dp,ord-funcs.sys_dirderi_provided,'nf',size(x,1),varargin{:});
    end
end
if funcs.tp_del && funcs.sys_dirdtau_provided<order
    derivbase=@(it,xx,p,dxx,dp)funcs.sys_dirdtau{funcs.sys_dirdtau_provided}(it,xx,p,dxx,dp);
    for ord=funcs.sys_dirdtau_provided+1:order
        funcs.sys_dirdtau{ord}=@(itau,x,p,dx,dp)dde_dirderiv(...
                @(xa,pa,dx,dp)derivbase(itau,xa,pa,dx,dp),...
                x,p,dx,dp,ord-funcs.sys_dirdtau_provided,...
                'nf',1,varargin{:});
    end
end
%% apply finite differences or user-provided sys_dirderi & sys_dirdtau
[y,y2,h]=nmfm_ftau_diff(funcs,order,point,pdevs,factors,varargin{:});
end
