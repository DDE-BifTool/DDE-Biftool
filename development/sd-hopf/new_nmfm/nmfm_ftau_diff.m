%% Differentiate h->sum cj f(xeq+h yj_t) for functional defined by funcs k times
%
% function y=nmfm_ftau_diff(funcs,order,pt,devs,fac)
%
%% Inputs
%
% * funcs: user-provided problem functions, set with set_funcs
% * order: order of derivative to be taken
% * pt: point of kind 'stst' (or local bifurcation), pt.x and pt.parameter
% are needed
% * devs: 1 x ndev array of history functions yj
% * fac; 1 x ndev array of complex numbers cj
%
% devs and fac are obtained from polarization identity in nmfm_mfderiv
% (nmfm_dev_group, nmfm_pol_from_devs) to compute an arbitrary mixed
% derivative.
%
%% Outputs
%
% * y: directional derivative d/dh [sum cj f(x0+h yj)] at h=0
% * y2: lower order estimate (when using finite differences)
% * h: interpolation point spacing used for finte difference approximation
%
% $Id$
%%
function [y,y2,h]=nmfm_ftau_diff(funcs,order,pt,devs,fac,varargin)
ndev=length(fac);
if isempty(funcs.sys_dirderi)
    %% user has provided directional derivatives
    for i=ndev:-1:1
        if ~funcs.tp_del
            taus = [0,par(funcs.sys_tau())];
            xdevs=real(nmfm_dev_call(devs(i),-taus));
            yi(:,i)=funcs.sys_dirderi(order,pt.x,xdevs,pt.parameter);
        else
            xdev=@(j,theta)real(nmfm_dev_call(devs(i),theta,j));
            yi(:,i)=funcs.sys_dirderi(order,pt.x,xdev,pt.parameter);
        end
    end
    facvec=repmat(reshape(fac,1,ndev),size(yi,1),1);
    y=sum(yi.*facvec,2);
    y2=y;
    h=0;
else
    %% single-deviation derivative by finite difference
    rhs_scalar=@(h)nmfm_ftau_dev(funcs,pt,devs,fac,h);
    [y,y2,h]=nmfm_dfdx_scalar(rhs_scalar,order,varargin{:});
end
end