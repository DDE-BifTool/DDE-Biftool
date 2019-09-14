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
% $Id: nmfm_ftau_diff.m 174 2017-03-10 21:59:17Z jansieber $
%%
function [y,y2,h]=nmfm_ftau_diff(funcs,order,pt,devs,fac,varargin)
ndev=length(fac);
if ~funcs.tp_del
    dirderi=funcs.sys_dirderi;
else
    dirderi=funcs.sys_dirmf;
end
if ~isempty(dirderi)
    %% user has provided directional derivatives
    for i=ndev:-1:1
        if ~funcs.tp_del
            taus = [0,pt.parameter(funcs.sys_tau())];
            xdevs=real(nmfm_dev_call(devs(i),-taus));
            xeq=pt.x(:,ones(1,length(taus)));
            yi(:,i)=dirderi(order,xeq,pt.parameter,xdevs,zeros(size(pt.parameter)));
        else
            xdev=@(j,theta)real(nmfm_dev_call(devs(i),theta,j));
            yi(:,i)=dirderi(order,pt.x,pt.parameter,xdev);
        end
    end
    facvec=repmat(reshape(fac,1,ndev),size(yi,1),1);
    y=sum(yi.*facvec,2);
    y2=y;
    h=0;
else
    %% single-deviation derivative by finite difference
    rhs_scalar=@(h)nmfm_ftau_dev_sav(funcs,pt,devs,fac,h);
    [y,y2,h]=nmfm_dfdx_scalar(rhs_scalar,order,varargin{:});
end
end