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
% $Id: nmfm_ftau_diff.m 309 2018-10-28 19:02:42Z jansieber $
%%
function [y,y2,h]=nmfm_ftau_diff(funcs,order,pt,devs,fac,varargin)
default={'debug',false};
options=dde_set_options(default,varargin,'pass_on');
ndev=length(fac);
ndim=size(pt.x,1);
p0=pt.parameter(1,:,ones(1,ndev));
if ~funcs.tp_del
        taus = [0,pt.parameter(funcs.sys_tau())];
        x0=pt.x(:,ones(1,length(taus)),ones(1,ndev));
        for i=ndev:-1:1
            devvals(:,:,i)=real(nmfm_dev_call(devs(i),-taus));
        end
        xdevs=devvals(1:ndim,:,:);
        pdevs=reshape(devvals(ndim+1:end,1,:),1,[],ndev);
        yi=funcs.sys_dirderi{order}(x0,p0,xdevs,pdevs);
else
    x0=pt.x(:,ones(1,ndev));
    yi=nmfm_dirmf_combined(funcs,x0,p0,devs,order,varargin{:});
    if options.debug
        for i=ndev:-1:1
            f0=@(h)nmfm_ftau_dev(funcs,pt,devs(i),1,h);
            yn(:,i)=nmfm_dfdx_scalar(f0,order);
        end
        assert(all(abs(yi(:)-yn(:))<eps^(1./order/2)));
    end
end
facvec=repmat(reshape(fac,1,ndev),size(yi,1),1);
y=sum(yi.*facvec,2);
y2=y;
h=0;
end