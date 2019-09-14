%% Evaluate function h->sum cj f(xeq+h yj_t) for functional defined by funcs
%
% function y=nmfm_ftau_dev(funcs,pt,devs,fac,h)
%
%% Inputs 
%
% * funcs: structure containing problem functions
% * pt: constant point to which deviations are added
% * devs (1 x ndev): array of deviation structures
% * fac (1 x ndev): factors cj in weighted sum
% * h (1 x nh): distance of deviation (vectorized)
%
%% Outputs
%
% * y: (ny x nh) array sum cj f(pt.x+h dev(j),pt.parameter)
%
% $Id: nmfm_ftau_dev.m 169 2017-03-03 22:06:16Z jansieber $
%%
function y=nmfm_ftau_dev(funcs,pt,devs,fac,h)
ndev=length(fac);
nh=length(h);
nvec=nh*ndev;
ndim=length(pt.x);
np=length(pt.parameter);
%% find values of history where to compute functional
% done differently if delays are state-dependent or constant
if ~funcs.tp_del % constant delays
   taus = pt.parameter(funcs.sys_tau());
   r = length(taus); % Number of delays
   taus = cat(3,zeros(nh,ndev),repmat(reshape(taus,[1,1,r]),[nh,ndev,1])); % First delay zero
   xx=NaN(ndim,r,nvec);
   par=NaN(1,np,nvec);
   for i=1:r+1
       [xx(:,i,:),par(1,:,:)]=nmfm_dev_callvec(pt.x,pt.parameter,devs,h,-taus(:,:,i));
   end
else % state-dependent delays
   r = funcs.sys_ntau(); % Number of delays
   taus = cat(1,zeros(1,nh,ndev),NaN(r,nh,ndev)); % zeroth delay is zero, others unknown
   theta=reshape(-taus(1,:),nh,ndev);
   [x0,par]=nmfm_dev_callvec(pt.x,pt.parameter,devs,h,theta);
   xx = cat(2,x0,NaN([ndim,r,nvec])); % first entry in xx known, others not
   for i = 1:r
       taus(i+1,:) = funcs.sys_tau(i,xx(:,1:i,:),par(1,:,:));
       theta=reshape(-taus(i+1,:),nh,ndev);
       xx(:,i+1,:)=nmfm_dev_callvec(pt.x,pt.parameter,devs,h,theta);
   end
end
%% evaluate functional sys_rhs in all deviations
yj=NaN(ndim,nvec);
for k=1:nvec
    yj(:,k)=funcs.sys_rhs(xx(:,:,k),par(1,:,k));
end
%% add up functional values, weighted with factors fac (cj)
yj=reshape(yj,[],nh,ndev);
fac=reshape(fac,[1,1,ndev]);
facvec=repmat(fac,[size(yj,1),nh,1]);
y=sum(yj.*facvec,3);
end
%% wrapper around call to dev structures
function [x,p]=nmfm_dev_callvec(x0,p0,devs,h,theta)
ndev=length(devs);
nh=length(h);
ndim=length(x0);
p0=p0(:);
np=length(p0);
% theta is assumed to be nh x ndevs
x=NaN(ndim,nh,ndev);
p=NaN(np,nh,ndev);
for i=1:ndev
    devvals=nmfm_dev_call(devs(i),theta(:,i)');
    xdevs=devvals(1:ndim,:).*h(ones(ndim,1),:);
    pdevs=devvals(ndim+(1:np),:).*h(ones(np,1),:);
    x(:,:,i)=xdevs+x0(:,ones(1,nh));
    p(:,:,i)=pdevs+p0(:,ones(1,nh));
end
x=reshape(real(x),ndim,1,nh*ndev);
p=reshape(real(p),1,np,nh*ndev);
end
