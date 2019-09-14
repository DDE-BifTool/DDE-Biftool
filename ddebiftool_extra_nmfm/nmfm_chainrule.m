function [zchain,cf]=nmfm_chainrule(df,dg,order,varargin)
%% apply chainrule to f(g(x0+t*dx)) to degree order. 
% calls of df are of the form z=df(k,y0,dy) where 
%
% $$z=\partial^k f(y_0)dy^k$$
%
% that is, df returns the kth derivative of f in y_0 in direction dy for
% k=1..order. dg is
%
% * either length(x) -by- (order+1) array where dg(:,k+1) is
%   (pre-computed) kth derivative of dg in x0, wrt to dx, or
% * a function: dg(k,x0,dx) returns partial^k g(x0,dx). For k=0, dg should
%  be (or return) g(x0).
%
% options.cf(k) are polarization coefficients as obtained by
% dde_pol_from_chain, if provided.
%
% $Id: nmfm_chainrule.m 309 2018-10-28 19:02:42Z jansieber $
%
%%
default={'cf',[],'x0',[],'dx',[]};
options=dde_set_options(default,varargin,'pass_on');
%% if polarization coefficients are not provided, re-compute them
if isempty(options.cf)
    cf0=dde_chainrule_combinatorics(order);
    cf=dde_pol_from_chain(cf0(end));
else
    cf=options.cf;
end
if ~isempty(options.x0)&&~isempty(options.dx)
    %% compute directional derivatives of g
    [ndim,nvec]=size(options.x0);
    for k=order:-1:0
        ga(:,:,k+1)=dg(k,options.x0,options.dx);
    end
elseif isnumeric(dg) && size(dg,3)==order+1
    ndim=size(dg,1);
    nvec=size(dg,2);
    ga=dg;
else
    error('nmfm_chainrule:call',...
        'nmfm_chainrule: inconsistent optional arguments');
end
zchain=0;
for k=1:order
    c=cf(k);
    cm=reshape(c.mats,[1,1,size(c.mats,1),size(c.mats,2)]);
    ydev=ga(:,:,2:size(cm,3)+1);
    cm=repmat(cm,[ndim,nvec,1,1]);
    for j=1:length(c.fac)
        ypol=sum(ydev.*cm(:,:,:,j),3);
        zchain=zchain+df(c.gx,ga(:,:,1),ypol)*c.fac(j);
    end
end    
end
