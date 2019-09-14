function df=dde_dirderiv(finp,xx,par,dx,dpar,order,varargin)
%% finite-difference directional derivatives approximations of finp of orbitrary order
%
% |dde_dirderiv(f,x,p,y,q,j)| returns approximately
%
% $$ \frac{\partial^j}{\partial h^j}f(x+hy,p+hq)\vert_{h=0}$
% 
%% Input
%
% * |finp| r.h.s, y=f(x,p), which gets differentiated
% * |xx| (n x ndelays+1 x nvec1 x ...)  states (current and delayed)
% * |par| (1 x npar) parameters
% * |dx| (n x ndelays+1 x nvec1 x ...)  state deviations (current and delayed)
% * |dpar| (1 x npar) parameter deviations
% * |order| order of requested derivative
%
%% Optional inputs
%
% * |hjac| (default 1e-2) initial scaling in spacing and output matrix
% * |fac| (1.99) ratio between successive deviations
% * |x_vectorized| can finp be called with nvec1 >1?
%
%% Output
%
% * |df| approx derivative of order |order| at (xx,par) in direction
% (dx,dpar)
%
% <html>
% $Id$
% </html>
%
%% process options
xsize=size(xx);
n=xsize(1);
default={'nf',n,'x_vectorized',false,'hjac',eps^(1/(order+2))};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
%% wrap for vectorization of finp
if ~options.x_vectorized && size(xx,3)>1
    f=@(x,p)wrap_f(finp,x,p);
else
    f=finp;
end
ndelays=xsize(2);
vecdim=[xsize(3:end),1];
nvec=prod(vecdim);
xx=reshape(xx,[n,ndelays,nvec]);
dx=reshape(dx,[n,ndelays,nvec]);
switch order
    %% use simple formulas for low orders
    case 1
        fp=f(xx+options.hjac*dx,par+options.hjac*dpar);
        fm=f(xx-options.hjac*dx,par-options.hjac*dpar);
        df=0.5*(fp-fm)/options.hjac;
    case 2
        fp=f(xx+options.hjac*dx,par+options.hjac*dpar);
        fm=f(xx-options.hjac*dx,par-options.hjac*dpar);
        f0=f(xx,par);
        df=(fp+fm-2*f0)/options.hjac^2;
    otherwise
        %% assert that dpar is zero
        assert(norm(dpar)==0);
        rhs=@(h)fdevh(f,h,xx,dx,par,options.nf);
        df=nmfm_dfdx_scalar(rhs,order,'isvectorized',true,'h',options.hjac,pass_on{:});
end
df=reshape(df,[options.nf,vecdim]);
end
function y=fdevh(f,h,xx,dx,par,nf)
[n,ndelays,nvec]=size(xx);
nh=length(h);
h=reshape(h,[1,1,1,nh]);
h=reshape(h(ones(n,1),ones(ndelays,1),ones(nvec,1),:),n,ndelays,[]);
xx=reshape(xx(:,:,:,ones(nh,1)),n,ndelays,[]);
dx=reshape(dx(:,:,:,ones(nh,1)),n,ndelays,[]);
y=f(xx+h.*dx,par);
y=reshape(y,[nf*nvec,nh]);
end
function y=wrap_f(f,x,p)
for i=size(x,3):-1:1
    y(:,i)=f(x(:,:,i),p);
end
end