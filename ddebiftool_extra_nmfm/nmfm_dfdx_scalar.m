%% Single-variable derivative using numerical finite difference
%
% $$ y=\frac{d^m}{dx^m}f(x)\vert_{x=0} $$
%
% arguments: |rhs=f|, |m=dord|. Additonal
% arguments for numerical differentiation:
% |fac| (1.99) ratio between distances of interpolation points from 0:
% interpolation points are:
% -fac^(n-1)*h,...-fac*h,-h,0,h,h*fac,...,h*fac^(n-1),
% |h| (default 1) scaling in spacing and output matrix
%
% $Id: nmfm_dfdx_scalar.m 309 2018-10-28 19:02:42Z jansieber $
%% Input
%
% * |rhs_inp| r.h.s, y=f(x), which gets differentiated
% * |dord| order of requested derivative
%
%% Optional inputs
%
% * |h| (default |eps^(1/(dord+1))|) initial scaling in spacing and output matrix
% * |fac| (1.99) ratio between successive deviations
% * |isvectorized| can rhs be called with several x?
% * |maxit| (20) maximal number of adjustments of h by fac
% * |output| (1): switch to 2 to swap order of first two outputs  
%
%% Output
%
% * |y| derivative of order |dord|
% * |y2| derivative of order |dord|, lower order estimate
% * |h| spacing/scaling used in final estimate
%%
function [y,y2,h]=nmfm_dfdx_scalar(rhs_inp,dord,varargin)
% process optional arguments
default={'h',eps^(1/(dord+1)),'fac',1.99,'isvectorized',true,'maxit',20,'output',1};
options=dde_set_options(default,varargin,'pass_on');
%% approximation order is 2,5,6,9,10,... for derivatives of order 1,2,3,4,5,..
fac=options.fac;
deriv=nmfm_Dbary(dord,fac);
deriv.x=deriv.x*options.h;
% wrapper permitting to assume vectorization
rhs=@(x)rhs_vec(rhs_inp,x,options.isvectorized); 
% apply linear combination to inerpolants to obtain approx derivative
D=@(ind,deriv,h)(deriv.f0*deriv.D0+deriv.fvals(:,2*ind+(1:length(deriv.D)))*deriv.D)/h^dord;
%% Compute initially 3 approximations for derivative for h/fac, h h*fac
% We first check which two approximations are closer to each other.
% Depending on that we decrease (or increase) h interatively until the
% closeness relation changes.
deriv.f0=rhs(0);
szfv=size(deriv.f0);
deriv.x=[...
    deriv.x(:,1)/fac, deriv.x, deriv.x(:,end)*fac];
deriv.fvals=rhs(deriv.x);
h=options.h;
df=[D(0,deriv,h/fac),D(1,deriv,h),D(2,deriv,h*fac)];
err=max(abs(diff(df,[],2)),[],1);
err_dir=err(1)<err(2);
for it=1:options.maxit
    if err_dir ~=(err(1)<err(2))
        break
    end
    if err(1)<err(2)
        deriv.x=[deriv.x(:,1)/options.fac,deriv.x(:,1:end-1)];
        deriv.fvals=[rhs(deriv.x(:,1)),deriv.fvals(:,1:end-2)];
        h=h/fac;
    elseif err(1)>=err(2)
        deriv.x=[deriv.x(:,2:end),deriv.x(:,end)*options.fac];
        deriv.fvals=[deriv.fvals(:,3:end),rhs(deriv.x(:,end))];
        h=h*fac;
    end
    df=[D(0,deriv,h/fac),D(1,deriv,h),D(2,deriv,h*fac)];
    err=max(abs(diff(df,[],2)),[],1);
end
%% Final output
% apply Richardson extrapolation to the approximations h,h*fac
% "lower order" estimate is either the approximation for h/fac or h*fac,
% depending on which one is further away from the final output.
rfac=1/fac^deriv.approx_order;
extrapmat=[[1;-rfac;0], [0;1;-rfac]]/(1-rfac);
df_extrap=df*extrapmat;
y=reshape(df_extrap(:,2),szfv);
err1=norm(df_extrap(:,2)-df(:,1),inf);
err3=norm(df_extrap(:,2)-df(:,3),inf);
if err1>err3
    y2=reshape(df(:,1),szfv);
else
    y2=reshape(df(:,3),szfv);
end
if options.output==2
    ytmp=y;
    y=y2;
    y2=ytmp;
end
end
%% wrapper for vectorized call to rhs
function f=rhs_vec(rhs,x,isvec)
if isvec
    f=rhs(x(:)');
    f=reshape(f,[],numel(x));
else
    for i=numel(x):-1:1
        ftmp=rhs(x(i));
        f(:,i)=ftmp(:);
    end
end
end
