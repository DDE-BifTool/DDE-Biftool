%% Single-variable derivative using numerical finite difference
%
% $$ y=\frac{d^m}{dt^m}f(x+tv,p) $$
%
% arguments: |rhs=f|, |m=dord|, |x=xx|, |p=par|, |v=xxdev|. Additonal
% arguments for numerical differentiation: |h| interval length, |co| order
% of interpolation polynomial (multiple of 4, |co>m|).
%
% $Id: df_chebnum.m 151 2017-02-14 11:38:46Z jansieber $
%%
function [y,y2]=df_chebnum(rhs,dord,xx,par,xxdev,h,co)
%% Input
%
% *|rhs| user-defined (vectorized) r.h.s, y=rhs(xx,par), which gets
% differentiated
% * |dord| order of requested derivative
% * |xx| argument to be differentiated
% * |par| system parameters
% * |xxdev| direction of derivative
% * |h| length of interval on which Chebyshev interpolation polynomial is
% computed
% * |co| order of interpolation polynomial
%
%% Output
%
% * |y| derivative of order |dord|
% * |y2| derivative of order |dord|, computed at half-order co/2 of
% polynomial
% 
%% Chebyshev interpolation polynom
% |tdev| are the Chebyshev points (2nd type) on [-1,1] (|tdevh| are the
% interpolation points in [-h,h]), Matrix |D| is the differentiation matrix
% for the 2nd-type Chebyshev points: if |y=[y_1,...,y_co]^T| are the
% function values at |xdev|, then |D*y| are the derivatives at |xdev|.
% |D{1}| is for order |co|, |D{2}| is for order |co/2|.
co=ceil(co/4)*4;  % order of interpolation (adj to multiple of 4)
[D{1},tdev]=dde_dchebary(co);
D{2}=dde_dchebary(co/2);
tdevh=tdev*h;
ncalls=size(xxdev,3);
%% Compute interpolation polynomial on [-h,h] for rhs
% This is vectorized. The dimensions of the rhs argument are |n x ntau+1 x
% ncalls x co+1|, (n is called xdim, co+1=nt), then reshaped.
% Output fv has format |n x ncalls x (co+1)|, so needs to be reshaped
% before multiplying it with D^m.
od=ones(1,ncalls);
nt=length(tdevh);
ot=ones(1,nt);
xdim=size(xx,1);
ntaup1=size(xx,2);
ox=ones(1,xdim);
otau=ones(1,ntaup1);
tdevh=reshape(tdevh,[1,1,1,nt]);
xxdev_scale=max(max(abs(reshape(xxdev,xdim*ntaup1,ncalls)),1));
xxdev_scale_vec=reshape(xxdev_scale,[1,1,ncalls]);
xxdev=xxdev./xxdev_scale_vec(ox,otau,:);
xxdt=xx(:,:,od,ot)+tdevh(ox,otau,od,:).*xxdev(:,:,:,ot);
xxdt=reshape(xxdt,[xdim,ntaup1,ncalls*nt]);
fv=rhs(xxdt,par);
fv=reshape(fv,[xdim*ncalls,nt]);
df=fv*(D{1}'^dord)/h^dord;
df2=fv(:,1:2:end)*(D{2}'^dord)/h^dord;
%% Output is derivative at midpoint of interpolation interval
% Combine derivatives of single-vaiable functions using factors output from
% |polarization_coeffs|.
y=reshape(df(:,(nt+1)/2),xdim,ncalls);
y2=reshape(df2(:,(nt+3)/4),xdim,ncalls);
y=y.*xxdev_scale(ox,:).^dord;
y2=y2.*xxdev_scale(ox,:).^dord;
end
%% Matrix corresponding to derivative of interpolant
function [D,x,w]=dde_dchebary(n)
j=0:n;
x=cos(j*pi/n);
on=ones(1,n+1);
w=on;
jodd=mod(j,2)==1;
w(jodd)=-w(jodd);
w([1,end])=0.5*w([1,end]);
wrep=w(on,:);
xreph=x(on,:);
denom=(xreph-xreph');
denom(1:n+2:end)=Inf;
D=-wrep./wrep'./denom;
Drsum=sum(D,2);
D(1:n+2:end)=-Drsum;
end
