%% Matrix corresponding to kth derivative of interpolant around 0
%%
function ret=nmfm_Dbary(order,fac,varargin)
%
%% Inputs:
% * order: integer, order k of derivative
% * fac: spacing of interpolation points is 
%       -fac^(n-1)*h,...-fac*h,-h,0,h,h*fac,...,h*fac^(n-1)
% * optional 'h' (default 1) scaling in spacing and output matrix
% * optional 'n' (default order-1) number of interpolation points to each
% side
%
%% Output structure ret. Fields:
% * D (1 x n) coefficients for non-zero interpolcation points
% * x (2 x n) interpolcation points to each side: x(:,k+1)=h*fac^k*[-1;1]
% * D0: coeficient at 0
% * deriv: order of derivative
% * approx_order: interpolcation order of formula (2*order for odd,
% 2*order+1 for even order derivatives)
%
% $Id: nmfm_Dbary.m 309 2018-10-28 19:02:42Z jansieber $
%
%% optional arguments
default={'h',1,'n',order};
options=dde_set_options(default,varargin,'pass_on');
%% interpolcation points
facpow=fac.^(0:options.n-1);
facpowsgn=[-facpow;facpow];
x=[0,facpowsgn(:)']*options.h;
nx=length(x);
%% create barycentric weights wi
[xi1,xi2]=ndgrid(x,x);
dxi=xi1-xi2+eye(nx);
w=1./prod(dxi,2)';
on=ones(1,nx);
wrep=w(on,:);
xreph=x(on,:);
%% Barycentric formula for 1st derivative at interpolation points
% gives matrix D (see paper by Berrut, Trefethen ,SIAM Review 46(3))
denom=(xreph-xreph');
denom(1:nx+1:end)=Inf;
D=-wrep./wrep'./denom;
Drsum=sum(D,2);
D(1:nx+1:end)=-Drsum;
%% D^order is (higher) order derivative interpolation
D=D^order;
D=D(1,:)';
ret=struct('D',D(2:end),'x',reshape(x(2:end),2,[]),'D0',D(1),...
    'deriv',order,'approx_order',order*2+1);
if mod(order,2)
    ret.D0=0;
    ret.approx_order=ret.approx_order-1;
end
end
