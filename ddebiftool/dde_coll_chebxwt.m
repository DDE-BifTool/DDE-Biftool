%% Nodes & weights for Chebyshev interpolation, scaled to [0,1] by default
%%
function [x,w]=dde_coll_chebxwt(n,sel,interval)
%
% inputs: n=degree, sel=which of the nodes and weights should be output
% (default: all of them), interval=what base interval should be used
% (default [0,1] 
% output: x=node positions, w=barycentric weights
% $Id: dde_coll_chebxwt.m 369 2019-08-27 00:07:02Z jansieber $
%%
if nargin<2
    sel=1:n+1;
end
if nargin<3
    interval=[0,1];
end
len=interval(2)-interval(1);
j=0:n;
x=cos(j*pi/n);
w=ones(1,n+1);
jodd=mod(j,2)==1;
w(jodd)=-w(jodd);
w([1,end])=0.5*w([1,end]);
x=(x+1)/2*len+interval(1);
w=w*(2/len)^n;
x=x(end:-1:1);
w=w(end:-1:1);
x=x(sel);
w=w(sel);
end
