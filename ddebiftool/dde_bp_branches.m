function [t1,t2]=dde_bp_branches(Jfunc,u0,v0,w0,primtan,h)
%% find tangents to branches in branch point of func(u)=0
%
% Inputs:
%
% * Jfunc: Jfunc(u) is jacobian if f(u) (n x n+1 dim)
% * u0: branch point (n+1 x 1 dim)
% * v0: kernel of Jfunc(u0) (n+1 x 2 dim)
% * w0: kernel of Jfunc(u0)' (n x1 dim)
% * primtan: approximate tangent to primary branch, t1 will be tangent for
% other branch, t2 for primary
% * h: deviation for finite difference of 2nd derivative of func()
%
% Outputs:
%
% * t1: tangent to secondary branch
% * t2: tangent most aligned to primary branch
%
%* branches are approx h -> u0+tj*h +O(h^2) (j=1,2)
%
% $Id: dde_bp_branches.m 308 2018-10-28 15:08:12Z jansieber $
%%
J0=Jfunc(u0);
%% coefficients of bifurcation equation
%
% [J0,w0;v0',[0;0]]*[d2beta_ij;d2gamma_ij]=-d2f(u0)v0_i v0_j where the
% solution expansion is v0*a+beta(a) where gamma(a)=0 (and a is 2d, beta
% orthogonal to v0, gamma(a) is 1d).
%
% The bifurcation equation gamma(a)=0 (single equation, 2 variables) has
% both first derivatives equal to zero, s owe need the second derivatives
% g_ij.
Jext=[J0,w0;v0',sparse(2,1)];
rhs=[Jfunc(u0+h*v0(:,1))*v0-Jfunc(u0-h*v0(:,1))*v0,...
    Jfunc(u0+h*v0(:,2))*v0-Jfunc(u0-h*v0(:,2))*v0];
d2f_vivj=-Jext\cat(1,rhs,zeros(2,4));
g_ij=reshape(d2f_vivj(end,:),2,2);
%% solve quadratic bifurcation equation
[Bg,lamg]=eig(g_ij);
lamg=diag(lamg);
assert(prod(lamg)<0);
[~,ix]=sort(lamg);
Bg=Bg(:,ix);
sg=sqrt(diag([-1,1])*lamg(ix));
tpair1=Bg*[sg(2)*[1,-1]; sg(1)*[1,1]];
%% the columns of tpair are the tangents to the branches, rescaled
tpair1_p=tpair1';
tpair=v0*tpair1_p;
%% find which column is closer to primtan
fac=tpair\(primtan/norm(primtan,'inf'));
[~,ix]=sort(abs(fac));
t1=tpair(:,ix(1));
t1=t1/norm(t1,'inf');
t2=tpair(:,ix(2));
t2=t2/norm(t2,'inf');
end