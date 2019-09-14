function BINV = nmfm_binv(funcs, xx, par, lambda, q, p, zeta, kappa)
%% Compute a bordered inverse
% INPUT:
%   xx: state vector
%   par: parameter vector
%	lambda: eigenvalue of the characteristic matrix Delta
%   q: normalized null vector of Delta
%   p: normalized null vector of Delta^T
%   zeta: n-vector
%   kappa: scalar
% OUTPUT:
%	BINV: solution *function* to the bordered inverse problem
%
% $Id: nmfm_binv.m 79 2015-01-02 18:42:50Z jan.sieber $
%
%%
DDelta = ch_matrix(funcs,xx,par,lambda,'deri',1);

A = [ch_matrix(funcs,xx,par,lambda), q; p, 0];
B = [zeta + kappa*DDelta*q; 0];
X = A\B;
xi = [X(1,1); X(2,1)];
gamma = -p*DDelta*xi + (1/2)*kappa*p*ch_matrix(funcs,xx,par,lambda,'deri',2)*q;

BINV = @(theta) exp(lambda*theta)*(xi+gamma*q - kappa*theta*q);


end

