function Delta=ch_matrix(funcs,xx,par,lambda,varargin)
%% Characteristic matrix and its derivatives
% function Delta=ch_matrx(funcs,x,par,l)
% INPUT:
%   funcs problem functions
%	x steady state solution in R^n (either n x (ntau+1), or n x 1)
%	par parameter values
%	lamba complex number at which  charactersitic matrix is computed
%   optional named argument: 'deri', integer (default 0) return derivative
%   of characteristic matrix wrt lambda
% OUTPUT: 
%	D characteristic matrix in C^(n x n)
%
% |xx| is assumed to be equilibrium such that |xx(:,2:end)| are ignored and
% |xx(:,ones(1,(ntau+1))| is used instead.
%
% (c) DDE-BIFTOOL v. 2.00, 23/11/2001
%
% $Id: ch_matrix.m 369 2019-08-27 00:07:02Z jansieber $
%
%%
[A,tau]=dde_stst_linearize_rhs(funcs,dde_stst_create('x',xx,'parameter',par));
Delta=dde_stst_ch_matrix(A,tau,lambda,...
    'lhs_matrix',funcs.lhs_matrix(size(A,1)),varargin{:});
end
