function [p,q,sg] = nmfm_border(A, p0, q0)
%% Construct null vectors by bordering technique
% INPUT:
%	 A: n by n matrix
%   p0: approximate null vector s.t. p0*A = 0
%   q0: approximate null vector s.t. A*q0 = 0
% OUTPUT:
%   p: null vector s.t. p*A = 0
%   q: null vector s.t. A*q = 0
%
% $Id: nmfm_border.m 309 2018-10-28 19:02:42Z jansieber $
%
%%
n = size(A,1);

if isempty(p0) && isempty(q0) % No previous null vectors specified
    [U,S,V]=svd(A); %#ok<ASGLU>
    pp0=U(:,end);
    qq0=V(:,end);
elseif isempty(p0) && ~isempty(q0)  % only right previous null vector specified
    pp0=q0;
    qq0=q0;
else
   pp0 = p0';
   qq0 = q0;
end
% Use bordering technique
Border = [A, pp0(:,1); qq0(:,1)', 0];
qsol = Border\[zeros(n,1) ;1];
q = qsol(1:n);
psol = [zeros(1,n), 1]/Border;
p = psol(1:n);
sg=rank(Border)<n;
end

