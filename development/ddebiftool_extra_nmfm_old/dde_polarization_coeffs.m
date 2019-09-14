%% Determine coefficients for polarization formula to compute higher-order derivatives
% The combinations in |mats| and |signs| help to reduce an arbitrary-order
% mixed derivative to a linear combination of single-directional
% derivatives.
%%
function [mats,factors]=dde_polarization_coeffs(order)
%% input
%
% * |order|: coefficients for mixed derivatives for degree |n=order| are
% computed
%
%% Output
%
% * |mats|: n x (2^(n-1)-1) matrix P of 1s and 0s (logical)
% * |factors|: 1 x (2^(n-1)-1) array S
%
%% Usage of output
% to obtain multilinear form B(x1,x2,...,xn), sum up S(j)*B(yj,yj,...,yj)
% over j where j=1..length(S),  and yj=P(1,j)*x1+...+P(n,j)*xn.
%%
mats=logical([1,0]);
factors=[1,-1];
for i=2:order
    nend=size(mats,2);
    mats=[mats,mats;true(1,nend),false(1,nend)]; %#ok<AGROW>
    factors=[factors,-factors]; %#ok<AGROW>
end
mats=mats(:,1:end-1);
factors=factors(1:end-1)/factorial(order);
end
