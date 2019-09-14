%% Determine coefficients for polarization formula to compute higher-order derivatives
% The combinations in |mats| and |signs| help to reduce an arbitrary-order
% mixed derivative to a linear combination of single-directional
% derivatives.
%%
function [mats,factors]=dde_polarization_coeffs(order,groups)
%% input
%
% * |order|: coefficients for mixed derivatives for degree |n=order| are
% computed
% * |groups|: cell array, forming a partition of 1:order (eg, {1:2,3} to
% indicate that linear deviation arguments 1 and 2 are equal)
%
%% Output
%
% * |mats|: n x nd matrix P of 1s and 0s or -1s
% * |factors|: 1 x nd array S
%
% if groups is missing then nd=2^n-1, otherwise, it can be shorter
%
%% Usage of output
% to obtain multilinear form B(x1,x2,...,xn), sum up S(j)*B(yj,yj,...,yj)
% over j where j=1..length(S),  and yj=P(1,j)*x1+...+P(n,j)*xn.
mats=[1,0];
factors=[1,-1];
for i=2:order
    nend=size(mats,2);
    mats=[mats,mats;ones(1,nend),zeros(1,nend)]; %#ok<AGROW>
    factors=[factors,-factors]; %#ok<AGROW>
end
mats=mats(:,1:end-1);
factors=factors(1:end-1)/factorial(order);
if nargin==1 || length(groups)==order
    return
end
%% if several deviations are flagged as equal, group them together
[mats,factors]=nmfm_pol_contract(mats,factors,groups,order);
end
