function [mats,factors]=nmfm_pol_contract(mats,factors,groups,order)
%% remove redundant polarization combinations assuming several args are equal
%
% [mats,factors]=nmfm_pol_contract(mats,factors,groups,order)
%
% The polarization identity assumes that all arguments are different, thus,
% creating a large number of polarization directions. if several of the
% directions. This collects equal polarization directions assuming that
% some original directions were equal (which is indicated by the cell array
% group).
%
%% Inputs
%
% * mats (ndev x nterms): linear combination coefficients of directions to
% get polarization directions.
% * factors (1 x nterms): linear combination coefficients to assemble
% mixed derivative from directional derivatives
% * groups: cell array splitting up 1:ndev into ng groups indicating which
% deviations are considered equal.
% * order: order of requested derivative.
%
%% Outputs
%
% * mats (ng x npolnew) potentially smaller array of polarization
% combinations
% * factors (1 x npolnew) coefficients for linear combination
%
%% Example
%
% dde_polarization_coeffs for a general 3rd-order mixed derivative
% $d1=d^3f(x)[p,q,r]=D[p,q,r]$ returns
% >> [mats,factors]=dde_polarization_coeffs(3,{1,2,3})
% mats =
%      1     0     1     0     1     0     1
%      1     1     0     0     1     1     0
%      1     1     1     1     0     0     0
% factors =
%       0.16667 -0.16667 -0.16667  0.16667 -0.16667  0.16667 0.16667
%
% meaning that
% $d1=1/6*( D[p+q+r]^3 - D[q+r]^3 - D[p+r]^3 + D[r]^3 - D[p+q]^3 + D[q]^3 +
% D[p]^3)$
%
% However, for the 3rd-order mixed derivative $d2=d^3f(x)[p,p,q]=D[p,p,q]$
% it returns
% >> [mats,factors]=dde_polarization_coeffs(3,{[1,2],3})
% mats =
%      0     1     1     2
%      1     0     1     1
% factors =
%       0.16667           -1     -0.33333      0.16667
%
% meaning that
% $d2= D[q]^3/6 - D[p]^3 - D[p+q]^3/3 + D[2p+q]^3/6
% thus, saving computational effort compared to the general case.
%
% $Id$
%% if several deviations are flagged as equal, group them together
ng=length(groups);
nterms=size(mats,2);
smats=zeros(ng,nterms);
for i=1:ng
    smats(i,:)=sum(mats(groups{i},:),1);
end
%% divide each column of smats by its gcd
[smats,gc]=nmfm_matscale(smats);
factors=factors.*gc.^order;
%% sort smats and detect repeated entries
[mats,factors]=nmfm_matcollect(smats,factors);
end
