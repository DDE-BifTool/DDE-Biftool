function [sign, product] = nmfm_hopfdet(roots)
%% Monitoring function for Hopf bifurcation along steady-state branch
% function sign = detect_hopf(point)
% Purpose:
%   This routine computes the product of all sums of the eigenvalues in
%   'roots'. A sign change indicates a Hopf bifurcation!
% INPUT:
%	roots: the eigenvalues of a suspected Hopf point
% OUTPUT:
%	sign: +1 or -1
%
% $Id: nmfm_hopfdet.m 309 2018-10-28 19:02:42Z jansieber $
%
%%
product = 1;
sign = 0;
n = length(roots);

if n == 0
    fprintf('NMFM_HOPFDET: no roots to operate on!\n');
    return;
end

for i = 1:n
    for j = i+1:n
        product = product * (roots(i) + roots(j));
    end
end

if product == 0
    fprintf('NMFM_HOPFDET: test function is zero!\n');
elseif ~isfinite(product)
    fprintf('NMFM_HOPFDET: test function is infinite or not a number!\n');
else
    sign = real(product)/abs(real(product));
end
end
