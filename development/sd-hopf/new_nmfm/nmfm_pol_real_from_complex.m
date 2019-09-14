function [mats,factors]=nmfm_pol_real_from_complex(order)
%% decompose D^kf(x)[x+iy]^k into linear combination of real directions
%
%% Input
%
% order (integer): power of the multi-linear form.
%% Output
%
% * mats: linear integer combinations of real and imaginary parts
% * factors: complex coefficients
%
%% Example
%
% >> [mats,factors]=nmfm_pol_real_from_complex(2)
% mats =
%      0     1     1
%      1     0     1
% factors = -1-1i   1-1i   1i
%
% because (x+1i*y)^2= (-1-1i) y^2 +(1-1i)* x^2 + 1i (x+y)^2
%
% >> [mats,factors]=nmfm_pol_real_from_complex(3)
% mats =
%      0     1     1     1     2
%      1     0     1     2     1
% factors = 3-0.5i  0.5-3i  1-1i -0.5   0.5i
%
% because (x+1i*y)^3=
%  (3-0.5i)x^3+(0.5-3i)y^3+(1-1i)(x+y)^3-0.5(x+2y)^3+0.5i(2x+y)^3
%
% Note that all nonlinear terms are real in the output.
%
% $Id$
%% 
matcell{1}=[1;0];
factorcell{1}=1;
for j=order-1:-1:1
    g={1:order-j,order-j+1:order};
    [m,f]=dde_polarization_coeffs(order,g);
    matcell{j+1}=m;
    factorcell{j+1}=f;
end
matcell{order+1}=[0;1];
factorcell{order+1}=1;
%% multiply each factor with corresponding multiple from binomial theorem
for j=0:order
    factorcell{j+1}=factorcell{j+1}*nchoosek(order,j)*1i^j;
end
%% collect equal combinations
[mats,factors]=nmfm_matcollect([matcell{:}],[factorcell{:}]);
end
    