function s=dde_stst_bdf_safety(alpha,beta,epsi,safe)
%% compute safety factor for given distance to imaginary axis 
% function s=dde_stst_bdf_safety(alpha,beta,epsi,safe)
% INPUT: 
%	alpha alpha-LMS parameters
%	beta beta-LMS parameters
%	epsi discretisation step in [0,2pi]
%	safe maximum distance to imaginary axis
% OUTPUT:
%	s safety factor for given distance to imaginary axis 

% (c) DDE-BIFTOOL v. 1.00, 08/04/2000
%
% $Id: dde_stst_bdf_safety.m 296 2018-09-24 21:36:56Z jansieber $
%
%%
xi=0:epsi:2*pi;

s1=0;
s2=0;
for j=1:length(alpha)
    s1=s1+alpha(j)*exp(1i*j*xi);
    s2=s2+beta(j)*exp(1i*j*xi);
end

rho=s1./s2;

s=abs(rho(floor(length(xi)/2)));

for j=1:length(xi)
    if abs(real(rho(j)))>safe
        s=min(s,abs(rho(j)));
    end
end
end
