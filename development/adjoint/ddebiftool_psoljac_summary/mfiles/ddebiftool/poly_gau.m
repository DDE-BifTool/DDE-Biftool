function c=poly_gau(m)

% function c=poly_gau(m)
% INPUT:
% 	m degree polynomial
% OUTPUT:
%	c roots of m-degree gauss-legendre polynomial in [0,1]^m

% (c) DDE-BIFTOOL v. 1.00, 15/03/2000

% gauss-legendre polynomials:

if m<=0,
  err=m
  error('POLY_GAU: degree too low!');
elseif m==1,
  p=[1 0];
elseif m==2
  p=0.5*[3 0 -1];
elseif m==3,
  p=0.5*[5 0 -3 0];
elseif m==4,
  p=0.125*[35 0 -30 0 3];
elseif m==5,
  p=0.125*[63 0 -70 0 15 0];
elseif m==6,
  p=0.0625*[231 0 -315 0 105 0 -5];
elseif m==7,
  p=0.0625*[429 0 -693 0 315 0 -35 0];
elseif m==8,
  p=0.0078125*[6435 0 -12012 0 6930 0 -1260 0 35];
end;

if m>8,
  err=m
  error('POLY_GAU: degree too high!');
end;

r=sort(roots(p));

c=0.5*(r+1);

return;
