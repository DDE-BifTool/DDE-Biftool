function gamma=time_nrd(epsi,r,s)

% function gamma=time_nrd(epsi,r,s)
% INPUT:
%	epsi in [0,1)
%	r number of points to the left (the past)
%	s number of points to the right (the present)
% OUTPUT:
%	gamma nordsieck interpolation vector (from past to present)

% (c) DDE-BIFTOOL v. 1.00, 15/03/2000

for j=-r:s
  g=1;
  for k=-r:s
    if k~=j,
      g=g*(epsi-k)/(j-k);
    end;
  end;
  gamma(j+r+1)=g;
end;

return;
