function [ab]=dde_stst_mxo_ellipse(order,delta,ellipse_nb)
%%  
% function ab=dde_stst_mxo_ellipse(order,delta,ellipse_nb)
% INPUT:
%       order order of the (maximal order) LMS method
%       delta 0.1 OR 0.01  
%       ellipse_nb (optional, default: 1) 
% OUTPUT:
%	ab=a+bi

% (c) DDE-BIFTOOL v. 2.03, 05/03/2007
% Added on 05/03/2007 
% $Id: dde_stst_mxo_ellipse.m 366 2019-07-14 21:52:54Z jansieber $
%%
if nargin<3
  ellipse_nb=1;
end
assert(delta==0.1)
switch order
    case 4
        ab=2.193376576735286+1.679419881849218i;
    case 6
        switch ellipse_nb
            case 1
                ab=1.345424070239944+2.441717171717172i;
            case 2
                ab=2.957575757575758+2.010046051921072i;
        end
    case 8
        ab=3.522874942435445+2.204791869425589i;
end
end