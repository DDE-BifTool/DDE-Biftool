function [a,b]=get_lms_ellipse_new(order,delta,ellipse_nb)
  
% function [a,b]=get_lms_ellipse_new(order,delta,ellipse_nb)
% INPUT:
%       order order of the (maximal order) LMS method
%       delta 0.1 OR 0.01  
%       ellipse_nb (optional, default: 1) 
% OUTPUT:
%	a,b

% (c) DDE-BIFTOOL v. 2.03, 05/03/2007
% Added on 05/03/2007 

if ~exist('ellipse_nb','var')  
  ellipse_nb=1;
end

switch delta
 case 0.1
  switch order
   case 4
    a=2.193376576735286e+00;
    b=1.679419881849218e+00;
   case 6
    switch ellipse_nb
     case 1
      a=1.345424070239944e+00;
      b=2.441717171717172e+00;
     case 2
      a=2.957575757575758e+00;
      b=2.010046051921072e+00;
    end     
   case 8
    a=3.522874942435445;
    b=2.204791869425589;
  end
 case 0.01
  error('GET_LMS_ELLIPSE_NEW: not yet implemented.');
  switch order
   case 4

   case 6
    switch ellipse_nb	
     case 1
      
     case 2
      
    end
   case 8
    
  end
end

return;