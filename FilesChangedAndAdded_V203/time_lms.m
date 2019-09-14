function [alpha,beta]=time_lms(method,order)

% function [alpha,beta]=time_lms(method,order)
% INPUT:
%	method some properly abbreviated lms-method
%       order step order of the method
% OUTPUT:
%       alpha alpha-LMS parameters ordered past to present
%       beta beta-LMS parameters ordered past to present
% COMMENT:
%	method can be 'bdf' (backwards differentiation methods), 
%	'adb' (Adams-Bashforth methods), 'adm' (Adams-Moulton methods)

% (c) DDE-BIFTOOL v. 1.00, 15/03/2000
% Update on 05/03/2007 ("mxo" ...)   

if method=='bdf'
  if order==1 % aka backward euler
    alpha=[-1 1];
    beta=[0 1];
  elseif order==2
    alpha=[1/3 -4/3 1];
    beta=[0 0 2/3];
  elseif order==3
    alpha=[-2/11 9/11 -18/11 1];
    beta=[0 0 0 6/11];
  elseif order==4
    alpha=[3/25 -16/25 36/25 -48/25 1];
    beta=[0 0 0 0 12/25];
  elseif order==5
    alpha=[-12/137 75/137 -200/137 300/137 -300/137 1];
    beta=[0 0 0 0 0 60/137];
  elseif order==6  
    alpha=[10/147 -72/147 225/147 -400/147 450/147 -360/147 1];
    beta=[0 0 0 0 0 0 60/147];
  else
    err=order
    error('TIME_LMS: requested bdf order not supported.');
  end;
elseif method=='adb'
  if order==1 % aka forward euler
    alpha=[-1 1];
    beta=[1 0];
  elseif order==2
    alpha=[0 -1 1];
    beta=[-1 3 0]/2;
  elseif order==3
    alpha=[0 0 -1 1];
    beta=[5 -16 23 0]/12;
  elseif order==4
    alpha=[0 0 0 -1 1];
    beta=[-9 37 -59 55 0]/24;
  elseif order==5
    alpha=[0 0 0 0 -1 1];
    beta=[251 -1274 2616 -2774 1901 0]/720;
  elseif order==6
    alpha=[0 0 0 0 0 -1 1];
    beta=[-475 2877 -7298 9982 -7923 4277 0]/1440;
  else
    err=order
    error('TIME_LMS: requested adb order not supported.');
  end;
elseif method=='adm'
  if order==1 % aka trapezium rule
    alpha=[-1 1];
    beta=[1 1]/2;
  elseif order==2
    alpha=[0 -1 1];
    beta=[-1 8 5]/12;
  elseif order==3
    alpha=[0 0 -1 1];
    beta=[1 -5 19 9]/24;
  elseif order==4
    alpha=[0 0 0 -1 1];
    beta=[-19 106 -264 646 251]/720;
  elseif order==5
    alpha=[0 0 0 0 -1 1];
    beta=[27 -173 482 -798 1427 475]/1440;
  else
    err=order
    error('TIME_LMS: requested adm order not supported.');
  end;
elseif method=='mxo'
  % maximal order LMS method
 switch order
  case 4
   % aka Milne-Simpson, order 4
   alpha=[-1 0 1];
   beta=[1 4 1]/3;
  case 6
   alpha=[-1 -27/11 27/11 1];
   beta=[3 27 27 3]/11;
  case 8
   alpha=[-1 -32/5 0 32/5 1];
   beta=[6 96 216 96 6]/25;
  otherwise
   err=order
   error('TIME_LMS: requested mxo order not supported.');
 end;
else
  err=method
  error('TIME_LMS: requested method not recognized.');
end;

return;
