function [mu,nL]=help_stst_stabil(AA,tau,hh,alpha,beta,interp_order)
  
% function [mu,nL]=help_stst_stabil(AA,tau,hh,alpha,beta,interp_order)
% INPUT:
%       AA,tau,hh,alpha,beta,interp_order
% OUTPUT:
%       mu
%       nL the size of the eigenvalue problem solved

[S_h_row,nL,rcM]=get_S_h_row(AA,tau,hh,alpha,beta,interp_order);

% check condition number (post factum)
if rcM<1e-6
  warning(['Matrix ''M'' is ill', ...
	   ' conditioned: rcond = ',num2str(rcM)]);    
end

n=size(S_h_row,1);
S_h=[S_h_row; ...
     eye(nL-n), zeros(nL-n,n)];

mu=eig(S_h);

return;
