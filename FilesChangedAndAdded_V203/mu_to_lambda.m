function lambda=mu_to_lambda(mu,h)
  
% lambda=mu_to_lambda(mu,h)
% INPUT:
%       mu (row vector!)
%       h 
% OUTPUT:
%	lambda (row vector)
% COMMENT:
%       Assumes that all(abs(mu)>0)
  
% (c) DDE-BIFTOOL v. 2.03, 05/03/2007
% Added on 05/03/2007 
  
len_mu=length(mu);  
lambda=zeros(1,len_mu);

for j=1:len_mu
  if real(mu(j))>=0
    lambda(j)=complex(log(abs(mu(j)))/h, ...
		      asin(imag(mu(j))/abs(mu(j)))/h);	
  elseif imag(mu(j))>=0
    lambda(j)=complex(log(abs(mu(j)))/h, ...
		      (pi-asin(imag(mu(j))/abs(mu(j))))/h);	      
  else
    lambda(j)=complex(log(abs(mu(j)))/h, ...
		      (-pi-asin(imag(mu(j))/abs(mu(j))))/h);
  end;
end;

return;
