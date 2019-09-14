function [S_h_row,nL,rcM]=get_S_h_row(AA,tau,hh,alpha,beta,interp_order)

% function [S_h_row,nL,rcM]=get_S_h_row(AA,tau,hh,alpha,beta,interp_order)
% INPUT:
%       AA,tau,hh,alpha,beta,interp_order
% OUTPUT:
%	S_h_row,nL,rcM
  
% (c) DDE-BIFTOOL v. 2.03, 05/03/2007
% Added on 05/03/2007   
  
n=size(AA{1},1);
m=length(tau); 
kk=length(alpha);
%
tau_h=tau/hh; 
interp_degree=interp_order-1;
s=ceil(interp_degree/2);
r=floor(interp_degree/2);
ll=ceil(tau_h);
ss=repmat(s,1,m);  
rr=repmat(r,1,m);
ids=find(s>ll);
if ~isempty(ids)
  ss(ids)=ll(ids);
  rr(ids)=r+s-ss(ids);
end;
epsi=ll-tau_h;

no_interp_idx=find(epsi==0);
ss(no_interp_idx)=deal(0);
rr(no_interp_idx)=deal(0);
Gamma=cell(1,m);
for j=1:m
  Gamma{j}=time_nrd(epsi(j),rr(j),ss(j));
end
nL=n*(max(ll+rr)+kk-1);
S_h_row=zeros(n,nL);
% present terms:
if beta(kk)
  M=alpha(kk)*eye(n)-(hh*beta(kk))*AA{1};
  for j=1:m
    if ll(j)==ss(j)
      M=M-(hh*beta(kk)*Gamma{j}(rr(j)+ss(j)+1))*AA{j+1};
      % trick:
      ss(j)=ss(j)-1;
    end
  end
  rcM=rcond(M);
  [L_fac,U_fac,P_fac]=lu(M);
  p_fac=P_fac*[1:n]';
else
  rcM=1;
  L_fac=eye(n);
  U_fac=alpha(kk)*eye(n);
  p_fac=1:n;
end;
for j=1:kk-1
  S_h_row(:,(j-1)*n+(1:n))=-alpha(kk-j)*eye(n)+hh*beta(kk-j)*AA{1};
end;
S_h_row(:,1:((kk-1)*n))=U_fac\(L_fac\(S_h_row(p_fac,1:((kk-1)*n))));
% past terms:
for j=1:m
  gamma_vect=Gamma{j};	
  if size(S_h_row,2)<(ll(j)+rr(j)+kk-1)*n       
    error('ERROR:  size(S_h_row,2)<(ll(j)+rr(j)+kk-1)*n');
    %%% SS_HH=[SS_HH,zeros(n,(ll(j)+rr(j)+kk-1)*n-size(SS_HH,2))];
  end; 
  for k=1:kk
    for o=-ss(j):rr(j)
      S_h_row(:,(ll(j)+o+kk-k-1)*n+(1:n))=S_h_row(:,(ll(j)+o+kk-k-1)*n+(1:n))+ ...
	  (hh*beta(k)*gamma_vect(rr(j)+1-o))*(U_fac\(L_fac\(AA{j+1}(p_fac,:))));
    end;
  end;  
end;

return;



