function stability=stst_stabil_nwt_corr(stability,AA,tau,method)
%%
% function stability=stst_stabil_nwt_corr(stability,AA,tau,method)
% INPUT:
%       stability
%       AA
%       tau
%       method
% OUTPUT:
%	stability

% (c) DDE-BIFTOOL v. 2.03, 05/03/2007
% Added on 05/03/2007 
%
% $Id: stst_stabil_nwt_corr.m 366 2019-07-14 21:52:54Z jansieber $
%
%% 
epsi=method.root_accuracy;
max_n=method.max_newton_iterations;

len_l0=length(stability.l0);
n1 = zeros(size(stability.l0));
l1=n1; % modification by JS to avoid error for empty l0
for j=1:len_l0,
  lambda=stability.l0(j);
  DD=hlp_get_only_DD(AA,tau,lambda);
  try
    [VV,EE]=eig(DD);
  catch
    n1(j)=-1;
    l1(j)=lambda;
    continue;
  end;
  [dummy,kk]=min(abs(diag(EE))); %#ok<ASGLU>
  vv=VV(:,kk);
  vv=vv/norm(vv);
  % normalisation cc for vv :
  cc=vv';

  while 1, 
    [DD,dD]=hlp_get_DD_dD(AA,tau,lambda);
    JJ=[DD, dD*vv; cc, 0];
    rhs=-[DD*vv; cc*vv-1];
    dd=JJ\rhs;
    
    vv=vv+dd(1:end-1);
    lambda=lambda+dd(end);
    n1(j)=n1(j)+1;
      
    corr=abs(dd(end));
    nrm=abs(lambda);    
    cvg_flag=(corr<=max(nrm,sqrt(length(dd)-1))*epsi & ...
		abs(dd(end))<=(1+abs(lambda))*epsi*1e3);    
    if cvg_flag || n1(j)>=max_n || corr>nrm
      break;
    end;    
  end;    
  
  if ~cvg_flag,
    n1(j)=-1;
  end;
  l1(j)=lambda;
end;
  
if method.remove_unconverged_roots
  % remove n1 and sort roots
  l1=l1(n1~=-1);
  n1=[];
  if ~isempty(l1)
    [dummy,idx]=sort(real(l1)); %#ok<ASGLU>
    l1=l1(idx(end:-1:1)).';    
  end
end

stability.l1=l1;
stability.n1=n1;
  
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DD,dD]=hlp_get_DD_dD(AA,tau,lambda)
  n=size(AA{1},1);
  m=length(tau);
  DD=lambda*eye(n)-AA{1};
  dD=eye(n);
  for j=1:m
    B=AA{j+1}*exp(-lambda*tau(j));
    DD=DD-B;
    dD=dD+tau(j)*B;
  end;  
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DD=hlp_get_only_DD(AA,tau,lambda)
  n=size(AA{1},1);
  m=length(tau);
  DD=lambda*eye(n)-AA{1};
  for j=1:m
    DD=DD-AA{j+1}*exp(-lambda*tau(j));
  end;  
return;
