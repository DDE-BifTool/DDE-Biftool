function stability=dde_stst_lms_nwt(stability,AA,tau,method)
%%
% function stability=dde_stst_lms_nwt(stability,AA,tau,method)
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
% $Id: dde_stst_lms_nwt.m 366 2019-07-14 21:52:54Z jansieber $
%
%% 
epsi=method.root_accuracy;
max_n=method.max_newton_iterations;

len_l0=length(stability.l0);
n1 = zeros(size(stability.l0));
l1=n1; % modification by JS to avoid error for empty l0
res1=NaN(1,len_l0);
stability.v=NaN(size(AA,1),len_l0);
for j=1:len_l0
  lambda=stability.l0(j);
  DD=dde_stst_ch_matrix(AA,tau,lambda);
  try
    [VV,EE]=eig(DD);
  catch
    n1(j)=-1;
    l1(j)=lambda;
    continue
  end
  [dummy,kk]=min(abs(diag(EE))); %#ok<ASGLU>
  vv=VV(:,kk);
  vv=vv/norm(vv);
  stability.err(j)=norm(DD*vv,'inf');
  % normalisation cc for vv :
  cc=vv';
  while 1
      DD=dde_stst_ch_matrix(AA,tau,lambda);
      dD=dde_stst_ch_matrix(AA,tau,lambda,'deri',1);
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
          break
      end
  end
  
  if ~cvg_flag
      n1(j)=-1;
  end
  l1(j)=lambda;
  res1(j)=norm(DD*vv,'inf');
  stability.v(:,j)=vv;
end
lrem=false(size(n1));
if method.remove_unconverged_roots
    % remove n1 and sort roots
    lrem=n1==-1;
    stability.discarded=struct('l0',stability.l0(lrem),...
        'err',stability.err(lrem));
    l1=l1(~lrem);
    n1=[];
    if ~isempty(l1)
        [dummy,idx]=sort(real(l1)); %#ok<ASGLU>
        l1=l1(idx(end:-1:1)).';
    end
end
stability.l1=l1(:);
stability.err=res1(~lrem);
stability.v=stability.v(:,~lrem);
stability.n1=n1;
end