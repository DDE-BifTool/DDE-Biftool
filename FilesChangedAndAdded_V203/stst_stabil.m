function stability=stst_stabil(stst,method) 
  
% function stability=stst_stabil(stst,method)
% INPUT:
%	stst steady state point
%	method method parameters 
% OUTPUT:
%	stability stability information
% COMMENT:
%       Assumes (imag(method.lms_parameter_rho)~=0)
%       This condition is tested in the code.
  
% (c) DDE-BIFTOOL v. 2.03, 05/03/2007
% Added on 05/03/2007 

if imag(method.lms_parameter_rho)==0,
  error(['STST_STABIL: imag(method.lms_parameter_rho)==0 :' ...
	 ' hence, method contains inapproriate paramters for this function. ' ...
	 'Use method=df_mthod(''stst'',1);stst_stabil(stst,' ...
	 ' method.stability), to obtain an appropriate method.stability' ...
	 ' structure for this function.']);
end;

nb_nu = 20;
delta_real_min = 0.1;

if nargin('sys_tau')==0
  tau=stst.parameter(sys_tau);
  m=length(tau);
else
  d_ac=method.delay_accuracy;
  m=sys_ntau;
  tau=zeros(1,m);
  xx=repmat(stst.x,[1,m]);
  for j=1:m
    tau(j)=sys_tau(j,xx,stst.parameter);
  end;
  for j=1:m
    if (tau(j)<d_ac),
      s=strcat(['WARNING: delay number_', num2str(j), ...
		' is negative, no stability computed.']);
      disp(s);
      stability=[]; 
      return;
    end;
  end;
end;
taumin=min(tau);
taumax=max(tau);

xx=repmat(stst.x, 1, m+1);
AA=cell(1,m+1);
for j=0:m
  AA{j+1}=sys_deri(xx,stst.parameter,j,[],[]);
end

n=size(AA{1},2);

alpha=method.lms_parameter_alpha;
beta=method.lms_parameter_beta;
k_lms = length(alpha) - 1;

interp_order=method.interpolation_order;
order_method=interp_order;
s_min = floor((interp_order-1)/2);
s_plus = (interp_order-1) - s_min;
real_min=method.minimal_real_part;
if isempty(real_min) | (real_min==-Inf) 
  if taumax>0 & taumin>0.1
    real_min=-1/taumin;
  else
    real_min=-1;
  end
end

a_ellipse=real(method.lms_parameter_rho);
b_ellipse=imag(method.lms_parameter_rho);
  
% 
eigA0 = eig(AA{1}).';
pts=get_pts_h_new(AA,tau,real_min,delta_real_min,nb_nu,eigA0);
hh=0.9/sqrt(max((real(pts)/a_ellipse).^2+(imag(pts)/b_ellipse).^2));

hh=min([hh,min(tau(tau>=hh*interp_order*1e-2))/s_plus,taumax*method.maximal_time_step]);
hh=max([hh,taumax*method.minimal_time_step]);
nn=n*(k_lms + ceil(taumax/hh) + s_min);

[mu,nL]=help_stst_stabil(AA,tau,hh,alpha,beta,interp_order);   

% Note: nn==nL,
% except if some interpolation can be avoided ...

% Throw away mu too close to the origin
ss=exp(real_min*hh);
mu=mu(ss<=abs(mu));

% % Throw away zeros
% mu=find(abs(mu));
% % Zeros were already discarded in the previous step

lambda=mu_to_lambda(mu,hh);

% "A posteriori safeguard", 
% but here we approximate (0.9/h)*LMS^{-1}(T_delta) 
% by (0.9/h)*[ellipse(a_ellipse,b_ellipse)] ...:
lambda=lambda((real(lambda)/a_ellipse).^2+(imag(lambda)/b_ellipse).^2<=(0.9/hh)^2);

[dummy,idx]=sort(real(lambda));
lambda=lambda(idx(end:-1:1));

if length(lambda)>method.max_number_of_eigenvalues
  lambda=lambda(1:method.max_number_of_eigenvalues);
end;

stability=struct('h',hh,'l0',lambda,'l1',[],'n1',[]);

stability=stst_stabil_nwt_corr(stability,AA,tau,method);

return;