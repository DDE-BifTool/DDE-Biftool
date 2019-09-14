function stability=dde_stst_eig_mxo(A,tau,varargin) 
%% compute spectrum of linear const coeff DDE
% function stability=stst_stabil(A,tau,method)
% INPUT:
%   A (n x n x ntau+1) coefficient matrices
%   tau (1 x ntau) delays
%	method method parameters 
% OUTPUT:
%	stability stability information
% COMMENT:
%       Assumes (imag(method.lms_parameter_rho)~=0)
%       This condition is tested in the code.
%
% (c) DDE-BIFTOOL v. 2.03, 05/03/2007
% Added on 05/03/2007 
%
% $Id: dde_stst_eig_mxo.m 331 2019-03-19 01:48:53Z jansieber $
%
%%
default={'method',[]};
options=dde_set_options(default,varargin,'pass_on');
if ~isempty(options.method)
    method=options.method;
else
    order=4;
    delta_region=0.1;
    [alpha,beta]=dde_stst_time_lms('mxo',4);
    default={...
        'minimal_real_part',-0.1,...
        'lms_parameter_alpha',alpha,...
        'lms_parameter_beta',beta,...
        'lms_parameter_rho',mxo_ellipse(order,delta_region),...
        'interpolation_order',order,...
        'minimal_time_step',1e-2,...
        'maximal_time_step',0.1,...
        'max_newton_iterations',6,...
        'root_accuracy',1e-6,...
        'remove_unconverged_roots',1,...
        'newheuristics_tests',2000,...
        'max_number_of_eigenvalues',100};
    method=dde_set_options(default,varargin,'pass_on');
end
if imag(method.lms_parameter_rho)==0
  error(['dde_stst_eig_mxo: imag(method.lms_parameter_rho)==0 :' ...
	 ' hence, method contains inapproriate parameters for this function. ' ...
	 'Use method=df_mthod(''stst'',1);p_stabil(stst,' ...
	 ' method.stability), to obtain an appropriate method.stability' ...
	 ' structure for this function.']);
end
%%
if length(tau)==size(A,3)
    tau=tau(2:end);
end
tau=tau(:)';
m=length(tau);
taumin=min(tau);
taumax=max(tau);
n=size(A,2);
%%
alpha=method.lms_parameter_alpha;
beta=method.lms_parameter_beta;
%k_lms = length(alpha) - 1;

interp_order=method.interpolation_order;
s_min = floor((interp_order-1)/2);
s_plus = (interp_order-1) - s_min;
real_min=method.minimal_real_part;
%% modification by JS to cope with special case of zero delay
if taumax==0 || norm(reshape(A(:,:,2:end),n,[]),'inf')==0
    A=sum(A,3);
    l0=eig(A).';
    if ~isempty(real_min)
        l0=l0(real(l0)>real_min);
    end
    if length(l0)>method.max_number_of_eigenvalues
        l0=l0(1:method.max_number_of_eigenvalues);
    end
    l1=l0.';
    stability=struct('h',NaN,'l0',l0,'l1',l1,'n1',[]);
    return
end
% end of modification (JS)
%%
if isempty(real_min) || (real_min==-Inf) 
  if taumax>0 && taumin>0.1
    real_min=-1/taumin;
  else
    real_min=-1;
  end
end

a_ellipse=real(method.lms_parameter_rho);
b_ellipse=imag(method.lms_parameter_rho);
%% change by JS: adapt nb_nu to keep effort for heuristics constant
% effort ~ nb_nu^ntau/2 * (1+n^3/100) assuming that eig is 100x faster than
% matlab
% approx number of points to create for new heuristics
if isfield(method,'newheuristics_tests')
    npts_n3=method.newheuristics_tests;
else
    npts_n3=2000;
end
nb_nu = floor((npts_n3/(1+n^3/100))^(1/length(tau))); % number of points on unit circle tested
if nb_nu<4 % too few points to approximate circle
    use_newheuristics=false;
else
    use_newheuristics=true;
end
delta_real_min = 0.1; % desired estimated accuracy for roots
%% change by js to avoid error message when heuristics gives no restriction
% (~use_newheuristics | isempty(pts))
h_upperbound=taumax*method.maximal_time_step;
h_lowerbound=taumax*method.minimal_time_step;
if use_newheuristics
    %% point cloud estimating possible spectrum, 
    % original code by K Verheyden, generalized to arbitrary number of
    % delays by JS
    eigA0 = eig(A(:,:,1)).';
    pts=get_pts_h_new(A,tau,real_min,delta_real_min,nb_nu,eigA0);
else
    pts=[];
end
%% case distinction by JS
if ~isempty(pts)
    %% heuristics gives restriction
    hh=0.9/sqrt(max((real(pts)/a_ellipse).^2+(imag(pts)/b_ellipse).^2));
else
    %% pts can be empty: use old heuristics, according to Verheyden etal
    % J. Comp. Appl. Math. 214 (2007)
    radLMS=min(a_ellipse,b_ellipse); %??
    denominator=abs(real_min)+norm(A(:,:,1),'inf');
    for j=2:m+1
        denominator=denominator+...
            norm(A(:,:,j),'inf')*exp(-abs(real_min)*tau(j-1));
    end
    hh=0.9*radLMS/denominator;
end
%%
hh=min([hh,min(tau(tau>=hh*interp_order*1e-2))/s_plus,h_upperbound]);
hh=max([hh,h_lowerbound]);
% nn=n*(k_lms + ceil(taumax/hh) + s_min);

[mu,nL]=help_stst_stabil(A,tau,hh,alpha,beta,interp_order);     %#ok<ASGLU>

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

[dummy,idx]=sort(real(lambda)); %#ok<ASGLU>
lambda=lambda(idx(end:-1:1));

if length(lambda)>method.max_number_of_eigenvalues
  lambda=lambda(1:method.max_number_of_eigenvalues);
end

stability=struct('h',hh,'l0',lambda(:),'l1',[],'n1',[]);
stability=dde_stst_lms_nwt(stability,A,tau,method);
end
%%
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
end
%%
function [S_h_row,nL,rcM]=get_S_h_row(AA,tau,hh,alpha,beta,interp_order)

% function [S_h_row,nL,rcM]=get_S_h_row(AA,tau,hh,alpha,beta,interp_order)
% INPUT:
%       AA,tau,hh,alpha,beta,interp_order
% OUTPUT:
%	S_h_row,nL,rcM
  
% (c) DDE-BIFTOOL v. 2.03, 05/03/2007
% Added on 05/03/2007   
  
n=size(AA(:,:,1),1);
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
  Gamma{j}=dde_stst_lms_nrd(epsi(j),rr(j),ss(j));
end
nL=n*(max(ll+rr)+kk-1);
S_h_row=zeros(n,nL);
% present terms:
if beta(kk)
  M=alpha(kk)*eye(n)-(hh*beta(kk))*AA(:,:,1);
  for j=1:m
    if ll(j)==ss(j)
      M=M-(hh*beta(kk)*Gamma{j}(rr(j)+ss(j)+1))*AA(:,:,j+1);
      % trick:
      ss(j)=ss(j)-1;
    end
  end
  rcM=rcond(M);
  [L_fac,U_fac,P_fac]=lu(M);
  p_fac=P_fac*(1:n)';
else
  rcM=1;
  L_fac=eye(n);
  U_fac=alpha(kk)*eye(n);
  p_fac=1:n;
end;
for j=1:kk-1
  S_h_row(:,(j-1)*n+(1:n))=-alpha(kk-j)*eye(n)+hh*beta(kk-j)*AA(:,:,1);
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
	  (hh*beta(k)*gamma_vect(rr(j)+1-o))*(U_fac\(L_fac\(AA(p_fac,:,j+1))));
    end;
  end;  
end;

end
%%
function pts=get_pts_h_new(AA,tau,real_min,delta_real_min,nb_nu,points)
%% creates cloud of points indicating area where eigenvalue can be
% function pts=get_pts_h_new(AA,tau,real_min,delta_real_min,nb_nu,eigA0)
% INPUT:
%	AA,tau,real_min,delta_real_min,nb_nu,eigA0
% OUTPUT:
%	pts
%   
% COMMENT: is implemented for arbitrary m>=1 in new version but
% computational cost increases exponentially with m m is number of delays)
%
% for formulas and motivations see 
%
% Verheyden, Luzyanina, Roose: Efficient computation of 
% characteristic roots of delay differential equations using LMS methods,
% Journal of Computational and Applied Mathematics 214 (2008) 209 – 226.
%
%
% (c) DDE-BIFTOOL v. 2.03, 05/03/2007
% Added on 05/03/2007
%
% $Id: dde_stst_eig_mxo.m 331 2019-03-19 01:48:53Z jansieber $
%
if iscell(AA)
    AA=cat(3,AA{:});
end
theta=linspace(0,2*pi,nb_nu+1);
theta=theta(1:end-1);
nu=complex(sin(theta),cos(theta));

r_min_d=real_min-delta_real_min;
ss=exp(-real_min*tau);

%% (JS) call to generalized version map_pts of hlp_get_pts_h_new (both fcns below)
%E=hlp_get_pts_h_new(AA,nu,ss); %original
E=map_pts(AA,nu,ss); 
%%
points=[points,E];

rp=real(points);
ip=imag(points);
inds_p=find(rp>=r_min_d);
rp=rp(inds_p);
ip=ip(inds_p);
%

pts=complex(rp,ip);
end


function E=map_pts(CC,nu,ss)
%% maps points on unit circle to eigenvalues
m=size(CC,3)-1;           % number of delays
nb_nu=length(nu);         % number of points on unit circle
count=BaseFromNum(nb_nu,m);
count=count(:,count(end,:)<=ceil((nb_nu+1)/2));
nc=size(count,2);
dim=size(CC(:,:,1),1);
E=NaN(dim,nc); % array of approximate Eigenvalue locations
ss=[1,ss(:)'];
Mlist=reshape(CC,dim*dim,m+1).*ss(ones(dim*dim,1),:);
nu=nu(:);
for k=1:size(count,2)
    nusel=[1;nu(count(:,k))];
    M=Mlist*nusel;
    E(:,k)=eig(reshape(M,dim,dim));
end
E=E(:).';
end
%%
function count=BaseFromNum(npoints,ndelays)
%% create counts for npoints points and ndelays delays
k=0:npoints^ndelays-1;
count=NaN(ndelays,length(k));
count(1,:)=mod(k,npoints);
for j=2:ndelays
    k=k-count(j-1,:);
    k=k/npoints;
    count(j,:)=mod(k,npoints);
end
count=count+1;
end
%%
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
    end
end
end
