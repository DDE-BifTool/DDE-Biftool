function stability=dde_stst_eig_bdf(A,tau,varargin)
%% compute stability information for point
% function stability=p_stabil(funcs,point,method)
% INPUT:
%   A (n x n x ntau+1) coefficient matrices
%   tau (1 x ntau) delays
%	method method parameters 
% OUTPUT:
%	stability stability information

% (c) DDE-BIFTOOL v. 2.00, 30/11/2001
% Update on 05/03/2007 ("flag_newhheur"  <=> (imag(method.lms_parameter_rho)~=0) )   
%
% $Id: dde_stst_eig_bdf.m 296 2018-09-24 21:36:56Z jansieber $
%
%%
default={'method',[]};
options=dde_set_options(default,varargin,'pass_on');
if ~isempty(options.method)
    method=options.method;
else
    order=4;
    [alpha,beta]=dde_stst_time_lms('bdf',order);
    default={...
        'minimal_real_part',-0.1,...
        'lms_parameter_alpha',alpha,...
        'lms_parameter_beta',beta,...
        'lms_parameter_rho',dde_stst_bdf_safety(alpha,beta,0.01,0.01),...
        'interpolation_order',order,...
        'minimal_time_step',1e-2,...
        'maximal_time_step',0.1,...
        'max_newton_iterations',6,...
        'root_accuracy',1e-6,...
        'remove_unconverged_roots',1,...
        'newheuristics_tests',false,...
        'max_number_of_eigenvalues',100};
    method=dde_set_options(default,varargin,'pass_on');
end
real_part=method.minimal_real_part;
a=method.lms_parameter_alpha;
b=method.lms_parameter_beta;
ord=method.interpolation_order;
h_min=method.minimal_time_step;
h_max=method.maximal_time_step;
rho=method.lms_parameter_rho;
if length(tau)==size(A,3)
    tau=tau(2:end);
end
tau=tau(:)';
[l0,h]=root_app(A,tau,real_part,a,b,ord,h_min,h_max,rho);
if length(l0)>method.max_number_of_eigenvalues
    l0=l0(1:method.max_number_of_eigenvalues);
end
stability=dde_stst_stability('h',h,'l0',l0(:),'err',NaN(1,numel(l0)));
max_n=method.max_newton_iterations;
if max_n>0 && ~isempty(l0)
    stability=dde_stst_lms_nwt(stability,A,tau,method);
end
end
%%
function [l,h]=root_app(A,tau,real_part,alpha,beta,interp_order,h_min,h_max,rho)
%% approximate rightmost characteristic roots for equilibrium
% function l=root_app(funcs,x,real_part,alpha,beta,interp_order,par,h_min,h_max,rho,d_ac)
% INPUT:
%   A (n x n x ntau+1) coefficient matrices
%   tau (1 x ntau) delays
%	real_part compute roots with real part >= real_part
%	alpha alpha-LMS parameters in R^k
%	beta beta-LMS parameters in R^k
%	interp_order order of interpolation in the past
%       par current parameter values in R^p
%	h_min minimal stepsize of discretisation h (relative to max(tau))
%	h_max maximal stepsize of discretisation h (relative to max(tau))
%	rho safety radius of LMS-method
% OUTPUT:
%	l approximations of rightmost characteristic roots
%	h stepsize used in discretisation (relative to max(tau))

% (c) DDE-BIFTOOL v. 2.00, 23/11/2001
%
% $Id: dde_stst_eig_bdf.m 296 2018-09-24 21:36:56Z jansieber $
%
%%
%% set boundary for real part depending on max(tau) if not specified
taumax=max(tau);
if isempty(real_part)
    if taumax>0
        real_part=-1/taumax;
    else
        real_part=-1000;
    end
end
%% time step for discrete time map
h=time_h(A,tau,alpha,beta,real_part,h_min,h_max,rho);
%% determine roots
n=size(A,1);
if h>0 % delay case:
    hh=h*taumax;
    S_h=root_int(A,tau,alpha,beta,hh,interp_order);
    nL=size(S_h,2);
    S_h(n+1:nL,1:nL-n)=eye(nL-n);
    s=eig(S_h);
    [dummy,I]=sort(abs(s)); %#ok<ASGLU>
    e=s(I(length(I):-1:1));
    r=exp(real_part*hh);
    l=[];
    for j=1:length(e)
        if abs(e(j))<r
            re=e(1:max(1,j-1));
            l=(log(abs(re))+1i*asin(imag(re)./abs(re)))/hh;
            break
        end
    end
else % ode case: tau(1:m)=0 or D_(1:m)=zeros(n)
    D=sum(A,3);
    le=eig(D);
    l=le(real(le)>=real_part);
end
end
%%
function S=root_int(A,tau,alpha,beta,h,interp_order)
%% first n rows of integration operator S(h) in R^(n x L)
% function S=root_int(funcs,x,alpha,beta,h,interp_order,par)
% INPUT:
%	alpha alpha-LMS parameters in R^k
%	beta beta-LMS parameters in R^k
%	h absolute stepsize > 0
%	interp_order order of interpolation in the past
%       par current parameter values in R^p
% 
% OUTPUT:
%	S first n rows of integration operator S(h) in R^(n x L)

% (c) DDE-BIFTOOL v. 2.00, 23/11/2001
%
% $Id: dde_stst_eig_bdf.m 296 2018-09-24 21:36:56Z jansieber $
%
%%
n=size(A,1);
m=length(tau);
k=length(alpha);

% present terms: A(:,:,1)

if abs(beta(k))>0
  fac=inv(alpha(k)*eye(n)-h*beta(k)*A(:,:,1));
else
  fac=eye(n)/alpha(k);
end;
S=zeros(n,n*(k-1));
for j=1:k-1
  S(1:n,(j-1)*n+(1:n))=fac*(-alpha(k-j)*eye(n)+h*beta(k-j)*A(:,:,1));
end;

% past terms:

interp_degree=interp_order-1;

s=ceil(interp_degree/2);
r=floor(interp_degree/2);

for p=1:m
    B=A(:,:,p+1);
    l=ceil(tau(p)/h);
    if s>l-2
        ss=l-2;
        rr=interp_order-ss;
    else
        ss=s;
        rr=r;
    end
    epsi=l-tau(p)/h;
    gamma_vect=dde_stst_lms_nrd(epsi,rr,ss);
    if size(S,2)<(l+rr+k-1)*n
        S=[S zeros(n,(l+rr+k-1)*n-size(S,2))];
    end
    for j=1:k
        for o=-ss:rr
            S(1:n,(l+o-j+k-1)*n+(1:n))=S(1:n,(l+o-j+k-1)*n+(1:n)) + ...
                h*beta(j)*gamma_vect(rr+1-o)*fac*B;
        end
    end
end
end
%%
function h=time_h(A,tau,alpha,beta,real_part,h_min,h_max,rho_safety)
%% recommend steplength for stability calculation (relative to max(tau))
% function h=time_h(A,tau,alpha,beta,real_part,h_min,h_max,rho_safety)
% INPUT:
%   A (n x n x ntau+1) coefficient matrices
%   tau (1 x ntau) delays
%	alpha alpha-LMS parameters in R^k
%	beta beta-LMS parameters in R^k
%	real_part to compute roots with real part >= real_part  
%	h_min minimal h (relative to max(tau))
%	h_max maximal h (relative to max(tau))
%	rho_safety safety radius of LMS-method
% OUTPUT: 
%       h recommended steplength (relative to max(tau))

% (c) DDE-BIFTOOL v. 2.00, 23/11/2001
%
% $Id: dde_stst_eig_bdf.m 296 2018-09-24 21:36:56Z jansieber $
%
%%
taumax=max(tau);
k=length(alpha);
m=length(tau);

D=A(:,:,1);
n0=norm(D);
normAB=n0+abs(real_part);
for i=1:m
    D=A(:,:,i+1);
    normAB=normAB+norm(D)*exp(-real_part*tau(i));
end;
ode=0;

% h1 for stability:

if normAB>0
    h1=0.9*rho_safety/normAB;
elseif normAB==0
    ode=1;
else % if normAB == inf
    h1=h_min*taumax/2;
end;

% h2 for invertibility:

if n0>0
    h2=0.9*alpha(k)/(beta(k)*n0);
elseif n0==0
    h2=h1;
end;

% h:

if ~ode
    if taumax>0
        h=min(min([h1 h2])/taumax,h_max);
    else
        ode=1;
    end;
end;

if ode
    disp('TIME_H warning: dde has become ode.');
    h=0;
elseif h<h_min
    h=h_min;
    if h_max>h_min
        disp('TIME_H warning: h_min is reached.');
    elseif h_max<h_min
        error('TIME_H: h_min>h_max not allowed (h_min=%g, h_max=%g).',h_min,h_max);
    end
end
end
