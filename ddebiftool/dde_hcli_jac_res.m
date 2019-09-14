function [J,res]=dde_hcli_jac_res(funcs,pt,free_par,method,varargin)
%% Jacobian and residual for connecting orbits
% function [J,res]=dde_hcli_jac_res(funcs,pt,free_par,method,pref,varargin)
% INPUT:
%   funcs problem functions
%	method: method parameters
%   point components
%	T period 
%	profile profile in R^(n x m*l+1)
%	t representation points in [0,1]^(m*l+1)
%	deg degree piecewise polynomial
%	par current parameter values in R^p
%	free_par free parameters numbers in N^d 
%	ph use phase condition or not (1 or 0)
%       lambda_v unstable eigenvalues of x1 (in R^s1)
%       lambda_w unstable eigenvalues of x2 (in R^s2)
%       v unstable eigenvectors of x1 (in R^n x s1)
%       w unstable eigenvectors of x2 (in R^n x s2)
%       alpha coefficients of initial functionsegment (in R^s1)
%       epsilon global coefficient of initial function segment (in R)
%       x1 steady state solution at t=-inf (in R^n)
%       x2 steady state solution at t=+inf (in R^n)
%       (optional) previous: previous solution, used in phase condition
% OUTPUT: 
%	J   jacobian in
%            R^(n*m*l+3*n+(s1+s2)*n+s1+2*s2+1 x n*m*l+3*n+(s1+s2)*n+2*s1+s2+d)
%	res residual in R^(n*m*l+3*n+(s1+s2)*n+s1+2*s2+1)
%
% (c) DDE-BIFTOOL v. 2.02, 16/6/2002 
%
% $Id: dde_hcli_jac_res.m 305 2018-10-05 20:31:58Z jansieber $
%
%#ok<*AGROW>
%% inputs/components
default={'pref',pt,'collocation_parameters',method.collocation_parameters,...
    'phase_condition',method.phase_condition};
options=dde_set_options(default,varargin,'pass_on');
previous=options.pref;
if isempty(options.pref)
    previous=pt;
end
c=options.collocation_parameters;
T=pt.period;
profile=pt.profile;
t=pt.mesh;
deg=pt.degree;
par=pt.parameter;
ph=options.phase_condition;
lambda_v=pt.lambda_v;
lambda_w=pt.lambda_w;
v=pt.v;
w=pt.w;
alpha=pt.alpha;
epsilon=pt.epsilon;
x1=pt.x1;
x2=pt.x2;
%% functions
sys_tau=funcs.sys_tau;
sys_rhs=funcs.sys_rhs;
sys_deri=funcs.sys_deri;
%% if previous point not present
if isempty(previous)
    previous=pt;
else
    previous=p_remesh(previous,pt.degree,pt.mesh);
end
m=deg;
n=size(profile,1); % system dimension

tp_del=funcs.tp_del;
if tp_del==0,
    n_tau=sys_tau(); % delay numbers
    tau=par(n_tau); % delay values
    nb_tau=length(n_tau); % number of delays
else
    error('HCLI_JAC: computing connected orbits is not implemented for equations with state-dependent delays');
end;

l=(length(t)-1)/m; % number of intervals
s1=length(lambda_v); % number of unstable eigenvalues of x1
s2=length(lambda_w); % idem x2
nb_par=length(free_par); % number of free parameters

% check:

if l~=floor(l)
    err=[length(t) m];
    error('HCLI_JAC: t (length=%d) does not contain l intervals of m=%d points!',...
        err(1),err(2));
end;

if length(c)~=m && ~isempty(c)
    err=[length(c) m];
    error('HCLI_JAC: wrong number (=%d) of collocation parameters (degree=%d)!',...
        err(1),err(2));
end;

if T<0,
    err=T;
    error('HCLI_JAC: period T=%g, became smaller than zero!',err);
end;

if max(tau)>T,
    err=[max(tau) T];
    error('HCLI_JAC: period %g became smaller than maximal delay %g!',err(2),err(1));
end;

% init J, res:

nml=n*m*l;
ml=m*l;
nml_n_1=nml+n+1;

J=zeros(nml+3*n+(s1+s2)*n+s1+2*s2+1+1,nml_n_1+2*n+(s1+s2)*n+2*s1+s2+nb_par);
res=zeros(nml+3*n+(s1+s2)*n+s1+2*s2+1+1,1);

% initialize gauss points

if isempty(c)
    c=poly_gau(m);
    gauss_c=c;
    non_gauss=0;
else
    gauss_c=poly_gau(m);
    non_gauss=1;
end;

% determination of gaussian weights:

gauss_abs=ones(1,m);
g=poly_gau(m-1);
for k=1:m
    for j=1:m-1
        gauss_abs(k)=gauss_abs(k)/(gauss_c(k)-g(j));
    end;
    for j=1:m
        if j~=k
            gauss_abs(k)=gauss_abs(k)/(gauss_c(k)-gauss_c(j));
        end;
    end;
end;
gauss_abs=gauss_abs/sum(gauss_abs);

% more initialisation:

all_dPa=poly_dla((0:m)/m,c(:)')';
all_Pa=poly_lgr((0:m)/m,c(:)')';


% for all collocation points, make equation:

tTT=tau/(T*T);

for l_i=1:l
    % determine index for a in profile:
    index_a=(l_i-1)*m+1;
    i_index_a=(index_a-1)*n;
    t_start=t(index_a);
    hhh=t(index_a+m)-t_start;
    for m_i=1:m
        i=index_a+m_i-1;
        i_range=(i-1)*n+1:i*n;
        
        % determine c
        
        col=t_start+c(m_i)*hhh;
        
        % determine dPa for a:
        
        dPa=all_dPa(m_i,:)/hhh;
        
        % add dPa for da in J:
        
        for k=0:m,
            kk=k*n;
            J(i_range,i_index_a+kk+1:i_index_a+kk+n)= ...
                J(i_range,i_index_a+kk+1:i_index_a+kk+n)+eye(n)*dPa(k+1);
        end;
        
        % add sum a*dPa to res:
        
        u_prime=profile(:,index_a:index_a+m)*dPa';
        u_prime_previous=previous.profile(:,index_a:index_a+m)*dPa';
        
        res(i_range,1)=u_prime;
        
        % determine Pa for a:
        
        Pa=all_Pa(m_i,:);
        
        % phase_condition:
        
        if ph && ~non_gauss
            fup=gauss_abs(m_i)*hhh*u_prime_previous';
            i_l_i=(l_i-1)*m*n;
            for q=0:m
                qq=q*n;
                J(nml+2*n+s1*(n+1)+s2*(n+2)+n+1,i_l_i+1+qq:i_l_i+qq+n)= ...
                    J(nml+2*n+s1*(n+1)+s2*(n+2)+n+1,i_l_i+1+qq:i_l_i+qq+n) + Pa(q+1)*fup;
            end;
        end;
        
        % compute x:
        
        x=profile(:,index_a:index_a+m)*Pa';
        
        for tau_i=1:nb_tau
            
            % determine c_tau:
            c_tau_i=col-tau(tau_i)/T;
            c_tau(tau_i)=c_tau_i;
            if c_tau_i<0
                initial_function_segment(tau_i)=1;
            else
                initial_function_segment(tau_i)=0;
            end;
            
            if (~initial_function_segment(tau_i))
                
                % determine index for b in profile:
                
                index_b_i=length(t)-m;
                while (c_tau_i<t(index_b_i))
                    index_b_i=index_b_i-m;
                end;
                
                index_b(tau_i)=index_b_i;
                
                % c transformed to [0,1] and hhh_tau
                
                hhh_tau(tau_i)=t(index_b_i+m)-t(index_b_i);
                c_tau_trans(tau_i)=(c_tau_i-t(index_b_i))/hhh_tau(tau_i);
                
                % determine Pb for b:
                
                Pb(tau_i,:)=poly_elg(m,c_tau_trans(tau_i));
                
                % compute x_tau:
                
                x_tau(:,tau_i)=profile(:,index_b_i:index_b_i+m)*Pb(tau_i,:)';
            else
                x_tau(:,tau_i)=x1+epsilon*v*(alpha.*exp(T*lambda_v*c_tau_i));
            end;
            
        end; % for tau_i=1:d
        
        % determine A0:
        
        T_A0=T*sys_deri([x x_tau],par,0,[],[]);
        
        % add -T*A0*Pa for da in J:
        
        for k=0:m
            kk=k*n;
            J(i_range,i_index_a+kk+1:i_index_a+kk+n)= ...
                J(i_range,i_index_a+kk+1:i_index_a+kk+n) ...
                -T_A0*Pa(k+1);
        end;
        
        % determine parameter derivatives:
        
        for p_i=1:nb_par
            df=sys_deri([x x_tau],par,[],free_par(p_i),[]);
            J(i_range,nml_n_1+p_i)=-T*df;
        end;
        
        % compute f(x,x_tau):
        
        f=sys_rhs([x x_tau],par);
        
        % add -f for dT in J:
        
        J(i_range,nml_n_1)=J(i_range,nml_n_1)-f;
        
        % add Tf in res:
        
        res(i_range,1)=res(i_range,1)-T*f;
        
        for t_i=1:nb_tau
            
            % determine A1:
            
            T_A1=T*sys_deri([x x_tau],par,t_i,[],[]);
            
            if (~initial_function_segment(t_i))
                
                index_b(t_i);
                % add -T*A1*Pb for db in J:
                
                i_index_b=(index_b(t_i)-1)*n;
                for k=0:m
                    kk=n*k;
                    J(i_range,i_index_b+kk+1:i_index_b+kk+n)= ...
                        J(i_range,i_index_b+kk+1:i_index_b+kk+n) - T_A1*Pb(t_i,k+1);
                end;
                
                % determine dPb for b:
                
                dPb=poly_del(m,c_tau_trans(t_i))/hhh_tau(t_i);
                
                % add -T*A1*sum b*dP*dc_tau for dT in J:
                
                J(i_range,nml_n_1)=J(i_range,nml_n_1)+ ...
                    -T_A1*(profile(:,index_b(t_i):index_b(t_i)+m)*dPb')*tTT(t_i);
                
                % special case where delay is the free parameter
                for p_i=1:length(free_par)
                    if free_par(p_i)==n_tau(t_i)
                        J(i_range,nml_n_1+p_i)=J(i_range,nml_n_1+p_i) ...
                            +T_A1/T*(profile(:,index_b(t_i):index_b(t_i)+m)*dPb');
                    end;
                end;
                
            else
                % derivatives to alpha_k
                
                for k=1:s1
                    J(i_range,nml_n_1+nb_par+2*n+s1*(n+1)+s2*(n+1)+k)=...
                        J(i_range,nml_n_1+nb_par+2*n+s1*(n+1)+s2*(n+1)+k)-...
                        T_A1*epsilon*v(:,k)*exp(T*lambda_v(k)*c_tau(t_i));
                end;
                
                % derivatives to x1
                
                J(i_range,nml_n_1+nb_par+1:nml_n_1+nb_par+n)=...
                    J(i_range,nml_n_1+nb_par+1:nml_n_1+nb_par+n)-T_A1;
                
                % derivatives to v(k)
                
                for k=1:s1
                    J(i_range,nml_n_1+nb_par+2*n+(k-1)*n+1:nml_n_1+nb_par+2*n+k*n)=...
                        J(i_range,nml_n_1+nb_par+2*n+(k-1)*n+1:nml_n_1+nb_par+2*n+k*n)-...
                        T_A1*alpha(k)*epsilon*exp(T*lambda_v(k)*c_tau(t_i));
                end;
                
                % derivatives to lambda_v(k)
                
                for k=1:s1
                    J(i_range,nml_n_1+nb_par+2*n+s1*n+k)=...
                        J(i_range,nml_n_1+nb_par+2*n+s1*n+k)...
                        -T_A1*alpha(k)*v(:,k)*epsilon*...
                        T*c_tau(t_i)*exp(T*lambda_v(k)*c_tau(t_i));
                end;
                
                % derivative to T
                for k=1:s1
                    J(i_range,nml_n_1)=J(i_range,nml_n_1)-...
                        T_A1*alpha(k)*v(:,k)*exp(T*lambda_v(k)*c_tau(t_i))*...
                        lambda_v(k)*epsilon*(c_tau(t_i)+tau(t_i)/T);
                end;
                
                % special case where delay is a free parameter
                for p_i=1:length(free_par)
                    if free_par(p_i)==n_tau(t_i)
                        for k=1:s1
                            J(i_range,nml_n_1+p_i)=J(i_range,nml_n_1+p_i)...
                                +T_A1*epsilon*...
                                v(:,k)*alpha(k)*exp(T*lambda_v(k)*c_tau(t_i))*lambda_v(k);
                        end;
                    end;
                end;
                
            end; % if (index_b(t_i))
            
        end; % for t_i=1:nb_tau
        
    end;  % for m_i
    
end;  % for l_i

% end of collocation equation section

% steady state point condition x1

xx1=x1(:,ones(nb_tau+1,1));

res(nml+1:nml+n,1)=sys_rhs(xx1,par);

% derivatives to x(t-tau(i))

for i=0:nb_tau
    J(nml+1:nml+n,nml_n_1+nb_par+1:nml_n_1+nb_par+n)=...
        J(nml+1:nml+n,nml_n_1+nb_par+1:nml_n_1+nb_par+n)...
        +sys_deri(xx1,par,i,[],[]);
end;

% derivatives to eta

for j=1:nb_par
    J(nml+1:nml+n,nml_n_1+j)=...
        J(nml+1:nml+n,nml_n_1+j)+sys_deri(xx1,par,[],free_par(j),[]);
end;


% steady state point condition x2

% derivatives to x(t-tau(i))

xx2=x2(:,ones(nb_tau+1,1));

res(nml_n_1:nml+2*n)=sys_rhs(xx2,par);

for i=0:nb_tau
    J(nml_n_1:nml+2*n,nml_n_1+nb_par+n+1:nml_n_1+nb_par+2*n)=...
        J(nml_n_1:nml+2*n,nml_n_1+nb_par+n+1:nml_n_1+nb_par+2*n)...
        +sys_deri(xx2,par,i,[],[]);
end;

% derivatives to eta

for j=1:nb_par
    J(nml_n_1:nml+2*n,nml_n_1+j)=...
        J(nml_n_1:nml+2*n,nml_n_1+j)+sys_deri(xx2,par,[],free_par(j),[]);
end;


% eigenvalues and eigenvector equations

for k=1:s1
    res(nml+2*n+(k-1)*n+1:nml+2*n+k*n,1)=...
        eye(n)*lambda_v(k)*v(:,k)-sys_deri(xx1,par,0,[],[])*v(:,k);
    for i=1:nb_tau
        res(nml+2*n+(k-1)*n+1:nml+2*n+k*n,1)=...
            res(nml+2*n+(k-1)*n+1:nml+2*n+k*n,1)...
            -sys_deri(xx1,par,i,[],[])*exp(-lambda_v(k)*tau(i))*v(:,k);
    end;
end;

% entries for the v_k variables

for k=1:s1
    J(nml+2*n+(k-1)*n+1:nml+(k+2)*n,...
        nml_n_1+nb_par+2*n+(k-1)*n+1:nml_n_1+nb_par+2*n+k*n)=...
        J(nml+2*n+(k-1)*n+1:nml+2*n+k*n,...
        nml_n_1+nb_par+2*n+(k-1)*n+1:nml_n_1+nb_par+2*n+k*n)+...
        eye(n)*lambda_v(k)-sys_deri(xx1,par,0,[],[]);
    for i=1:nb_tau
        J(nml+2*n+1+(k-1)*n:nml+2*n+k*n,...
            nml_n_1+nb_par+2*n+(k-1)*n+1:nml_n_1+nb_par+2*n+k*n)=...
            J(nml+2*n+1+(k-1)*n:nml+2*n+k*n,...
            nml_n_1+nb_par+2*n+(k-1)*n+1:nml_n_1+nb_par+2*n+k*n)-...
            sys_deri(xx1,par,i,[],[])*exp(-lambda_v(k)*tau(i));
    end;
end;

% entry for the lambda_v_k variables

for k=1:s1
    J(nml+2*n+(k-1)*n+1:nml+(k+2)*n,nml_n_1+nb_par+2*n+s1*n+k)=...
        J(nml+2*n+(k-1)*n+1:nml+(k+2)*n,nml_n_1+nb_par+2*n+s1*n+k)+v(:,k);
    for i=1:nb_tau
        J(nml+2*n+(k-1)*n+1:nml+(k+2)*n,nml_n_1+nb_par+2*n+s1*n+k)=...
            J(nml+2*n+(k-1)*n+1:nml+(k+2)*n,nml_n_1+nb_par+2*n+s1*n+k)+...
            sys_deri(xx1,par,i,[],[])*v(:,k)*tau(i)*exp(-lambda_v(k)*tau(i));
    end;
end;

% entries for the free parameters

for k=1:s1
    for j=1:nb_par
        J(nml+2*n+(k-1)*n+1:nml+(k+2)*n,nml_n_1+j)=...
            J(nml+2*n+(k-1)*n+1:nml+(k+2)*n,nml_n_1+j)-...
            sys_deri(xx1,par,0,free_par(j),[])*v(:,k);
        for i=1:nb_tau
            J(nml+2*n+(k-1)*n+1:nml+(k+2)*n,nml_n_1+j)=...
                J(nml+2*n+(k-1)*n+1:nml+(k+2)*n,nml_n_1+j)...
                -sys_deri(xx1,par,i,free_par(j),[])*v(:,k)*exp(-lambda_v(k)*tau(i));
            % special case if delay is free parameter
            if free_par(j)==n_tau(i)
                J(nml+2*n+(k-1)*n+1:nml+(k+2)*n,nml_n_1+j)=...
                    J(nml+2*n+(k-1)*n+1:nml+(k+2)*n,nml_n_1+j)...
                    +lambda_v(k)*sys_deri(xx1,par,i,[],[])*v(:,k)*...
                    exp(-lambda_v(k)*tau(i));
            end;
        end;
    end;
end;

% entry for the steady state point

for k=1:s1
    for i=0:nb_tau
        J(nml+2*n+(k-1)*n+1:nml+(2+k)*n,nml_n_1+nb_par+1:nml_n_1+nb_par+n)=...
            J(nml+2*n+(k-1)*n+1:nml+(2+k)*n,nml_n_1+nb_par+1:nml_n_1+nb_par+n)-...
            sys_deri(xx1,par,[0 i],[],v(:,k));
        for j=1:nb_tau
            J(nml+2*n+(k-1)*n+1:nml+(2+k)*n,nml_n_1+nb_par+1:nml_n_1+nb_par+n)=...
                J(nml+2*n+(k-1)*n+1:nml+(2+k)*n,nml_n_1+nb_par+1:nml_n_1+nb_par+n)-...
                sys_deri(xx1,par,[j i],[],v(:,k))*exp(-lambda_v(k)*tau(j));
        end;
    end;
end;

% entries for the normalization condition: sum(vk_i^2=1);

for k=1:s1
    res(nml+2*n+s1*n+k,1)=v(:,k)'*v(:,k)-1;
    J(nml+2*n+s1*n+k,nml_n_1+nb_par+2*n+(k-1)*n+1:nml_n_1+nb_par+2*n+k*n)=...
        2*v(:,k)';
end;

% characteristic equations for w_k

for k=1:s2
    res(nml+2*n+s1*(n+1)+(k-1)*n+1:nml+2*n+s1*(n+1)+k*n,1)=...
        lambda_w(k)'*w(:,k)-sys_deri(xx2,par,0,[],[])'*w(:,k);
    for i=1:nb_tau
        res(nml+2*n+s1*(n+1)+(k-1)*n+1:nml+2*n+s1*(n+1)+k*n,1)=...
            res(nml+2*n+s1*(n+1)+(k-1)*n+1:nml+2*n+s1*(n+1)+k*n,1)-...
            (sys_deri(xx2,par,i,[],[])*exp(-lambda_w(k)*tau(i)))'*w(:,k);
    end;
end;

% entries for the w_k variables

for k=1:s2
    J(nml+2*n+s1*(n+1)+(k-1)*n+1:nml+2*n+s1*(n+1)+k*n,...
        nml_n_1+nb_par+2*n+s1*(n+1)+(k-1)*n+1:nml_n_1+nb_par+2*n+s1*(n+1)+k*n)=...
        J(nml+2*n+s1*(n+1)+(k-1)*n+1:nml+2*n+s1*(n+1)+k*n,...
        nml_n_1+nb_par+2*n+s1*(n+1)+(k-1)*n+1:...
        nml_n_1+nb_par+2*n+s1*(n+1)+k*n)+...
        eye(n)*lambda_w(k)'-(sys_deri(xx2,par,0,[],[]))';
    
    for j=1:nb_tau
        J(nml+2*n+s1*(n+1)+(k-1)*n+1:nml+2*n+s1*(n+1)+k*n,...
            nml_n_1+nb_par+2*n+s1*(n+1)+(k-1)*n+1:...
            nml_n_1+nb_par+2*n+s1*(n+1)+k*n)=...
            J(nml+2*n+s1*(n+1)+(k-1)*n+1:nml+2*n+s1*(n+1)+k*n,...
            nml_n_1+nb_par+2*n+s1*(n+1)+(k-1)*n+1:...
            nml_n_1+nb_par+2*n+s1*(n+1)+k*n)-...
            (sys_deri(xx2,par,j,[],[])*exp(-lambda_w(k)*tau(j)))';
    end;
end;

% entry for the lambda_w_k variables

for k=1:s2
    J(nml+2*n+s1*(n+1)+(k-1)*n+1:...
        nml+2*n+s1*(n+1)+k*n,nml_n_1+nb_par+2*n+s1*(n+1)+s2*n+k)=...
        J(nml+2*n+s1*(n+1)+(k-1)*n+1:...
        nml+(k+2)*n+s1*(n+1),nml_n_1+nb_par+2*n+s1*(n+1)+s2*n+k)+w(:,k);
    for j=1:nb_tau
        J(nml+2*n+s1*(n+1)+(k-1)*n+1:...
            nml+(k+2)*n+s1*(n+1),nml_n_1+nb_par+2*n+s1*(n+1)+s2*n+k)=...
            J(nml+2*n+s1*(n+1)+(k-1)*n+1:...
            nml+(k+2)*n+s1*(n+1),nml_n_1+nb_par+2*n+s1*(n+1)+s2*n+k)+...
            (sys_deri(xx2,par,j,[],[])*exp(-lambda_w(k)*tau(j)))'*w(:,k)*tau(j);
    end;
end;
% entries for the free parameters

for k=1:s2
    for j=1:nb_par
        J(nml+2*n+s1*(n+1)+(k-1)*n+1:nml+(2+s1+k)*n+s1,nml_n_1+j)=...
            J(nml+2*n+s1*(n+1)+(k-1)*n+1:nml+(2+s1+k)*n+s1,nml_n_1+j)-...
            (sys_deri(xx2,par,0,free_par(j),[]))'*w(:,k);
        for i=1:nb_tau
            J(nml+2*n+s1*(n+1)+(k-1)*n+1:nml+s1*(n+1)+(k+2)*n,nml_n_1+j)=...
                J(nml+2*n+s1*(n+1)+(k-1)*n+1:nml+(k+2)*n+s1*(n+1),nml_n_1+j)...
                -(sys_deri(xx2,par,i,free_par(j),[])*exp(-lambda_w(k)*tau(i)))'*w(:,k);
            % special case if delay is free parameter
            if free_par(j)==n_tau(i)
                J(nml+2*n+s1*(n+1)+(k-1)*n+1:nml+2*n+s1*(n+1)+k*n,nml_n_1+j)=...
                    J(nml+2*n+s1*(n+1)+(k-1)*n+1:nml+2*n+s1*(n+1)+k*n,nml_n_1+j)...
                    +lambda_w(k)*(sys_deri(xx2,par,i,[],[])*exp(-lambda_w(k)*tau(i)))'*w(:,k);
            end;
        end;
    end;
end;

% entry for the steady state point
In=eye(n);
for k=1:s2
    for i=0:nb_tau
        for j=1:n
            J(nml+2*n+(k-1)*n+s1*(n+1)+j,nml_n_1+nb_par+n+(1:n))=...
                J(nml+2*n+(k-1)*n+s1*(n+1)+j,nml_n_1+nb_par+n+(1:n))-...
                w(:,k)'*sys_deri(xx2,par,[i 0],[],In(:,j));
        end;
        for j=1:nb_tau
            for l=1:n
                J(nml+2*n+(k-1)*n+s1*(n+1)+l,nml_n_1+nb_par+n+(1:n))=...
                    J(nml+2*n+(k-1)*n+s1*(n+1)+l,nml_n_1+nb_par+n+(1:n))-...
                    w(:,k)'*sys_deri(xx2,par,[i j],[],In(:,l))*...
                    exp(-lambda_w(k)*tau(j));
            end ;
        end;
    end;
end;

% entries for the normalization condition: sum(wk_i^2=1);

for k=1:s2
    res(nml+2*n+s1*(n+1)+s2*n+k,1)=w(:,k)'*w(:,k)-1;
    J(nml+2*n+s1*n+s1+s2*n+k,...
        nml_n_1+nb_par+2*n+s1*(n+1)+(k-1)*n+1:nml_n_1+nb_par+2*n+s1*(n+1)+k*n)=...
        2*w(:,k)';
end;

% integral orthogonality conditions

for k=1:s2
    res(nml+2*n+s1*n+s1+s2*n+s2+k,1)=w(:,k)'*(profile(:,ml+1)-x2);
end;

% first term of derivative to w_k
for k=1:s2
    J(nml+2*n+s1*n+s1+s2*n+s2+k,...
        nml_n_1+nb_par+2*n+s1*n+s1+(k-1)*n+1:nml_n_1+nb_par+2*n+s1*n+s1+k*n)=...
        J(nml+2*n+s1*n+s1+s2*n+s2+k,...
        nml_n_1+nb_par+2*n+s1*n+s1+(k-1)*n+1:...
        nml_n_1+nb_par+2*n+s1*n+s1+k*n)+...
        (profile(:,ml+1)-x2)';
end;

% derivative to y(t_l)

for k=1:s2
    J(nml+2*n+s1*n+s1+s2*n+s2+k,nml+1:nml+n)=w(:,k)';
end;

for l=1:nb_tau
    if (1-tau(l)/T)~=1  % complicated here because of machine precision
        % problems with tau=0 and variable
        A=sys_deri(xx2,par,l,[],[]);
        % recall that gauss_c=poly_gau(m);
        
        tau_T=1-tau(l)/T;
        index=length(t)-m;
        %save tauT tau_T t T tau;
        while tau_T<t(index)
            index=index-m;
        end;
        K=index;
        gauss_c_trans(1:m)=gauss_c*(t(index+m)-tau_T)+tau_T;
        for i=1:m
            Pa(index-K+i,:)=...
                poly_elg(m,(gauss_c_trans(index-K+i)-tau_T)/(t(index+m)-tau_T));
            gauss_y(1:n,index-K+i)=profile(:,index:index+m)*Pa(index-K+i,:)'-x2;
            gauss_weight(index-K+i)=gauss_abs(i)*T*(t(index+m)-tau_T);
        end;
        
        for index=K+m:m:ml-m+1
            % transform gauss_c to points in [t(index),t(index+m)]
            gauss_c_trans(index-K+1:index-K+m)=...
                gauss_c*(t(index+m)-t(index))+t(index);
            % evaluate the profile in those transformed gauss points
            for i=1:m
                Pa(index-K+i,:)=poly_elg(m,(gauss_c_trans(index-K+i)-...
                    t(index))/(t(index+m)-t(index)));
                gauss_y(1:n,index-K+i)=profile(:,index:index+m)*Pa(index-K+i,:)'-x2;
                gauss_weight(index-K+i)=gauss_abs(i)*T*(t(index+m)-t(index));
            end;
        end;
        theta=(gauss_c_trans*T-T+tau(l));
        
        for k=1:s2
            for i=1:length(gauss_weight)
                res(nml+2*n+s1*(n+1)+s2*n+s2+k,1)=...
                    res(nml+2*n+s1*(n+1)+s2*n+s2+k,1)+...
                    gauss_weight(i)*w(:,k)'*exp(-lambda_w(k)*theta(i))* ...
                    A*gauss_y(:,i);
            end;
        end;
        
        for k=1:s2
            % derivative to w_k
            for i=1:length(gauss_weight)
                J(nml+2*n+s1*n+s1+s2*n+s2+k,...
                    nml_n_1+nb_par+2*n+s1*n+s1+(k-1)*n+1:...
                    nml_n_1+nb_par+2*n+s1*n+s1+k*n)=...
                    J(nml+2*n+s1*n+s1+s2*n+s2+k,...
                    nml_n_1+nb_par+2*n+s1*n+s1+(k-1)*n+1:...
                    nml_n_1+nb_par+2*n+s1*n+s1+k*n)+...
                    gauss_weight(i)*exp(-lambda_w(k)*theta(i))*(A*gauss_y(:,i))';
            end;
        end;
        
        % derivative to lambda_w_k
        
        for k=1:s2
            for i=1:length(gauss_weight)
                J(nml+2*n+s1*n+s1+s2*n+s2+k,nml_n_1+nb_par+2*n+s1*n+s1+s2*n+k)=...
                    J(nml+2*n+s1*n+s1+s2*n+s2+k,nml_n_1+nb_par+2*n+s1*n+s1+s2*n+k)-...
                    gauss_weight(i)*w(:,k)'*theta(i)*...
                    exp(-lambda_w(k)*theta(i))*A*gauss_y(:,i);
            end;
        end;
        
        % derivative to x2
        
        for k=1:s2
            J(nml+2*n+s1*n+s1+s2*n+s2+k,nml_n_1+nb_par+n+1:...
                nml_n_1+nb_par+2*n)=-w(:,k)';
            for i=1:length(gauss_weight)
                dA1_dx2= sys_deri(xx2,par,[l 0],[],gauss_y(:,i));
                for j=1:nb_tau
                    dA1_dx2= dA1_dx2+sys_deri(xx2,par,[l j],[],gauss_y(:,i));
                end;
                J(nml+2*n+s1*n+s1+s2*n+s2+k,nml_n_1+nb_par+n+1:nml_n_1+nb_par+2*n)=...
                    J(nml+2*n+s1*n+s1+s2*n+s2+k,nml_n_1+nb_par+n+1:nml_n_1+nb_par+2*n)+...
                    gauss_weight(i)*w(:,k)'*exp(-lambda_w(k)*theta(i))*(dA1_dx2-A);
            end;
        end;
        
        % derivative to free parameters
        
        for j=1:nb_par
            dA1_deta = sys_deri(xx2,par,l,free_par(j),[]);
            for k=1:s2
                for i=1:length(gauss_weight)
                    J(nml+2*n+s1*n+s1+s2*n+s2+k,nml_n_1+j)=...
                        J(nml+2*n+s1*n+s1+s2*n+s2+k,nml_n_1+j)+...
                        gauss_weight(i)*w(:,k)'*...
                        exp(-lambda_w(k)*theta(i))*dA1_deta*gauss_y(:,i);
                end;
            end;
        end;
        
        % derivative to T is zero
        
        % derivatives to representation points
        for k=1:s2
            for i=1:length(gauss_weight)
                for teller=1:m+1
                    J(nml+2*n+s1*n+s1+s2*n+s2+k,...
                        (K-1+m*floor((i-1)/m)+teller-1)*n+(1:n))=...
                        J(nml+2*n+s1*n+s1+s2*n+s2+k,...
                        (K-1+m*floor((i-1)/m)+teller-1)*n+(1:n))+...
                        gauss_weight(i)*w(:,k)'*exp(-lambda_w(k)*theta(i))*A*Pa(i,teller);
                end;
            end;
        end;
        
    end;
    
end;

% continuity condition in the initial point

res(nml+2*n+s1*n+s1+s2*n+2*s2+1:nml+2*n+s1*(n+1)+s2*n+2*s2+n,1)=...
    profile(:,1)-x1;
for k=1:s1
    res(nml+2*n+s1*n+s1+s2*n+2*s2+1:nml+2*n+s1*(n+1)+s2*n+2*s2+n,1)=...
        res(nml+2*n+s1*n+s1+s2*n+2*s2+1:nml+2*n+s1*(n+1)+s2*n+2*s2+n,1)-...
        epsilon*alpha(k)*v(:,k);
end;

% derivative to x1

J(nml+2*n+s1*(n+1)+s2*(n+2)+1:nml+2*n+s1*(n+1)+s2*(n+2)+n,...
    nml_n_1+nb_par+1:nml_n_1+nb_par+n)=-eye(n);

% derivative to alpha(k)

J(nml+2*n+s1*(n+1)+s2*(n+2)+1:nml+2*n+s1*(n+1)+s2*(n+2)+n,...
    nml_n_1+nb_par+2*n+s1*n+s1+s2*n+s2+1:...
    nml_n_1+nb_par+2*n+s1*n+s1+s2*n+s2+s1)=-epsilon*v;

% derivative to v(:,k)

for k=1:s1
    
    J(nml+2*n+s1*(n+1)+s2*(n+2)+1:nml+2*n+s1*(n+1)+s2*(n+2)+n,...
        nml_n_1+nb_par+2*n+(k-1)*n+1:nml_n_1+nb_par+2*n+k*n)=...
        -eye(n)*alpha(k)*epsilon;
    
end;

% derivative to y(0)

J(nml+2*n+s1*(n+1)+s2*(n+2)+(1:n),1:n)=eye(n);

% phase condition:

if ph && non_gauss,
    for l_i=1:l
        index_a=(l_i-1)*m+1;
        for k=1:m
            fac=gauss_abs(k)*(t((l_i-1)*m+1)-t(l_i*m+1));
            dPa=poly_dla(t(index_a:index_a+m),gauss_c(k));
            Pa=poly_lgr(t(index_a:index_a+m),gauss_c(k));
            u_prime_previous=previous.profile(:,index_a:index_a+m)*dPa';
            for q=1:m+1
                J(nml+2*n+s1*(n+1)+s2*(n+2)+n+1,...
                    (l_i-1)*m*n+1+(q-1)*n:(l_i-1)*m*n+q*n)= ...
                    J(nml+2*n+s1*(n+1)+s2*(n+2)+n+1,...
                    (l_i-1)*m*n+1+(q-1)*n:(l_i-1)*m*n+q*n) + ...
                    fac*Pa(q)*u_prime_previous';
            end
        end
    end
end

if ph
    res(nml+2*n+s1*(n+1)+s2*(n+2)+n+1,1)=0;
end;

% extra condition on the alpha(k)'s (their squares have to sum
%  up to 1 -- with a fixed epsilon, this adds the unknown T)
res(nml+2*n+s1*n+s1+s2*n+s2+s2+n+2,1)=alpha'*alpha-1;
J(nml+2*n+s1*n+s1+s2*n+s2+s2+n+2,...
    nml_n_1+nb_par+2*n+s1*(n+1)+s2*(n+1)+1:...
    nml_n_1+nb_par+2*n+s1*(n+1)+s2*(n+1)+s1)=2*alpha';

end
