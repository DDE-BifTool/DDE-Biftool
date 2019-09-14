function [J,res,tT,extmesh,Jstruc]=psol_jac_simple(funcs,psol,free_par,varargin)
%% simplified residual & Jacobian of collocation problem for periodic orbits
% function [J,res,tT,extmesh]=psol_jac_simple(c,T,profile,t,deg,par,free_par,phase,varargin)
% INPUT:
%   funcs problem functions
%	c collocation parameters in [0,1]^m
%   psol: solution point (type psol)
%	free_par free parameters numbers in N^np
%   wrapJ (optional key-value pair, default true) 
%        wrap time points periodically into [0,1]
% OUTPUT:
%	J jacobian in R^(n*deg*l+n+s x n*deg*l+1+n+np)
%	res residual in R^(n*deg*l+n+s)
%   tT delays, scaled by period
%   extmesh mesh of time points, extended back to -max(tT(:))
%
% modified to permit arbitrary nesting of state-dependent delays,
% vectorisation and optional re-use for computation of Floquet multipliers
% and extra conditions for state-dependent delay (p_tau, etc)
%
% Optional inputs:
%
% If 'wrapJ' (default true) the returned jacobian is augmented with
% derivative wrt period and free_par and wrapped around periodically inside
% [0,1] if ~wrapJ the returned jacobian puts its entries into the interval
%  [-min(delay,0),max(1-[0,delays])], no augmentation is done. This
%  Jacobian can be used for computation of Floquet multipliers and modes
%
% If 'c_is_tvals' (default false) then the values in c (must be not empty) are
% taken as those values of time in [0,1] of the entire period, where
% residual and Jacobian are calculated.
%
% The argument 'Dtmat' (default=Id) is multiplied with the time
% derivative. Dtmat==zeros() permits algebraic equations.
%
% The argument 'bc' (default true) controls whether boundary conditions
% should be attached.

% (c) DDE-BIFTOOL v. 2.00, 30/11/2001
%
% $Id$
%
%% 
psol_prof=psol.profile;
tmesh=psol.mesh;
period=psol.period;
deg=psol.degree;
par=psol.parameter;
%% optional
default={'wrapJ',true,'c',[],'c_is_tvals',false,'Dtmat',eye(size(psol_prof,1)),'period',true};
options=dde_set_options(default,varargin,'pass_on');
%% define problem & problem dimensions
n=size(psol_prof,1);      % dimension of x
nf=size(options.Dtmat,1); % dimension of f
np=length(free_par);      % number of free parameters
%% constant delay
n_tau=funcs.sys_tau(); % delay numbers
tau=par(n_tau);  % delay values
tau=[0;tau(:)];  % incl tau=0
d=length(n_tau)+1; % number of delays
nint=(length(tmesh)-1)/deg; % number of collocation intervals
if ~options.c_is_tvals
    neqs=deg*nint;          % number of equations from collocation points
else
    neqs=length(options.c); % only evaluate res,Jac at times c
end
%% check array sizes:
if nint~=floor(nint)
    error('PSOL_JAC: t (length%d) does not contain %d intervals of %d points!',...
        length(tmesh),nint,deg);
end
if ~options.c_is_tvals && length(options.c)~=deg && ~isempty(options.c)
    error('PSOL_JAC: wrong number of collocation parameters length(c)=%d, m=%d!',...
        length(options.c),deg);
end
tcoarse=tmesh(1:deg:end); % boundaries of collocation intervals
h_int=diff(tcoarse); % lengths of collocation intervals
%% obtain collocation parameters and times
if ~options.c_is_tvals
    if isempty(options.c)
        c=poly_gau(deg);
        c=c(:)';
    else
        c=options.c(:)';
    end
    t_c=tcoarse(ones(length(c),1),1:end-1)+...
        c(ones(length(h_int),1),:)'.*h_int(ones(length(c),1),:);
else
    t_c=options.c;
end
t_c=t_c(:)'; % evaluation points for DE residuals
%% pre-allocate needed temporary arrays
Pb{1}=zeros(d,deg+1,neqs);  % Lagrange stamp, mapping delayed values onto profile
Pb{2}=Pb{1}; % Lagrange stamp, mapping delayed derivatives onto profile
Pb{3}=Pb{1}; % Lagrange stamp, mapping delayed 2nd derivatives onto profile
xx=zeros(n,d,neqs);  % array [x(t),x(t-tau1),...]
dtx=zeros(n,d,neqs);  % array [x'(t),x'(t-tau1),...]
ddtx=zeros(n,d,neqs);  % array [x''(t),x''(t-tau1),...]
c_tau=zeros(d,neqs); % time points t-tau(i-1) for current colloc point t
c_tau_mod=zeros(d,neqs); % time points t-tau(i-1) for current colloc point t (mod[0,1])
index_b=zeros(d,neqs); % starting indices of collocation intervals for t-tau(i-1)
index_b_mod=zeros(d,neqs); % starting indices of collocation intervals for t-tau(i-1) (mod[0,1])
h_int_del=zeros(d,neqs); % length of collocation interval in which t-tau(i-1) lies
c_tau_trans=zeros(d,neqs); % position of t-tau(i) in its colloc interval, rescaled to [0,1]
tT=tau(:)/period;
tT=repmat(tT,1,neqs);
%% for all collocation points find delays, delayed profiles and derivatives
% (vectorization in l_i & m_i possible)
for t_i=1:d
    c_tau(t_i,:)=t_c-tT(t_i,:);
    c_tau_mod(t_i,:)=mod(c_tau(t_i,:),1);
    [~,ibase,c_tau_trans(t_i,:)]=psol_eva(tmesh,tmesh,c_tau_mod(t_i,:),deg); 
    index_b_mod(t_i,:)=(ibase-1)*deg+1;
    index_b(t_i,:)=index_b_mod(t_i,:)+round(c_tau(t_i,:)-c_tau_mod(t_i,:))*nint*deg;
    h_int_del(t_i,:)=h_int(ibase);
    Pb{1}(t_i,:,:)=poly_elg(deg,c_tau_trans(t_i,:));
    Pb{2}(t_i,:,:)=poly_del(deg,c_tau_trans(t_i,:))./...
        h_int_del(t_i*ones(deg+1,1),:);
    Pb{3}(t_i,:,:)=poly_d2l(deg,c_tau_trans(t_i,:))./...
        h_int_del(t_i*ones(deg+1,1),:).^2;
    for k=1:neqs
        xx(:,t_i,k)=  psol_prof(:,index_b_mod(t_i,k)+(0:deg))*Pb{1}(t_i,:,k)';
        dtx(:,t_i,k)= psol_prof(:,index_b_mod(t_i,k)+(0:deg))*Pb{2}(t_i,:,k)';
        ddtx(:,t_i,k)=psol_prof(:,index_b_mod(t_i,k)+(0:deg))*Pb{3}(t_i,:,k)';
    end
end
%% determine index ranges for wrapped or unwrapped Jacobian
if ~options.wrapJ % compute Jde_dx only for Floquet multipliers
    indmin=min(index_b(:));
    indshift=1-indmin;
    indmax=max(index_b(:))+deg;
    doJcomb=false;
    % set up extended mesh
    indrg=indmin:indmax;
    t_ind=mod(indrg-1,nint*deg)+1;
    t_shift=floor((indrg-1)/(nint*deg));
    extmesh=tmesh(t_ind)+t_shift;
    index_b=index_b+indshift;
else % compute augmented jacobian
    indshift=0;
    indmax=deg*nint+1;
    index_b=index_b_mod;
    doJcomb=true;
    extmesh=tmesh;
end
%% init J, res:
Jstruc.cols.x=(0:indshift+indmax-1)*n+1;
Jstruc.cols.xrg=1:Jstruc.cols.x(end)+n-1;
Jstruc.rows.de=(0:neqs-1)*nf+1;
Jstruc.cols.nx=n;
Jstruc.rows.nf=nf;
Jstruc.cols.np=np;
Jstruc.free_par=free_par;
[~,par_is_delay]=ismember(free_par,n_tau);
[~,delay_is_par]=ismember([0,n_tau],free_par);
Jstruc.delay_is_par=delay_is_par;
Jstruc.par_is_delay=par_is_delay;
resde=zeros(nf,neqs);
Jde_dx=zeros(nf,neqs,n,indshift+indmax);
if doJcomb
    if options.period
        Jde_dT=zeros(nf,neqs);
        Jstruc.cols.T=Jstruc.cols.xrg(end)+1;
        Jstruc.cols.nT=1;
    else
        Jstruc.cols.T=[];
        Jstruc.cols.nT=0;
    end
    if np>0
        Jde_dp=zeros(nf,neqs,np);
        Jstruc.cols.p=Jstruc.cols.xrg(end)+Jstruc.cols.nT+(1:np);
    else
        Jstruc.cols.p=[];
    end
else
    Jstruc.cols.T=[];
    Jstruc.cols.nT=0;
    Jstruc.cols.p=[];
    Jstruc.cols.np=0;
end
%% obtain all values of sys_rhs, sys_deriv and sys_dxtau
vals=psol_sysvals(funcs,xx,par,free_par,'fdim',nf);
%%
for i=1:neqs
    %% insert all derivatives into large Jacobians
    %% add dtx for x' in Jx and res:
    resde(:,i)=options.Dtmat*dtx(:,1,i);
    Jde_dx(:,i,:,index_b(1,i)+(0:deg))=Jde_dx(:,i,:,index_b(1,i)+(0:deg))...
        +reshape(kron(Pb{2}(1,:,i),options.Dtmat),[nf,1,n,deg+1]);
    if Jstruc.cols.np>0
        %% add parameter in Jde_dp:
        for p_i=1:np
            Jde_dp(:,i,p_i)=Jde_dp(:,i,p_i)-period*vals.dfdp(:,p_i,i);
        end
    end
    if Jstruc.cols.nT>0
        %% add -f for dT in J:
        Jde_dT(:,i)=Jde_dT(:,i)-vals.f(:,i);
    end
    %% add Tf in res:
    resde(:,i)=resde(:,i)-period*vals.f(:,i);
    %% delayed values (incl tau=0): insert Jacobians into J
    for t_i=1:d
        %% add -T*Pb{1}*dfdx(:,t_i)
        Jde_dx(:,i,:,index_b(t_i,i)+(0:deg))=...
            Jde_dx(:,i,:,index_b(t_i,i)+(0:deg))-...
            reshape(kron(Pb{1}(t_i,:,i),period*vals.dfdx(:,:,t_i,i)),[nf,1,n,deg+1]);
        if Jstruc.cols.nT>0
            %% add -T*A1*sum b*Pb{2}*dc_tau for dT in J:
            Jde_dT(:,i)=Jde_dT(:,i)-vals.dfdx(:,:,t_i,i)*dtx(:,t_i,i)*tT(t_i,i);
        end
        %% parameter derivative if delay is free par
        if Jstruc.cols.np>0 && delay_is_par(t_i)>0
            Jde_dp(:,i,delay_is_par(t_i))=Jde_dp(:,i,delay_is_par(t_i))+...
                vals.dfdx(:,:,t_i,i)*dtx(:,t_i,i);
        end
    end
end
%% assemble overall residual & Jacobian
Jde_dx=reshape(Jde_dx,[nf*neqs,n*(indshift+indmax)]);
res=resde(:);
J=Jde_dx;
if Jstruc.cols.nT>0
    J=[J,Jde_dT(:)];
end
if Jstruc.cols.np
    J=[J,reshape(Jde_dp,[nf*neqs,np])];
end
end