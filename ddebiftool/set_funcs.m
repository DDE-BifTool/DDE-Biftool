function funcs=set_funcs(varargin)
%% fill in funcs structure with user-defined functions for use with DDE-Biftool functions
%
% possible named arguments: 'sys_rhs' (mandatory),'sys_ntau', 'sys_cond',
% 'sys_deri', 'sys_dtau', 'sys_mfderi', 'sys_dirderi', 'x_vectorized' (logical)
%
% Examples
% for DDE x'=-x(t-tau)+a*x^3-x^5 with parameters [tau,a]:
% funcs=set_funcs(...
%   'sys_rhs',@(x,p)-x(1,2)+p(2)*x(1,1)^3-x(1,1)^5,...
%   'sys_tau',@()1)
% uses default dde_gen_deriv for partial derivatives
%
% $Id: set_funcs.m 369 2019-08-27 00:07:02Z jansieber $
%
%% Process options
defaults={...
    'sys_rhs',[],...             % user-defined r.h.s.
    'sys_ntau',@()0,...          % number of delays (SD-DDEs only)
    'sys_tau',[],...             % index of delays as parameters or dependence of delays
    'sys_cond',@dummy_cond,...   % extra conditions
    'sys_deri',[],...            % Jacobians and second-order derivatives
    'sys_dtau',[],...            % Jacobians and second-order derivatives for delay (SD-DDEs only)
    'sys_mfderi',{},...          % higher-order derivatives (up to 5) for constant delay stst normal forms
    'sys_dirderi',[],...         % directional derivatives of sys_rhs
    'sys_dirdtau',[],...         % directional derivatives of sys_tau
    'x_vectorized',false,...     % can r.h.s, tau and their Jacobians be called with x=n x (nd+1) x nvec?
    'p_vectorized',false,...     % --"-- with x=n x (nd+1) x nvec and p= 1 x np x nvec?
    'hjac',@(ord)eps^(1/(2+ord)),... % deviation to be used for finite differences in numerical Jacobians
    'sys_cond_reference',false,...% does sys_cond need second (reference) input?
    'lhs_matrix',[]};            % for DDAEs an optional lhs_matrix ca nbe provided
if isstruct(varargin{1})
    funcs=dde_set_options(varargin{1},varargin(2:end));
else
    funcs=dde_set_options(defaults,varargin);
end
if isempty(funcs.sys_rhs)
    file=exist('sys_rhs','file');
    if file==2
        funcs.sys_rhs=@sys_rhs;
    else
        error('sys_rhs undefined');
    end
end
if isempty(funcs.sys_tau)
    file=exist('sys_tau','file');
    if file==2
        funcs.sys_tau=@sys_tau;
    else
        error('sys_tau undefined');
    end
end
%% convert user-provided (or empty) lhs matrix into function
% we do not know the dimension of the problem if the user does not specify
% lhs_matrix. Hence, funcs.lhs_matrix is a function that takes as its
% argument the dimension of x.
funcs.lhs_matrix=@(sz)lhs_matrix(sz,funcs.lhs_matrix);
%% test for state-dependent delay
funcs.tp_del=true;
try
    dummytau=funcs.sys_tau(); %#ok<NASGU>
    funcs.tp_del=false;
catch %#ok<CTCH>
    funcs.tp_del=true;
end
%% pseudo-vectorize and make sure that the output has expected shape
funcs.sys_rhs=@(x,p)dde_wrap_rhs(x,p,funcs.sys_rhs,funcs.x_vectorized,funcs.p_vectorized);
if funcs.tp_del
    funcs.sys_tau=@(itau,x,p)dde_wrap_rhs(x,p,funcs.sys_tau,funcs.x_vectorized,funcs.p_vectorized,itau);
end
%% if sys_deri is provided and vectorization is on, wrap output to ensure right format
if any(strcmp('sys_deri',varargin))||...
        (isfield(funcs,'sys_deri_provided') && funcs.sys_deri_provided)
        if ~isfield(funcs,'sys_deri_provided')
            funcs.sys_deri_provided=2;
        end
        funcs.sys_deri=@(x,p,nx,np,v)dde_wrap_deri(x,p,nx,np,v,funcs.sys_deri,...
            funcs.x_vectorized,funcs.p_vectorized);
else
    funcs.sys_deri_provided=0;
end
%% Directional derivatives provided?
if any(strcmp('sys_dirderi',varargin))
    funcs.sys_dirderi_provided=length(funcs.sys_dirderi);
    for order=funcs.sys_dirderi_provided:-1:1
        funcs.sys_dirderi{order}=@(xx,p,dx,dp)dde_wrap_dirderi(xx,p,dx,dp,funcs.sys_dirderi{order},...
            funcs.x_vectorized,funcs.p_vectorized);
    end
end
if ~isfield(funcs,'sys_dirderi_provided')
    funcs.sys_dirderi_provided=0;
end
%% if not enough derivatives provided, add finite differences
if funcs.sys_dirderi_provided<2 && ~funcs.sys_deri_provided
    if funcs.sys_dirderi_provided>0
        %% derivative orders are given as cell array of functions up to some orer
        derivbase=@(xx,p,dxx,dp)funcs.sys_dirderi{end}(xx,p,dxx,dp);
    elseif funcs.sys_dirderi_provided==0
        derivbase=@(xx,p,dxx,dp)funcs.sys_rhs(xx,p);
    end
    ldirderi=length(funcs.sys_dirderi);
    for order=ldirderi+1:2
        funcs.sys_dirderi{order}=@(x,p,dx,dp)dde_dirderiv(derivbase,...
        x,p,dx,dp,order-ldirderi,'nf',size(x,1),'hjac',funcs.hjac(order));
    end
end
%% funcs.sys_deri not provided (but possibly dirderi)
% combine Jacobians from directional derivatives
if funcs.sys_deri_provided<2
    funcs.sys_deri=...
        @(x,p,nx,np,v)dde_gen_deriv(funcs.sys_dirderi,x,p,nx,np,v,1);
end
%% dervatives of tau only for sd-DDEs
if funcs.tp_del
    %% if sys_dtau is provided and vectorization is on, wrap output to ensure right format
    if  any(strcmp('sys_dtau',varargin)) ||...
            (isfield(funcs,'sys_dtau_provided') && funcs.sys_dtau_provided)
        if ~isfield(funcs,'sys_dtau_provided')
            funcs.sys_dtau_provided=2;
        end
        funcs.sys_dtau=@(itau,x,p,nx,np)dde_wrap_deri(x,p,nx,np,[],funcs.sys_dtau,...
            funcs.x_vectorized,funcs.p_vectorized,itau);
    else
        funcs.sys_dtau_provided=0;
    end
    %% Directional derivatives provided?
    if any(strcmp('sys_dirdtau',varargin))
        funcs.sys_dirdtau_provided=length(funcs.sys_dirdtau);
        for order=funcs.sys_dirdtau_provided:-1:1
            funcs.sys_dirdtau{order}=...
                @(it,xx,p,dx,dp)dde_wrap_dirderi(xx,p,dx,dp,funcs.sys_dirdtau{order},...
                funcs.x_vectorized,funcs.p_vectorized,it);
        end
    end
    if ~isfield(funcs,'sys_dirdtau_provided')
        funcs.sys_dirdtau_provided=0;
    end
    %% if not enough derivatives provided, add finite differences
    if funcs.sys_dirdtau_provided<2 && ~funcs.sys_dtau_provided
        if funcs.sys_dirdtau_provided>0
            %% derivative orders are given as cell array of functions up to some orer
            derivbase=@(it,xx,p,dxx,dp)funcs.sys_dirdtau{end}(it,xx,p,dxx,dp);
        elseif funcs.sys_dirdtau_provided==0
            derivbase=@(it,xx,p,dxx,dp)funcs.sys_tau(it,xx,p);
        end
        %% not enough derivatives provided
        ldirdtau=length(funcs.sys_dirdtau);
        for order=ldirdtau+1:2
            funcs.sys_dirdtau{order}=@(itau,x,p,dx,dp)dde_dirderiv(...
                @(xx,p,dxx,dp)derivbase(itau,xx,p,dxx,dp),x,p,dx,dp,order-ldirdtau,...
                'nf',1,'hjac',funcs.hjac(order));
        end
    end
    %% sys_dtau not provided (but possibly dirdtau)
    if funcs.sys_dtau_provided<2
        dfunc=@(itau){...
            @(xa,pa,dx,dp)funcs.sys_dirdtau{1}(itau,xa,pa,dx,dp),...
            @(xa,pa,dx,dp)funcs.sys_dirdtau{2}(itau,xa,pa,dx,dp)};
        funcs.sys_dtau=@(itau,x,p,nx,np)dde_gen_deriv(dfunc(itau),x,p,nx,np,[]);
    end
    if isfield(funcs,'hjac')
        funcs=rmfield(funcs,'hjac');
    end
end
end
%% dummy condition
function [resi,condi]=dummy_cond(point,pref) %#ok<INUSD>
resi=zeros(0,1);
condi=repmat(point,0,1);
end
%% uniform format for lhs_matrix
function L=lhs_matrix(sz,arg)
if isempty(arg)
    L=eye(sz);
else
    L=arg;
end
end
