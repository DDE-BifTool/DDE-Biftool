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
% uses defaults df_deriv and df_derit for partial derivatives and no extra conditions
%
% $Id: set_funcs.m 190 2017-04-25 09:15:51Z jansieber $
%
%% Process options
defaults={...
    'sys_rhs',[],...             % user-defined r.h.s.
    'sys_ntau',@()0,...          % number of delays (SD-DDEs only)
    'sys_tau',[],...             % index of delays as parameters or dependence of delays
    'sys_cond',@dummy_cond,...   % extra conditions
    'sys_deri',[],...            % Jacobians and second-order derivatives
    'x_vectorized',false,'maxorder',5};
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
%% for vectorization make sure that the output has expected shape
if funcs.x_vectorized
    funcs.sys_rhs=@(x,p)reshape(funcs.sys_rhs(x,p),[size(x,1),1,size(x,3)]);
end
%% if sys_deri is provided and vectorization is on, wrap output to ensure right format
if any(strcmp('sys_deri',varargin)) ||...
        (isfield(funcs,'sys_deri_provided') && funcs.sys_deri_provided)
        funcs.sys_deri_provided=true;
        if funcs.x_vectorized
            funcs.sys_deri=@(x,p,nx,np,v)wrap_deri(x,p,nx,np,v,funcs.sys_deri);
        end
else
    funcs.sys_deri_provided=false;
end
%% if sys_deri is not provided
if ~funcs.sys_deri_provided 
    %% but directional derivatives are in sys_dirderi
    % compose sys_deri of dde_gen_deriv and sys_dirderi, otherwise of
    % dde_gen_deriv and dde_dirderiv
    if any(strcmp('sys_dirderi',varargin))
        funcs.sys_deri_provided=true;
        %% if derivative orders are given as cell array of functions, convert
        if iscell(funcs.sys_dirderi)
            funcs.sys_dirderi=cell2derifunc(funcs.sys_dirderi,funcs.maxorder,funcs.x_vectorized);
        end
    else
        %% no derivatives provided at all
        funcs.sys_dirderi=@(order,x,p,dx,dp)dde_dirderiv(funcs.sys_rhs,...
            x,p,dx,dp,order,'nf',size(x,1),'x_vectorized',funcs.x_vectorized);
    end
    funcs.sys_deri=@(x,p,nx,np,v)dde_gen_deriv(funcs.sys_dirderi,x,p,nx,np,v);
end
if isfield(funcs,'maxorder')
    funcs=rmfield(funcs,'maxorder');
end
end
%% dummy condition
function [resi,condi]=dummy_cond(point) %#ok
resi=[];
condi=[];
end
%% convert cell of directional derviatives to single functions,
% appending finite-difference approximations 
function dirderi=cell2derifunc(dirderi_cell_inp,maxorder,x_vectorized)
nderivs=length(dirderi_cell_inp);
if ~x_vectorized
    for i=nderivs:-1:1
        dirderi_cell{i}=@(xx,p,dx,dp)wrap_dirderi(xx,p,dx,dp,dirderi_cell_inp{i});
    end
else
    dirderi_cell=dirderi_cell_inp;
end
if nderivs<maxorder
    %% if only first derivatives are given append numerical approximations
    for i=nderivs+1:maxorder
        dirderi_cell{i}=@(xx,p,dx,dp)...
            dde_dirderiv(@(xa,pa)dirderi_cell{nderivs}(xa,pa,dx,dp),xx,p,dx,dp,i-nderivs,...
            'x_vectorized',true);
    end
end
dirderi=@(order,xx,p,dx,dp)dirderi_cell{order}(xx,p,dx,dp);
end
%%

%% reshape deri for vectorization
function J=wrap_deri(x,p,nx,np,v,sys_deri)
n=size(x,1);
nvec=size(x,3);
J=sys_deri(x,p,nx,np,v);
J=reshape(J,[n,numel(J)/(n*nvec),nvec]);
end
%% pseudo-vectorize sys_dirderi
function J=wrap_dirderi(x,p,dx,dp,dirderi)
xsize=size(x);
vecdim=[xsize(3:end),1];
x=reshape(x,xsize(1),xsize(2),[]);
dx=reshape(dx,size(x));
for i=size(x,3):-1:1
    J(:,i)=dirderi(x(:,:,i),p,dx(:,:,i),dp);
end
J=reshape(J,[xsize(1),vecdim]);
end
