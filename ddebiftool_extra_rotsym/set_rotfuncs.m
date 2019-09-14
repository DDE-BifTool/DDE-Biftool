%% set_rotfuncs - Fill in funcs structure for use with DDE-Biftool, rotational symmetric case
%%
function funcs=set_rotfuncs(varargin)
%% Named arguments
%
% * |'sys_rhs'| (mandatory), 
% * |'rotation'| (mandatory),
% * |'exp_rotation'| (mandatory),
% * |'sys_ntau'|, 
% * |'sys_cond'|,
% * |'sys_deri'| (unused),
% * |'sys_dtau'| (unused).
%
% See rotsym_demo for example usage. Uses default finite differences for
% partial derivatives of order>1.
%
% $Id: set_rotfuncs.m 369 2019-08-27 00:07:02Z jansieber $
%

%% Process options
defaults={'rotation',[],'exp_rotation',[],'rot_tol',1e-8,'hjac',...
    @(ord)eps^(1/(2+ord))};
[options,pass_on]=dde_set_options(defaults,varargin,'pass_on');
rot_tol=options.rot_tol;
funcs=set_funcs(pass_on{:},'hjac',options.hjac,'p_vectorized',false);
funcs.p_vectorized=false;
%% check rotation matrices
if isempty(options.rotation)
    error('rotation matrix must be given as argument ''rotation''');
else
    funcs.rotation=options.rotation;
end
if isempty(options.exp_rotation)
    funcs.exp_rotation=@(phi)expm(funcs.rotation*phi);
else
    funcs.exp_rotation=options.exp_rotation;
end
% test comaptibility
phi=linspace(0,2*pi,10);
for i=1:length(phi)
    err=expm(funcs.rotation*phi(i))-funcs.exp_rotation(phi(i));
    if max(abs(err))>rot_tol
        error('exp(rot*phi)~=exp_rotation(phi) for phi=%g',phi(i));
    end
end
err=funcs.exp_rotation(2*pi)-eye(size(funcs.rotation));
if max(abs(err))>rot_tol
    error('exp(rot*2*pi)~=Id');
end
err=funcs.rotation+funcs.rotation.';
if max(abs(err))>rot_tol
    error('rot~=-rot^T');
end
funcs.lhs_num=funcs.lhs_matrix(size(funcs.rotation,1));
if norm(funcs.lhs_num*funcs.rotation-funcs.rotation*funcs.lhs_num,'inf')>rot_tol
    error('lhs_matrix and rotation do not commute');
end
%% modify sys_rhs to include shift to rotating coordinates
funcs.orig_rhs=funcs.sys_rhs;
rot_sys_rhs=@(xx,p)rot_rhs(xx,p,funcs.rotation,funcs.exp_rotation,...
    funcs.orig_rhs,funcs.sys_tau,funcs.lhs_num);
funcs.sys_rhs=@(x,p)dde_wrap_rhs(x,p,rot_sys_rhs,funcs.x_vectorized,funcs.p_vectorized);
%% set derivatives to defaults if not given by user
dfbase=@(xx,p,dxx,dp)funcs.sys_rhs(xx,p);
base_order=0;
orig_deri_provided=funcs.sys_deri_provided;
orig_dirderi_provided=funcs.sys_dirderi_provided;
if funcs.sys_deri_provided || funcs.sys_dirderi_provided>1
    funcs.sys_dirderi_provided=1;
    funcs.orig_dirderi{1}=@(xx,p,dx,dp)...
        dde_wrap_dirderi(xx,p,dx,dp,funcs.sys_dirderi{1},...
        funcs.x_vectorized,funcs.p_vectorized);
    dirderi=@(xx,p,dx,dp)rot_dirderi(xx,p,dx,dp,...
        funcs.rotation,funcs.exp_rotation,funcs.orig_dirderi{1},...
        funcs.sys_tau,funcs.lhs_num);
    funcs.sys_dirderi{1}=@(xx,p,dx,dp)...
        dde_wrap_dirderi(xx,p,dx,dp,dirderi,...
        funcs.x_vectorized,funcs.p_vectorized);
    dfbase=funcs.sys_dirderi{1};
    base_order=1;
    funcs.sys_deri_provided=1;
    if orig_deri_provided ||orig_dirderi_provided>1
        warning('set_rotfuncs:sys_deri',...
            ['set_rotfuncs: No support yet for ',...
            'passing on user-defined derivatives of degree>1.']);
    end
end
funcs.add_deriv=@add_deriv;
funcs=add_deriv(funcs,dfbase,base_order);
%% modify sys_cond to include fixing of rotational phase
funcs.orig_cond=funcs.sys_cond;
funcs.orig_cond_reference=funcs.sys_cond_reference;
funcs.sys_cond_reference=true;
funcs.sys_cond=@(p,pref)rot_cond(p,funcs,pref);
%% state-dependent delays currently not supported
if funcs.tp_del
    error('State-dependent delays are not supported.')
end
if isempty(funcs.sys_dtau) % unused
    %funcs.sys_dtau=@(x,p,nx,np,v)df_derit(funcs,x,p,nx,np,v);
end
end
%%
function funcs=add_deriv(funcs,derivbase,base_order)
for order=base_order+1:2
    funcs.sys_dirderi{order}=@(x,p,dx,dp)dde_dirderiv(derivbase,...
        x,p,dx,dp,order-base_order,'nf',size(x,1),'hjac',funcs.hjac(order-base_order));
end
funcs.sys_deri=...
    @(x,p,nx,np,v)dde_gen_deriv(funcs.sys_dirderi,x,p,nx,np,v,1);
end
