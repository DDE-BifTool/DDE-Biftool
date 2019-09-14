%% set_rotfuncs - Fill in funcs structure for use with DDE-Biftool, rotational symmetric case
%%
function funcs=m_set_rotfuncs(varargin)
%% Named arguments
%
% * |'sys_rhs'| (mandatory), 
% * |'rotation'| (mandatory),
% * |'exp_rotation'| (mandatory),
% * |'sys_ntau'|, 
% * |'sys_cond'|, (DOES NOT SUPPORT sys_cond UNTIL m_rot_cond IS UPDATED)
% * |'sys_deri'| (unused),
% * |'sys_dtau'| (unused).
%
% See rotsym_demo for example usage
% uses defaults df_deriv for partial derivatives
%
% (c) DDE-BIFTOOL v. 3.1.1(29), 15/04/2014
%

%% Process options
defaults={'rotation',[],'exp_rotation',[],'rot_tol',1e-8,'hjac',1e-6};
[options,pass_on]=dde_set_options(defaults,varargin,'pass_on');
rot_tol=options.rot_tol;
funcs=set_funcs(pass_on{:});
%% rotation matrices
% Matrix rotation
if isempty(options.rotation)
    error('rotation matrix must be given as argument ''rotation''');
else
    funcs.rotation=options.rotation;
end

% Exp_Matrix rotation
if isempty(options.exp_rotation) && isa(funcs.rotation,'double')
    funcs.exp_rotation=@(phi)expm(funcs.rotation*phi);
elseif isempty(options.exp_rotation) && isa(funcs.rotation,'cell')
    error(['When using multiple rotation matricies, ', ...
        'exp rotation matrix must be given as argument ''exp_rotation''']);
else
    funcs.exp_rotation=options.exp_rotation;
end

%% Compatibility Test

if isa(funcs.rotation,'double')
    % test comaptibility using rot_sym method for a SINGLE rotation matrix.
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
elseif isa(funcs.rotation,'cell')
    % test compatibility given MULTIPLE rotation matricies.
    warning(['Compatibility testing is currently not supported for ', ...
        'multiple rotation matricies.'])
end


%% modify sys_rhs to include shift to rotating coordinates
funcs.orig_rhs=funcs.sys_rhs;
funcs.sys_rhs = @(xx,p)m_rot_rhs(xx,p, ...
    funcs.rotation,funcs.exp_rotation, ...
    funcs.orig_rhs,funcs.sys_tau,funcs.x_vectorized);
%% set derivatives to defaults if not given by user
if ~funcs.sys_deri_provided;
    funcs.orig_deri=funcs.sys_deri;
    funcs.sys_deri=@(x,p,nx,np,v)df_deriv(funcs,x,p,nx,np,v,options.hjac);
else
    %funcs.orig_deri=funcs.sys_deri;
    %funcs.sys_deri=@(x,p,nx,np,v)rot_deriv(funcs,x,p,nx,np,v);
    error('No support yet for passing on user-defined derivatives.');
end
%% modify sys_cond to include fixing of rotational phase
funcs.orig_cond=funcs.sys_cond;
% funcs.sys_cond=@(p)rot_cond(p,funcs.rotation,funcs.orig_cond);
% Notice that the new version of m_rot_cond does not currently support user
% defined extra conditions!
funcs.sys_cond = @(p)m_rot_cond(p,funcs.rotation);
%% state-dependent delays currently not supported
if funcs.tp_del
    error('State-dependent delays are not supported.')
end
if isempty(funcs.sys_dtau) % unused
    funcs.sys_dtau=@(x,p,nx,np,v)df_derit(funcs,x,p,nx,np,v);
end
end
