function [ubp,suc,v0,w0,ulin]=dde_bp_solve(f,u0,varargin)
% solve for branch point accurately using Newton iteration of extended
% system:
%
% $$ 0=f(u+v0 a)+w0 b, v0'(u-u0)=0, 0=df(u)[v1,v2], v0'[v1,v2]=-Id_2 $$
%
% Variables are $u,v1,v2$ (n+1=nu dimensional), $a$ (two-dim), $b$
% (one-dim). v0 (n+1 x 2) is approx Nullspace of df(u_*) and w0 is approx
% nullspace of df(u_*0'
%
% Inputs: 
%
% * f (n+1-dimensional input), returning r (n dim), J (n x (n+1) dim) 
% * u0 (n+1) dim: initial guess
% * solver: Newton type solver, accepting [r,J]=F(u), u0, returning u_* and
% success flag
% * hdev: deviation used for finite-difference approximation of derivative
% of Jacobian 
%
% $Id: dde_bp_solve.m 324 2019-02-03 01:54:58Z jansieber $
mth=getfield(df_mthod('fold'),'point');
default={'solver',[],'hdev',1e-4,'print',0,'v0',[],'w0',[]};
options=dde_set_options(default,varargin,'pass_on');
mth.print_residual_info=options.print;
if isempty(options.solver)
    options.solver=@(fa,xa)dde_nsolve(fa,xa,mth);
end
%% initialize initial guess and residual function
nu=length(u0);
if isempty(options.v0) || isempty(options.w0)
    %% initial Jacobian
    [r0,J0]=f(u0); %#ok<ASGLU>
    %% its approximate nullspaces (right: 2d, left: 1d)
    [v0,w0]=dde_svdspaces_lr(J0,2);
else
    v0=options.v0;
    w0=options.w0;
end
%% regularized f: [f(u+v0 a)+w b; v0'(u-u0)] has dimension deficit 2 and its Jacobian
reg=@(u,action)f_reg_residual(u,u0,v0,w0,f,action);
f_reg=@(u)reg(u,'res');
J_reg=@(u)reg(u,'J');
%% r.h.s for BP problem
fbph=@(uex)f_bp(uex,nu,f_reg,J_reg,options.hdev);
%% initial guess
uex0=[u0;zeros(3,1)];
J0=J_reg(uex0);
Jub=J0(:,[1:end-3,end]);
Jva=J0(:,end-2:end-1);
uex12=-Jub\Jva;
uexini=[uex0(:);reshape(uex12(1:end-1,:),[],1)];
%% apply Newton iteration
[uexbp,suc]=options.solver(fbph,uexini);
ubp=reg(uexbp(1:nu+3),'uend');
ulin=reshape(uexbp(nu+3+(1:2*nu)),nu,2);
[r0,J0]=f(ubp); %#ok<ASGLU>
[v0,w0]=dde_nullspaces_lr(J0);
end
%% residual and Jacobian of defining system for BP
function [re,Je]=f_bp(uex,nu,f_reg,J_reg,hjac)
u0=uex(1:nu);
alpha=uex(nu+(1:2));
beta=uex(nu+3);
u12=reshape(uex(nu+3+(1:2*nu)),nu,2);
u0ext=[u0;alpha;beta];
u12ext=[u12;eye(2);zeros(1,2)];
[r0,Jex]=f_reg(u0ext);
r12=Jex*u12ext;
re=[r0;r12(:)];
%% Jacobian
nex=length(r0);
J0=[Jex,sparse(nex,2*nu)];
dJdu1=(J_reg(u0ext+hjac*u12ext(:,1))-J_reg(u0ext-hjac*u12ext(:,1)))/(2*hjac);
dJdu2=(J_reg(u0ext+hjac*u12ext(:,2))-J_reg(u0ext-hjac*u12ext(:,2)))/(2*hjac);
J1=[dJdu1, Jex(:,1:nu),    sparse(nex,nu)];
J2=[dJdu2, sparse(nex,nu), Jex(:,1:nu)];
Je=cat(1,J0,J1,J2);
end
%% residual and derivative of extended (regularised f)
function [re,Je]=f_reg_residual(uex,u0,v0,w0,f,output)
u=uex(1:end-3);
alpha=uex(end-2:end-1);
beta=uex(end);
[r,J]=f(u+v0*alpha);
re=[r+w0*beta;v0'*(u-u0)];
Je=cat(1,[J,J*v0,w0],[v0',zeros(2,3)]);
switch output
    case 'J'
        [Je,re]=deal(re,Je);
    case 'uend'
        re=u+v0*alpha;
end
end