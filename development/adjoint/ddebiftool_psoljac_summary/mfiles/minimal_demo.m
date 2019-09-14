%% minimal demo of simplified Jacobian for DDEs (w constant delays)
clear
addpath('ddebiftool');
%%
ntau=1;
parnames={'tau','a','b','d'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
%% Duffing oscillator
% The example is the Duffing oscillator with delayed feedback discussed in
% the large-delay limit by Yanchuk & Perlikowski in (PRE79,0462211,2009):
%
% $$x''(t)+d*x'(t)+a*x(t)+x^3+b*[x(t)-x(t-\tau)]=0$$
%
% The parameters are $(\tau,a,b,d)$ (used in this order in the |parameter|
% vector).
rhs=@(x,p)[x(2,1,:);...
    -p(ind.d)*x(2,1,:)-p(ind.a)*x(1,1,:)-x(1,1,:).^3-p(ind.b)*(x(1,1,:)-x(1,2,:))];
funcs=set_funcs('sys_rhs',rhs,'sys_tau',@()1,'x_vectorized',true);
%% load family of periodic orbits continued in delay
s=load('minimal_demo_results.mat','per_orb');
per_orb=s.per_orb;
clear s
%% Demo of W matrices
pt=per_orb.point(100);
for tau=0:0.01:pt.period
    W=dde_coeffmat_tau_simple(tau,pt);
    spy(W{1});
    drawnow
end
%% check Jacobians and residuals along family
% parameters with respect to which derivatives are computed. Note that
% ScJacobian is a simple finite-difference routine to create
% (non-vectorized) Jacobians.
free_par=[ind.tau,ind.a];
% conversion fro mpoint to x vector
x_from_psol=@(pt)[pt.profile(:);pt.period;pt.parameter(free_par)'];
for i=1:length(per_orb.point)
    [J1,r1]=dde_jac_simple(funcs,per_orb.point(i),free_par);
    assert(norm(r1,'inf')<1e-8);
    template=per_orb.point(i);
    p_from_x=@(x)psol_from_x(x,template,free_par);
    J0=ScJacobian(@(x)dde_jac_simple(funcs,p_from_x(x),free_par,'output','residual'),...
        x_from_psol(per_orb.point(i)),1e-4);
    fprintf('i=%d of %d: norm(J0-J1)=%g\n',i,length(per_orb.point),norm(J0-J1,'inf'));
end

