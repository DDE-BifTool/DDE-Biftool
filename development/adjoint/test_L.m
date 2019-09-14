%% test matrix mapping profile to list of delayed values
% clear
% load('minimal_demo_results.mat')
pt=per_orb.point(3);
f1=set_funcs('sys_rhs',@(xx,p)xx(:,1,:)+xx(:,2,:),'sys_tau',@()funcs.sys_tau(),'x_vectorized',true);
f0=set_funcs('sys_rhs',@(xx,p)0*xx(:,1,:),'sys_tau',@()funcs.sys_tau(),'x_vectorized',true);
L0=psol_jac_simple(f0,pt,[],'period',false);
L1=psol_jac_simple(f1,pt,[],'period',false);
Lid=L1-L0;