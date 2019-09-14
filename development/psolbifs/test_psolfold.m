%% test Psolfold
clear
base=[pwd(),'/../../'];
addpath([base,'ddebiftool/'],...
    [base,'ddebiftool_extra_psol/'],...
    [base,'ddebiftool_utilities/'],...
    [base,'ddebiftool_extra_nmfm']);
format compact
format short g
load('minimal_demo_parnames.mat')
load('minimal_demo_stst_psol_results.mat')
%%
ipf=find(nunst_per==0,1,'first');
pf0=per_orb.point(ipf);
pfoldini=dde_psolfold_init(funcs,pf0,[],per_orb.method);
method=dde_psolfold_method(per_orb.method);
method.point.print_residual_info=1;
[f,x0,data]=p_correc_setup(funcs,pfoldini,ind.b,method.point,'previous',pfoldini);
sm=setfield(method.point,'numderiv',true);
f0=p_correc_setup(funcs,pfoldini,ind.b,sm,'previous',pfoldini);
[r,J]=f(x0);
[r1,J1]=f0(x0);
method.point.hdev=0;
[pfold_res,suc]=p_correc(funcs,pfoldini,ind.b,[],method.point,0,pfoldini);
%%
figure(1);clf;
[pfbranch,suc]=SetupPsolfold(funcs,per_orb,ipf,'contpar',[ind.tau,ind.b],...
    'print_residual_info',1,'dir',ind.tau,'max_bound',[ind.b,0.61;ind.tau,20],...
    'max_step',[]);
pfbranch=br_contn(funcs,pfbranch,300);
