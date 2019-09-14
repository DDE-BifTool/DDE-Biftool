clear
addpath('../../ddebiftool',...
    '../../ddebiftool_extra_psol',...
    '../../ddebiftool_extra_nmfm/',...
    '../../ddebiftool_utilities',...
    '../tools');
format compact
format short g
%% test some basic examples for fold and Hopf bifurcations
ntau=2;
parnames={'p','tau1','tau2'};
cip=[parnames;num2cell(1:length(parnames))];
ip=struct(cip{:});
funcs=set_symfuncs(@sym_basic_fold,'sys_tau',[ip.tau1,ip.tau2]);
x0=2;
par0([ip.p,    ip.tau1,ip.tau2])=...
    [  x0^3-x0,  2,      1];
[br0,suc]=SetupStst(funcs,'contpar',ip.p,'step',-0.01,'x',x0*[1;1],'parameter',par0,...
    'print_residual_info',1,'max_bound',[ip.p,10],'max_step',[0,0.2]);
disp(suc)
figure(1);clf;
br1=br_contn(funcs,br0,100);
%%
[br1_wbifs,br1_tests,br1_bifs,br1_types]=LocateSpecialPoints(funcs,br1,'debug',true,'newton_max_iterations',6);
%%
[brf0,suc]=SetupFold(funcs,br1_wbifs,br1_bifs(2),'contpar',[ip.tau1,ip.p],'dir',ip.tau1,...
    'print_residual_info',1,'max_bound',[ip.p,10;ip.tau1,3],'step',-0.01);
clf;
brf=br_contn(funcs,brf0,100);
brf=br_rvers(brf);
brf=br_contn(funcs,brf,100);
[brf_wbifs,f_tests,f_indices,f_types]=LocateSpecialPoints(funcs,brf,'debug',true);
%%
[brh0,suc]=SetupHopf(funcs,br1_wbifs,br1_bifs(1),'contpar',[ip.tau1,ip.p],'dir',ip.tau1,...
    'print_residual_info',1,'max_bound',[ip.p,10;ip.tau1,3],'step',-0.01);
brh=br_contn(funcs,brh0,100);
brh=br_rvers(brh);
brh=br_contn(funcs,brh,100);
[brh_wbifs,h_tests,h_indices,h_types]=LocateSpecialPoints(funcs,brh,'debug',true);
