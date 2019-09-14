%% Tutorial V
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: tutorial_V.m 134 2016-09-12 11:10:44Z mmbosschaert $
%% Define system and load the workspace from the previous Tutorial
clear variables
close all
addpath(...
    '../../../ddebiftool',...
    '../../../ddebiftool_utilities',...
    '../../../ddebiftool_extra_psol',...
    '../../../ddebiftool_extra_nmfm',...
    '../IV');
plot_progress=true;
% load results from tutorial IV
load('../IV/neural.mat');
% need to define system again
g=@(x) (tanh(x-1)+tanh(1))*cosh(1)^2;
neuron_sys_rhs = @(xx,par) [
        -xx(1,1,:)-par(1)*g(par(2)*xx(1,2,:))+par(3)*g(par(4)*xx(2,3,:));
        -xx(2,1,:)-par(1)*g(par(2)*xx(2,2,:))+par(3)*g(par(4)*xx(1,3,:))];
tau_ind=[5 6];
funcs=set_funcs(...
    'sys_rhs', neuron_sys_rhs,...
    'sys_tau', @() tau_ind,...
    'sys_deri', @neural_deri,...
    'x_vectorized',true...
);
%% Verify that MATLAB is used
if dde_isoctave()
    error('This demo is not compatible with GNU Octave.')
end
%% continue fold point LPC1 in (c,a) (takes 62.994834 seconds on my comp)
ind_fold=ind_bif(2);
figure(1); clf;
% set max setsize for parameter c in per_orb
% then SetupPOfold and SetupPeriodDoubling will also use this value
per_orb.parameter.max_step =[indc 0.005];
[LPCfuncs,LPCbranch]=SetupPOfold(funcs,per_orb,ind_fold,...
    'contpar',[indc,inda],'dir',inda,'step',0.9e-3,'plot',plot_progress);
% set bounds
LPCbranch.parameter.max_bound=[[indc 1]; [inda 1]];
% contiue fold branch
LPCbranch=br_contn(LPCfuncs,LPCbranch,47);
LPCbranch=br_rvers(LPCbranch);
LPCbranch=br_contn(LPCfuncs,LPCbranch,150);
%% show structure of LPCfuncs, the extended system and auxiliary parameters
LPCfuncs
LPCfuncs.sys_rhs
LPCbranch.parameter
%% covert first point on the branch LPCbranch to generalized Hopf point
geho=p_tohopf(funcs,LPCfuncs.get_comp( LPCbranch.point(1),'solution'));
geho=nmfm_hopf(funcs,geho);
fprintf('First Lyapunov coefficient:%s\n',geho.nmfm.L1)
fprintf('Eigenvalues:\n')
method=df_mthod(funcs,'hopf');
geho.stability=p_stabil(funcs,geho,method.stability);
geho.stability.l1(1:6)'
geho=p_togenh(geho);
geho=nmfm_genh(funcs,geho);
fprintf('Second Lyapunov coefficient:%s\n',geho.nmfm.L2)
%% plot LPCbranch
figure(1); clf; hold on
getpars=@(br,ind) arrayfun(@(p)p.parameter(ind),br);
pl_LPCbranch=plot(getpars(LPCbranch.point,indc),...
    getpars(LPCbranch.point,inda));
plot(getpars(LPCbranch.point([1 23 68]),indc),...
    getpars(LPCbranch.point([1 23 68]),inda),'.r','MarkerSize',7)
text(0.51,0.24,'$GH$','Interpreter','LaTex');
text(0.41,0.385,'$CP_1$','Interpreter','LaTex');
text(0.43,-0.01,'$CP_2$','Interpreter','LaTex');
xlabel('$c$','Interpreter','LaTex');
ylabel('$a$','Interpreter','LaTex');
axis([0.4 1 -0.05 0.4])
%% remesh bifurcation points PD1 and PD3 on the branch per_orb for accuracy 
per_orb2=per_orb;
intervals=40;
degree=5;
% refine
indx=[ind_bif(1) ind_bif(4)];
per_orb2.point(indx)=arrayfun(@(p)p_remesh(p,degree,intervals),...
    per_orb.point(indx));
% and correct
per_orb2.point(indx)=arrayfun(@(p)p_correc(funcs,p,[],[],...
    per_orb.method.point),per_orb2.point(indx));
%% continue period doubling point 1 in (c,a) 
figure(2); clf
[pd1_funcs,pd1_br]=SetupPeriodDoubling(funcs,per_orb2,ind_bif(1),...
    'contpar',[indc,inda],'dir',inda,'step',1e-3,'plot',plot_progress);
pd1_br=br_contn(pd1_funcs,pd1_br,200);
xlabel('$c$','Interpreter','LaTex');
ylabel('$a$','Interpreter','LaTex');
%% continue period doubling point 3 in (c,a)
[pd3_funcs,pd3_br]=SetupPeriodDoubling(funcs,per_orb2,ind_bif(4),...
    'contpar',[indc,inda],'dir',inda,'step',1e-3,'plot',plot_progress);
pd3_br.parameter.max_bound=[[indc 1];[inda 1]];
pd3_br=br_rvers(pd3_br);
pd3_br=br_contn(pd3_funcs,pd3_br,400);
pd3_br=br_rvers(pd3_br);
pd3_br=br_contn(pd3_funcs,pd3_br,200);
%% continue NS1 point in (c,a)
ind_ns=find(abs(diff(nunst_per))==2);
[ns1_pfuncs,ns1_br]=SetupTorusBifurcation(funcs,per_orb,ind_ns,...
    'contpar',[indc,inda],'dir',inda,'step',1e-3,'plot',plot_progress);
ns1_br=br_contn(ns1_pfuncs,ns1_br,58);
ns1_br=br_rvers(ns1_br);
ns1_br=br_contn(ns1_pfuncs,ns1_br,30);
%% plot pd1_br and pd3_br
figure(1)
pl_pd1_br=plot(getpars(pd1_br.point,indc),getpars(pd1_br.point,inda));
pl_pd3_br=plot(getpars(pd3_br.point,indc),getpars(pd3_br.point,inda));
pl_ns1_br=plot(getpars(ns1_br.point,indc),getpars(ns1_br.point,inda));
pl_LPCs=plot(getpars([LPC1 LPC2],indc),getpars([LPC1 LPC2],inda),'ro');
pl_PDs=plot(getpars([PD1 PD2 PD3 PD4],indc),getpars([PD1 PD2 PD3 PD4],inda),'ks');
pl_NS1=plot(getpars(NS1,indc),getpars(NS1,inda),'dm');
legend([pl_LPCbranch pl_pd1_br pl_pd3_br pl_ns1_br pl_LPCs pl_PDs pl_NS1],...
    {'LPCbranch','pd1\_br','pd3\_br','pl\_ns1\_br','LPC','PD','NS_1'});
%% find points near CP2
points_near_CP2=find(getpars(LPCbranch.point,inda)<4*10^-4)
%% continue near CP2
figure(3); clf
pbranch_nearCP2=LPCbranch;
pbranch_nearCP2.point=LPCbranch.point(points_near_CP2(1));
pbranch_nearCP2.point(2)=LPCbranch.point(points_near_CP2(1));
pbranch_nearCP2.point(2).parameter(indc)=...
    pbranch_nearCP2.point(2).parameter(indc)+1e-06;
pbranch_nearCP2.method.point.minimal_accuracy=0.5e-6;
pbranch_nearCP2.parameter.max_step=[[indc 0.5e-06]; [inda 0.5e-06]];
pbranch_nearCP2.parameter.max_bound=[inda 3.0544e-04];
pbranch_nearCP2=br_rvers(pbranch_nearCP2);
pbranch_nearCP2=br_contn(LPCfuncs,pbranch_nearCP2,210);
% continue with relative larger stepsize
pbranch_nearCP2.parameter.max_step=[[indc 1e-06]; [inda 1e-06]];
pbranch_nearCP2.method.continuation.plot_progress=plot_progress;
pbranch_nearCP2=br_contn(LPCfuncs,pbranch_nearCP2,100);
% get stability
pbranch_nearCP2=br_stabl(LPCfuncs,pbranch_nearCP2,0,0);
%% plot
figure(4); clf; hold on
plot(getpars(LPCbranch.point,indc),getpars(LPCbranch.point,inda))
plot(getpars(pbranch_nearCP2.point,indc),...
    getpars(pbranch_nearCP2.point,inda),'-')
xlabel('$c$','Interpreter','LaTex');
ylabel('$a$','Interpreter','LaTex');
axis([0.459909117   0.459941577   0.000274128   0.000296046])
legend('LPCbranch','pbranch\_nearCP2')
%% plot Fold-Flip points
figure(5); clf; hold on
plot(getpars(pbranch_nearCP2.point,indc),...
    getpars(pbranch_nearCP2.point,inda),'-')
ind_FF2=50;
ind_FF3=54;
FF_pl=plot(getpars(pbranch_nearCP2.point([ind_FF2 ind_FF3]),indc),...
    getpars(pbranch_nearCP2.point([ind_FF2 ind_FF3]),inda),...
    '.r','MarkerSize',7);
FF2_text=text(0.459923,0.0002862,...
    '$FF_2$','Interpreter','LaTex');
FF3_text=text(0.45992,0.000285,...
    '$FF_3$','Interpreter','LaTex');
axis([0.459909117   0.459941577   0.000274128   0.000296046])
xlabel('$c$','Interpreter','LaTex');
ylabel('$a$','Interpreter','LaTex');
%% plot of stabillity of points 50,52,53 and 55 on pbranch_nearCP2
% demonstration two flip-fold bifurcations
figure(6)
subplot(2,2,1);
p_splot(pbranch_nearCP2.point(50))
axis([-1.0003219  -0.9997447  -0.0000316   0.0000316])
title('Point 50')
subplot(2,2,2);
p_splot(pbranch_nearCP2.point(52))
axis([-1.0003219  -0.9997447  -0.0000314   0.0000314])
title('Point 52')
subplot(2,2,3);
p_splot(pbranch_nearCP2.point(53))
axis([-1.0003217  -0.9997441  -0.0000314   0.0000314])
title('Point 53')
subplot(2,2,4);
p_splot(pbranch_nearCP2.point(55))
axis([-1.0003217  -0.9997441  -0.0000314   0.0000314])
title('Point 55')
%% continue limit points on the pbranch in c only
psolbr_from_LP=df_brnch(funcs,indc,'psol');
psol_from_LP=LPCfuncs.get_comp(pbranch_nearCP2.point(32),'solution');
% plot(psol_from_LP.parameter(indc),psol_from_LP.parameter(inda),'*')
% refine mash
method=df_mthod(funcs,'psol');
intervals=80;
degree=5;
psol_from_LP=p_remesh(psol_from_LP,degree,intervals)
[psol_from_LP,s]=p_correc(funcs,psol_from_LP,indc,[],method.point)
psolbr_from_LP.point=psol_from_LP;
% add second point
psolbr_from_LP.point(2)=psol_from_LP;
psolbr_from_LP.point(2).parameter(indc)=...
    psolbr_from_LP.point(2).parameter(indc)+1e-10;
% correct
[psolbr_from_LP.point(2),success]=p_correc(funcs,...
    psolbr_from_LP.point(2),[],[],method.point)
% continuation settings
psolbr_from_LP.method.continuation.steplength_growth_factor=1;
psolbr_from_LP.method.continuation.plot_progress=plot_progress;
% continue branch of period orbits
figure(7); clf
psolbr_from_LP=br_contn(funcs,psolbr_from_LP,50);
%% calculate stability
psolbr_from_LP=br_stabl(funcs,psolbr_from_LP,0,0);
[nunst_per,floqper,triv_defect,persols]=GetStability(psolbr_from_LP,...
    'exclude_trivial',true);
nunst_per'
%% add to bifurcation diagram near CP2  and plot stability
figure(8); clf; hold on;
box on
plot(getpars(pbranch_nearCP2.point,indc),...
    getpars(pbranch_nearCP2.point,inda),'-');
% period doubling branches
plot(getpars(pd1_br.point,indc),getpars(pd1_br.point,inda));
plot(getpars(pd3_br.point,indc),getpars(pd3_br.point,inda));
xlabel('$c$','Interpreter','LaTex');
ylabel('$a$','Interpreter','LaTex')
plot(getpars(psolbr_from_LP.point,indc),...
    getpars(psolbr_from_LP.point,inda),'*');
legend('pbranch\_nearCP2','pd1\_br','pd3\_br','psolbr\_from\_LP',...
    'Location','NorthWest');
axis([0.459935 0.459938 0.000293 0.000294])
% plot stability
figure(9); clf
br_splot2(psolbr_from_LP,nunst_per,indc,'amplitude',1.0);
xlabel('$c$','Interpreter','LaTex');
ylabel('$a$','Interpreter','LaTex');
% axis([0.4599357   0.4599367 1.6439   1.6446])
legend('Location','NorthWest')
drawnow;
%% continue NS
figure(10); clf
ind_ns=find(abs(diff(nunst_per))==2);
[ns_pfuncs,ns_br]=SetupTorusBifurcation(funcs,psolbr_from_LP,ind_ns,...
'contpar',[indc,inda],'dir',indc,'step',1e-2,'plot',plot_progress);
ns_br.method.continuation.steplength_growth_factor=1.2;
ns_br=br_contn(ns_pfuncs,ns_br,25);
ns_br=br_rvers(ns_br);
ns_br=br_contn(ns_pfuncs,ns_br,15);
%% plot NS-curve
figure(11); clf; hold on
plot(getpars(pbranch_nearCP2.point,indc),...
    getpars(pbranch_nearCP2.point,inda),'-');
% period doubling branches
plot(getpars(pd1_br.point,indc),getpars(pd1_br.point,inda));
plot(getpars(pd3_br.point,indc),getpars(pd3_br.point,inda));
plot(getpars(ns_br.point,indc),getpars(ns_br.point,inda))
% add bifurcation points
plot(ns_br.point(1).parameter(indc),ns_br.point(1).parameter(inda),'.k')
plot(ns_br.point(end).parameter(indc),...
    ns_br.point(end).parameter(inda),'.k')
% remove FF points for legend in tikz file
FF_pl.delete()
axis([0.459935028521171   0.459937974745440   ...
    0.000292965887357   0.000295220037358])
legend('pbranch\_nearCP2','pd2\_br','pd3\_br',...
    'ns\_br','Location','NorthWest')
%% save workspace to file
save('tutorial_V.mat')
