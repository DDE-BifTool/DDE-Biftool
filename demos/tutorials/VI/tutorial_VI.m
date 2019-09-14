%%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: tutorial_VI.m 134 2016-09-12 11:10:44Z mmbosschaert $
%%
%% Tutorial VI
%% Define system
clear variables
close all
plot_progress=false;
addpath(...
    '../../../ddebiftool',...
    '../../../ddebiftool_utilities',...
    '../../../ddebiftool_extra_psol',...
    '../../../ddebiftool_extra_nmfm',...
    '../IV');
g=@(x) (tanh(x-1)+tanh(1))*cosh(1)^2;
neuron_sys_rhs = @(xx,par) [
        -xx(1,1,:)-par(1)*g(par(2)*xx(1,2,:))+par(3)*g(par(4)*xx(2,3,:));
        -xx(2,1,:)-par(1)*g(par(2)*xx(2,2,:))+par(3)*g(par(4)*xx(1,3,:))];
tau_ind=[5 6];
funcs=set_funcs(...
    'sys_rhs', neuron_sys_rhs,...
    'sys_tau', @() tau_ind,...
    'x_vectorized',true...
);
%% Set parameter names
% This could also be loaded from mat file |'FHN_parnames.mat'|.
parnames={'a','b','c','d','tau1','tau2'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
%% Set funcs structure
funcs=set_symfuncs(@sym_neuron_mf,'sys_tau',@() [ind.tau1, ind.tau2]);
%% continue steady state
x0=[0; 0];
% set parameters
a=0.069;
b=2;
c=0.4;
d=1.2;
tau1=12.99;
tau2=20.15;
par=[a b c d tau1 tau2];
% set continuation parameter
contpar=ind.c;
% setup steady-state branch
stst_branch=SetupStst(funcs,'x',x0,'parameter',par,'step',0.01,...
    'contpar',contpar,'max_step',[contpar,0.01],...
    'max_bound',[[ind.c,1]; [ind.a, 1]],...
    'min_bound',[[ind.c,-0.1]; [ind.a, -0.1]],...
    'plot',plot_progress);
disp('--- Trivial equilibria ---');
figure(1);clf
stst_branch=br_contn(funcs,stst_branch,1000);
% calculate stability of the steady-states
stst_branch=br_stabl(funcs,stst_branch,0,0);
nunst=GetStability(stst_branch);
% detect bifurcation points (Hopf and fold point) along the stst curve
[stst_branch,stst_testfuncs]=LocateSpecialPoints(funcs,stst_branch);
%% continue Hopf point in (a,c)
figure(2); clf;
hopf_ind=br_getflags(stst_branch,'hopf');
fprintf('----- Hopf branch 1 -----\n');
[hopf_branch1,suc] = SetupHopf(funcs, stst_branch, hopf_ind(1),...
    'contpar', [ind.c,ind.a],'dir', ind.c, 'step', 0.002,'plot',plot_progress);
figure(2);clf; hold on
hopf_branch1=br_contn(funcs,hopf_branch1,300);
hopf_branch1=br_rvers(hopf_branch1);
hopf_branch1=br_contn(funcs,hopf_branch1,300);
fprintf('----- Hopf branch 2 -----\n');
[hopf_branch2,suc] = SetupHopf(funcs, stst_branch, hopf_ind(2),...
    'contpar', [ind.c,ind.a],'dir', ind.c, 'step', 0.002,'plot',plot_progress);
hopf_branch2=br_contn(funcs,hopf_branch2,300);
hopf_branch2=br_rvers(hopf_branch2);
hopf_branch2=br_contn(funcs,hopf_branch2,300);
%% plot hopf branches
figure(2);clf; hold on
getpars=@(points,ind) arrayfun(@(p)p.parameter(ind),points);
% plot hopf branches
cm=colormap('lines');
plot(getpars(hopf_branch1.point,ind.c),getpars(hopf_branch1.point,ind.a));
plot(getpars(hopf_branch2.point,ind.c),getpars(hopf_branch2.point,ind.a));
% plot hopf points $H_1$ and $H_2$
H1=stst_branch.point(hopf_ind(1));
H2=stst_branch.point(hopf_ind(2));
plot(H1.parameter(ind.c),H1.parameter(ind.a),'s','color',cm(1,:))
plot(H2.parameter(ind.c),H2.parameter(ind.a),'s','color',cm(2,:))
xlabel('$c$','Interpreter','LaTex');
ylabel('$a$','Interpreter','LaTex');
legend({'hopf\_branch1','hopf\_branch2','$H_1$','$H_2$'},...
    'Interpreter','LaTex');
%% detect codim-2 points
fprintf('----- Codimension-two detection along the first Hopf branch -----\n');
[hopf_branch1,hopf1_testfuncs]=LocateSpecialPoints(funcs,hopf_branch1);
fprintf('----- Codimension-two detection along the second Hopf branch -----\n');
[hopf_branch2,hopf2_testfuncs]=LocateSpecialPoints(funcs,hopf_branch2);
%% plot hopf branches with codim-2 points
hoho_ind_br1=br_getflags(hopf_branch1,'hoho');
genh_ind_br1=br_getflags(hopf_branch1,'genh');
zeho_ind_br1=br_getflags(hopf_branch1,'zeho');
hoho_ind_br2=br_getflags(hopf_branch2,'hoho');
genh_ind_br2=br_getflags(hopf_branch2,'genh');
zeho_ind_br2=br_getflags(hopf_branch2,'zeho');
figure(3); clf; hold on
cm=colormap('lines');
pl_hopf_branch1 = plot(getpars(hopf_branch1.point,ind.c),...
    getpars(hopf_branch1.point,ind.a),'color',cm(1,:));
plot(getpars(hopf_branch2.point,ind.c),...
    getpars(hopf_branch2.point,ind.a),'color',cm(1,:));
pl_hoho=plot(getpars(hopf_branch1.point(hoho_ind_br1),ind.c),...
    getpars(hopf_branch1.point(hoho_ind_br1),ind.a),'ko','MarkerSize',6);
plot(getpars(hopf_branch2.point(hoho_ind_br2),ind.c),...
    getpars(hopf_branch2.point(hoho_ind_br2),ind.a),'ko','MarkerSize',6);
pl_genh=plot(getpars(hopf_branch1.point(genh_ind_br1),ind.c),...
    getpars(hopf_branch1.point(genh_ind_br1),ind.a),'bd','MarkerSize',6);
plot(getpars(hopf_branch2.point(genh_ind_br1),ind.c),...
    getpars(hopf_branch2.point(genh_ind_br1),ind.a),'bd','MarkerSize',6);
plot(getpars(hopf_branch1.point(zeho_ind_br1),ind.c),...
    getpars(hopf_branch1.point(zeho_ind_br1),ind.a),'rs','MarkerSize',6);
pl_zeho=plot(getpars(hopf_branch2.point(zeho_ind_br2),ind.c),...
    getpars(hopf_branch2.point(zeho_ind_br2),ind.a),'rs','MarkerSize',6);
xlabel('$c$','Interpreter','LaTex');
ylabel('$a$','Interpreter','LaTex');
legend([pl_hopf_branch1 pl_hoho pl_genh pl_zeho],...
    {'Hopf branhes','Hopf-Hopf','generalized Hopf','zero-Hopf'},...
    'Interpreter','LaTex');
HH1=hopf_branch1.point(hoho_ind_br1(2));
hh1_text=text(HH1.parameter(ind.c)+0.01,HH1.parameter(ind.a)+0.01,'$HH_1$',...
    'Interpreter','LaTex');
HH2=hopf_branch1.point(hoho_ind_br1(3));
hh2_text=text(HH2.parameter(ind.c)+0.01,HH2.parameter(ind.a)+0.01,'$HH_2$',...
    'Interpreter','LaTex');
ax1=gca;
%% Constructing an initial small-amplitude orbit near HH1 on hopf_branch1
par(ind.a)=0.305045919228062;
par(ind.c)=0.543129788993018;
stst.x=[0; 0];
stst.parameter=par;
stst.kind='stst';
stst.stability=p_stabil(funcs,stst,stst_branch.method.stability)
hopf_nearHH=p_tohopf(funcs,stst);
method=df_mthod(funcs,'hopf');
[hopf_nearHH,s]=p_correc(funcs,hopf_nearHH,ind.c,[],method.point)
figure(4); clf; maindiagram;
pl_hopf_nearHH=plot(hopf_nearHH.parameter(ind.c),...
    hopf_nearHH.parameter(ind.a),'*','color',cm(2,:));
axis([0.505671175858481   0.609740066682240...
    0.272963646099839   0.309983935265068])
set(hh1_text,'Position',[0.575003564583368 0.280259687296467 0]);
legend([pl_hopf_branch1 pl_hopf_nearHH],...
    {'Hopf branhes','hopf\_nearHH'},'Interpreter','LaTex');
%% continue periodic orbits emanating from hopf_nearHH
% calculate the First Lyapunov coefficient
hopf_nearHH=nmfm_hopf(funcs,hopf_nearHH);
fprintf('First Lyapunov coefficient: %f\n',hopf_nearHH.nmfm.L1)
% plot stability
hopf_nearHH.stability=p_stabil(funcs,hopf_nearHH,stst_branch.method.stability);
figure(5); p_splot(hopf_nearHH);
axis([-0.0035    0.0032   -1.0979    1.3353])
% consturct small periodic orbit near hopf_nearHH
intervals=20;
degree=4;
[psol_nearHH,stepcond]=p_topsol(funcs,hopf_nearHH,0.01,degree,intervals);
% correct periodic solution guess
method=df_mthod(funcs,'psol');
[psol_nearHH,success]=p_correc(funcs,psol_nearHH,ind.c,stepcond,method.point);
% plot stability
psol_nearHH.stability=p_stabil(funcs,psol_nearHH,method.stability);
figure(6); clf
p_splot(psol_nearHH)
% construct second orbit
[psol2_nearHH,stepcond]=p_topsol(funcs,hopf_nearHH,0.02,degree,intervals);
[psol2_nearHH,success]=p_correc(funcs,psol2_nearHH,ind.c,stepcond,method.point)
psol2_nearHH.stability=p_stabil(funcs,psol_nearHH,method.stability);
% manual construction and continuation of per_orb_nearHH
figure(7); clf
per_orb_nearHH=df_brnch(funcs,ind.c,'psol'); % empty branch:
per_orb_nearHH.method.continuation.plot_progress=plot_progress;
per_orb_nearHH.point=psol_nearHH;
per_orb_nearHH.point(2)=psol2_nearHH;
per_orb_nearHH.method.continuation.steplength_growth_factor=1.05;
per_orb_nearHH=br_rvers(per_orb_nearHH);
% compute periodic solutions branch
per_orb_nearHH=br_contn(funcs,per_orb_nearHH,58);
xlabel('c');ylabel('amplitude');
% plot stability
per_orb_nearHH=br_stabl(funcs,per_orb_nearHH,0,0);
figure(7);
[nunst_per,floqper,triv_defect,persols]=...
    GetStability(per_orb_nearHH,'exclude_trivial',true);
clf;br_splot2(per_orb_nearHH,nunst_per,ind.c,'');
% hLeg=findobj(gcf,'tag','legend');
% set(hLeg,'Location','NorthWest')
xlabel('$c$','Interpreter','LaTex');
ylabel('amplitude','Interpreter','LaTex');
% plot per_orb_nearHH in bifurcation diagram
figure(8); clf; maindiagram;
plot(getpars(per_orb_nearHH.point,ind.c),...
    getpars(per_orb_nearHH.point,ind.a),'*','color',cm(2,:))
legend('off')
%% plot periodic orbits
ppars=arrayfun(@(x)x.parameter(ind.c),per_orb_nearHH.point);
pmeshes=cell2mat(arrayfun(@(x)x.mesh(:),per_orb_nearHH.point,'uniformoutput',false));
pprofs1=cell2mat(arrayfun(@(x)x.profile(1,:)',per_orb_nearHH.point,'uniformoutput',false));
pprofs2=cell2mat(arrayfun(@(x)x.profile(2,:)',per_orb_nearHH.point,'uniformoutput',false));
figure(9);clf
hold on
for i=1:length(ppars)
    plot3(repmat(ppars(i),length(pmeshes(:,i)),1)',pprofs1(:,i),pprofs2(:,i),'color',cm(nunst_per(i)+1,:));
end
xlabel('$c$','Interpreter','LaTex');
ylabel('$x_1$','Interpreter','LaTex');
zlabel('$x_2$','Interpreter','LaTex');
view(153,20); drawnow
%% Constructing an initial small-amplitude orbit near HH1 on hopf_branch2
par(ind.a)=0.266488849203010;
par(ind.c)=0.591265246332908;
stst.x=[0; 0];
stst.parameter=par;
stst.kind='stst';
stst.stability=p_stabil(funcs,stst,stst_branch.method.stability)
hopf2_nearHH=p_tohopf(funcs,stst);
method=df_mthod(funcs,'hopf');
[hopf2_nearHH,s]=p_correc(funcs,hopf2_nearHH,ind.c,[],method.point)
figure(10); clf; maindiagram;
plot(hopf2_nearHH.parameter(ind.c),hopf2_nearHH.parameter(ind.a),'*r',...
    'HandleVisibility','off')
axis([0.542325356612850   0.598925883841189   0.260032176134403   0.304877020640333])
hh_text=findobj(gca,'Type','Text');
set(hh_text(2),'Position',[0.575003564583368 0.280259687296467 0])
%% continue periodic orbits emanating from hopf_nearHH
% calculate the First Lyapunov coefficient
hopf2_nearHH=nmfm_hopf(funcs,hopf2_nearHH);
fprintf('First Lyapunov coefficient: %f\n',hopf2_nearHH.nmfm.L1)
% plot stability
hopf2_nearHH.stability=p_stabil(funcs,hopf2_nearHH,stst_branch.method.stability);
figure(5); clf; p_splot(hopf2_nearHH);
axis([-0.004   0.004  -0.6   0.6])
% consturct small periodic orbit near hopf_nearHH
intervals=20;
degree=4;
[psol2_nearHH,stepcond]=p_topsol(funcs,hopf2_nearHH,0.01,degree,intervals);
% correct periodic solution guess
method=df_mthod(funcs,'psol');
[psol2_nearHH,success]=p_correc(funcs,psol2_nearHH,ind.c,stepcond,method.point)
% plot stability
psol2_nearHH.stability=p_stabil(funcs,psol2_nearHH,method.stability);
figure(11); clf
p_splot(psol2_nearHH)
% construct second orbit
[psol3_nearHH,stepcond]=p_topsol(funcs,hopf2_nearHH,0.02,degree,intervals);
[psol3_nearHH,success]=p_correc(funcs,psol3_nearHH,ind.c,stepcond,method.point)
psol3_nearHH.stability=p_stabil(funcs,psol3_nearHH,method.stability);
% manual construction and continuation of per_orb_nearHH
figure(12); clf
per_orb2_nearHH=df_brnch(funcs,ind.c,'psol'); % empty branch:
per_orb2_nearHH.method.continuation.plot_progress=plot_progress;
per_orb2_nearHH.point=psol2_nearHH;
per_orb2_nearHH.point(2)=psol3_nearHH;
per_orb2_nearHH.method.continuation.steplength_growth_factor=1.05;
per_orb2_nearHH=br_rvers(per_orb2_nearHH);
per_orb2_nearHH=br_contn(funcs,per_orb2_nearHH,26); % compute periodic solutions branch
xlabel('c');ylabel('amplitude');
% plot stability
per_orb2_nearHH=br_stabl(funcs,per_orb2_nearHH,0,0);
figure(12); clf
nunst_per2=GetStability(per_orb2_nearHH,'exclude_trivial',true);
clf;br_splot2(per_orb2_nearHH,nunst_per2,ind.c,'');
hLeg=findobj(gcf,'tag','legend');
set(hLeg,'Location','NorthWest')
xlabel('$c$','Interpreter','LaTex');
ylabel('amplitude');
drawnow
%% plot per_orb_nearHH in bifurcation diagram
figure(13); clf; maindiagram
plot(getpars(per_orb2_nearHH.point,ind.c),getpars(per_orb2_nearHH.point,ind.a),'*r')
axis([0.542325356612850   0.598925883841189   0.260032176134403   0.304877020640333])
hh_text=findobj(gca,'Type','Text');
set(hh_text(2),'Position',[0.575003564583368 0.280259687296467 0])
drawnow
%% plot periodic orbits
cm=colormap('lines');
figure(14); clf
ppars=arrayfun(@(x)x.parameter(ind.c),per_orb2_nearHH.point);
pmeshes=cell2mat(arrayfun(@(x)x.mesh(:),per_orb2_nearHH.point,'uniformoutput',false));
pprofs1=cell2mat(arrayfun(@(x)x.profile(1,:)',per_orb2_nearHH.point,'uniformoutput',false));
pprofs2=cell2mat(arrayfun(@(x)x.profile(2,:)',per_orb2_nearHH.point,'uniformoutput',false));
hold on
for i=1:length(ppars)
    plot3(repmat(ppars(i),length(pmeshes(:,i)),1)',pprofs1(:,i),pprofs2(:,i),'color',cm(nunst_per2(i)+1,:));
end
xlabel('$c$','Interpreter','LaTex');
ylabel('$x_1$','Interpreter','LaTex');
zlabel('$x_2$','Interpreter','LaTex');
view(153,20); drawnow;
%% continue NS1 point in (c,a)
ind_ns=find(abs(diff(nunst_per))==2);
figure(15); clf; maindiagram;
[pfuncs_ns1,ns1_br]=SetupTorusBifurcation(funcs,per_orb_nearHH,ind_ns(1),...
'contpar',[ind.c,ind.a],'dir',ind.c,'step',1e-2,'plot',plot_progress);
ns1_br.method.continuation.steplength_growth_factor=1.2;
ns1_br=br_contn(pfuncs_ns1,ns1_br,20);
ns1_br=br_rvers(ns1_br);
ns1_br=br_contn(pfuncs_ns1,ns1_br,45);
drawnow
%% continue NS2 point in (c,a)
ind_ns=find(abs(diff(nunst_per2))==2);
[pfuncs_ns2,ns2_br]=SetupTorusBifurcation(funcs,per_orb2_nearHH,ind_ns(1),...
'contpar',[ind.c,ind.a],'dir',ind.c,'step',1e-3,'plot',plot_progress);
ns2_br=br_rvers(ns2_br);
ns2_br=br_contn(pfuncs_ns2, ns2_br,50);
ns2_br=br_rvers(ns2_br);
ns2_br=br_contn(pfuncs_ns2, ns2_br,100);
%% plot ns curves with focus on ns1
figure(16); clf; maindiagram;
plot(getpars(ns1_br.point,ind.c),getpars(ns1_br.point,ind.a),'linewidth',1);
plot(getpars(ns2_br.point,ind.c),getpars(ns2_br.point,ind.a),'linewidth',1);
axis([-0.0122    0.6388    0.2436    0.5321]); drawnow
%% plot ns curves with focus on ns2
figure(17); clf; maindiagram;
plot(getpars(ns1_br.point,ind.c),getpars(ns1_br.point,ind.a),'linewidth',1);
plot(getpars(ns2_br.point,ind.c),getpars(ns2_br.point,ind.a),'linewidth',1);
axis([0.4887    0.5976    0.0851    0.4045]); drawnow
%% unstable torus (this takes almost 5 hours!)
% tic
% tfinal=3*30000;
% par(ind.a)=0.5598/b;
% par(ind.c)=0.6890/d;
% % par(ind.c)=0.573979725801627; % stable symetric periodic orbit
% % par(ind.c)=0.573985422125899; % stable asymetric periodic orbit
% % par(ind.c)=0.573993061564200; % stable symetric periodic orbit
% history1 = [0.05 0.075];
% options = ddeset('RelTol',1e-06,'AbsTol',1e-08);
% sol1 = dde23(@(t,y,Z) neuron_sys_rhs([y,Z],par),[tau1 tau2],history1,[0 tfinal],options)
% toc
% %% plot torus 
% t=linspace(30000-2500,30000,2500);
% y=deval(sol1,t);
% % plot simulation
% % figure(5); clf
% % plot(sol1.x,sol1.y(2,:))
% % plot simulation
% figure(12); clf;
% plot(y(1,:),y(2,:));
% xlabel('$x_1$','Interpreter','LaTex');
% ylabel('$x_2$','Interpreter','LaTex');
% % add second torus to plot
% % figure; clf;
% % plot(y(1,:),y(2,:));
% t=linspace(tfinal-2500,tfinal,2500);
% y=deval(sol1,t);
% hold on
% plot(y(1,:),y(2,:));
% % toc
% % figure(3)
% % plot(par(ind.c),par(ind.a),'*')
%% continue LPC1
ind_fold=find(abs(diff(nunst_per))==1);
[lpc1funcs,lpc1_br]=SetupPOfold(funcs,per_orb_nearHH,ind_fold,...
    'contpar',[ind.c,ind.a],'dir',ind.a,'step',1e-3,'plot',plot_progress);
figure(18); clf; maindiagram;
lpc1_br.parameter.max_bound=[ind.c 1];
lpc1_br=br_contn(lpc1funcs,lpc1_br,200);
lpc1_br=br_rvers(lpc1_br);
lpc1_br=br_contn(lpc1funcs,lpc1_br,90);
%% continue LPC2
ind_fold=find(abs(diff(nunst_per2))==1);
[lpc2funcs,lpc2_br]=SetupPOfold(funcs,per_orb2_nearHH,ind_fold,...
    'contpar',[ind.c,ind.a],'dir',ind.a,'step',1e-3,'plot',plot_progress);
%% plot lpc curves
figure(19); clf; maindiagram;
lpc2_br.parameter.max_bound=[ind.c 1];
lpc2_br=br_contn(lpc2funcs,lpc2_br,150);
lpc2_br=br_rvers(lpc2_br);
lpc2_br=br_contn(lpc2funcs,lpc2_br,55);
figure(20); clf; maindiagram;
plot(getpars(ns1_br.point,ind.c),getpars(ns1_br.point,ind.a));
plot(getpars(ns2_br.point,ind.c),getpars(ns2_br.point,ind.a));
plot(getpars(lpc1_br.point,ind.c),getpars(lpc1_br.point,ind.a),'color',cm(5,:));
plot(getpars(lpc2_br.point,ind.c),getpars(lpc2_br.point,ind.a),'color',cm(7,:));
%% save
save('tutorial_VI.mat')
