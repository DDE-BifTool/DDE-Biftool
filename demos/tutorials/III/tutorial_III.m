%%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: tutorial_III.m 134 2016-09-12 11:10:44Z mmbosschaert $
%%

clear;
close all;
addpath('../../../ddebiftool',...
    '../../../ddebiftool_extra_nmfm',...
    '../../../ddebiftool_utilities');
%% Right-hand side
funcs=set_symfuncs(@sym_HollingTanner_mf,'sys_tau',@()2);
%funcs=set_funcs(...
%    'sys_rhs', @HollingTanner_rhs,...
%    'sys_tau', @()2,...
%    'sys_deri', @HollingTanner_deri,...
%    'sys_mfderi',@HollingTanner_mfderi);
%% Sequence of parameters
indbeta=1;
indtau=2;
inda=3;
indm=4;
indh=5;
inddelta=6;
getpar=@(x,i)arrayfun(@(p)p.parameter(i),x.point);
getx=@(x,i)arrayfun(@(p)p.x(i),x.point);
bgetpar=@(x,i,bif)arrayfun(@(p)p.parameter(i),x.point(br_getflags(x,bif)));
bgetx=@(x,i,bif)arrayfun(@(p)p.x(i),x.point(br_getflags(x,bif)));
%% Continuation of steady state branch
fprintf('----- Steady-state branch -----\n');
beta=0.5;
a=0.5; 
m=(1/30)*(1-beta/(a*beta+1)); 
h=(1/4)*(beta/(a*beta+1)-1)^2+m*beta/(a*beta+1); 
tau=1/4*(a*beta+1)^2/beta; 
beta=0.4;
delta=0.5409;
stst.parameter=[beta,tau,a,m,h,delta];
parameter_bd={'max_bound',[indbeta,0.6; inddelta,0.7],...
    'min_bound',[indbeta,0.4;inddelta,0.4],...
    'max_step',[0,0.1;indbeta,5e-3;inddelta,5e-3]};
xster=-(1/2)*((beta/(a*beta+1)+2*m-1)+sqrt((1-2*m-beta/(a*beta+1))^2+4*(m*(1-m)-h)));
yster=beta*xster;

stst.x=[xster; yster];

% continue steady state in beta
contpar=indbeta;
stst_branch0 = SetupStst(funcs,'x',stst.x,'parameter',stst.parameter,...
    'contpar',contpar,'max_step',[0,0.005],'min_bound',...
    [contpar 0.4],'max_bound',[contpar 0.6],...
    'newheuristics_tests',0);
figure(1);clf
ax1=gca;

[stst_branch0] = br_contn(funcs,stst_branch0,100);
[stst_branch_wbifs,stst_testfuncs]=LocateSpecialPoints(funcs,stst_branch0);
nunst_stst=GetStability(stst_branch_wbifs);

%% Plot bifurcation diagram
beta_stst=getpar(stst_branch_wbifs,indbeta);
x1_stst=getx(stst_branch_wbifs,1);
figure(2);clf
ax2=gca;
cla(ax2);
plot(ax2,beta_stst(nunst_stst==0),x1_stst(nunst_stst==0),'g.',...
    beta_stst(nunst_stst==1),x1_stst(nunst_stst==1),'r.',...
    beta_stst(nunst_stst==2),x1_stst(nunst_stst==2),'b.');
hold;
plot(bgetpar(stst_branch_wbifs,indbeta,'hopf'),...
    bgetx(stst_branch_wbifs,1,'hopf'),'ks','MarkerSize',4)
plot(bgetpar(stst_branch_wbifs,indbeta,'fold'),...
    bgetx(stst_branch_wbifs,1,'fold'),'mo','MarkerSize',4)
    
stst_lgtext={'unstable=0','unstable=1','unstable=2','hopf','fold'};
legend(ax2,stst_lgtext,'location','west');
xlabel('$\beta$','Interpreter','LaTex');
ylabel('$x_1$','Interpreter','LaTex');

%% Continuation of Hopf bifurcation in two parameters
% after initialization of hopfbranch
fprintf('----- Hopf branch -----\n');
[hopf_branch0,suc] = SetupHopf(funcs, stst_branch_wbifs, br_getflags(stst_branch_wbifs,'hopf'),...
    'contpar', [inddelta,indbeta],...
    'dir', indbeta, 'step', 0.002,parameter_bd{:});
disp(['suc=',num2str(suc)]);
figure(3);clf
ax3=gca;
% title(ax3,'Hopf in beta-delta plane');
hopf_branch0=br_contn(funcs,hopf_branch0,300);

%% Continuation of fold bifurcation in two parameters
% after initialization of foldbranch
fprintf('----- fold branch -----\n');
[fold_branch0,suc]=SetupFold(funcs,stst_branch_wbifs,br_getflags(stst_branch_wbifs,'fold'),...
    'contpar', [inddelta,indbeta], 'dir', inddelta, 'step', 0.005, parameter_bd{:});
disp(['suc=',num2str(suc)]);
figure(3);
% title(ax2,'Fold in beta-delta plane');
fold_branch0=br_contn(funcs,fold_branch0,300);
fold_branch0 = br_rvers(fold_branch0);
fold_branch0=br_contn(funcs,fold_branch0,300);

%% Detect special points along Hopf curve
fprintf('----- Codimension-two detection along Hopf branch -----\n');
[hopf_branch_wbifs,hopftestfuncs]=LocateSpecialPoints(funcs,hopf_branch0);
nunst_hopf=GetStability(hopf_branch_wbifs,'exclude_trivial',true);

%% Detect special points along fold curve
fprintf('----- Codimension-two detection along fold branch -----\n');
[fold_branch_wbifs,foldtestfuncs]=LocateSpecialPoints(funcs,fold_branch0);

%% add Bogdanov-Takens bifurcation points the figure
beta_fold=getpar(fold_branch_wbifs,indbeta);
delta_fold=getpar(fold_branch_wbifs,inddelta);
beta_hopf=getpar(hopf_branch_wbifs,indbeta);
delta_hopf=getpar(hopf_branch_wbifs,inddelta);

delta_bt1=bgetpar(hopf_branch_wbifs,inddelta,'BT');
beta_bt1=bgetpar(hopf_branch_wbifs,indbeta,'BT');
delta_bt2=bgetpar(fold_branch_wbifs,inddelta,'BT');
beta_bt2=bgetpar(fold_branch_wbifs,indbeta,'BT');

figure(5); clf;
plot(beta_fold,delta_fold,'k','linewidth',1.2);
hold on;
plot(beta_hopf,delta_hopf,'b','linewidth',1.2);
plot(beta_bt1,delta_bt1,'r.','MarkerSize',5,'linewidth',1.2);
plot(beta_bt2,delta_bt2,'go','MarkerSize',5,'linewidth',1.2);

legend('fold','Hopf','BT from Hopf','BT from fold','Location','NorthWest')
xlabel('$\beta$','Interpreter','LaTex');
ylabel('$\delta$','Interpreter','LaTex');

axis([0.47   0.51   0.4   0.7])

%% Test functions for codimension-2 bifurcations along Hopf curve
figure(5);clf
plot(delta_hopf,[tanh(hopftestfuncs.genh(1,:));...
    hopftestfuncs.bt;...
    hopftestfuncs.zeho;...
    hopftestfuncs.hoho],'.');
grid on
legend({'tanh(genh)','BT','zero-Hopf','Hopf-Hopf'});
set(gca,'ylim',[-0.5,1.5]);
xlabel('$\delta$','Interpreter','LaTex')
ylabel('$\phi^H$','Interpreter','LaTex')

%% Test functions for codimension-2 bifurcations along fold curve
figure(6);clf
plot(delta_fold,[real(foldtestfuncs.cusp(1,:))/100;foldtestfuncs.bt;foldtestfuncs.zeho],'.');
grid on
legend({'cusp/100','BT','zero-Hopf'},'location','northwest');
set(gca,'ylim',[-1,1]);
xlabel('$\delta$','Interpreter','LaTex');
ylabel('$\phi^f$','Interpreter','LaTex');

%% save workspace to file
save('tutorial_III.mat')
