%%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: tutorial_II.m 134 2016-09-12 11:10:44Z mmbosschaert $
%%
clear;
close all;

addpath('../../../ddebiftool',...
    '../../../ddebiftool_utilities');

g = @(z) (tanh(z-1)+tanh(1))*cosh(1)^2;

neuron_sys_rhs = @(xx,par) [
    -xx(1,1)-par(1)*g(par(2)*xx(1,2))+par(3)*g(par(4)*xx(2,3));
    -xx(2,1)-par(1)*g(par(2)*xx(2,2))+par(3)*g(par(4)*xx(1,3))
];
    
% par = [a, b, c, d, tau1, tau2];
tau_ind=[5 6];
funcs=set_funcs(...
    'sys_rhs', neuron_sys_rhs,...
    'sys_tau', @() tau_ind,...
    'sys_deri', @neural_deri...
);

%% initial guess
% set parameter values
a=0.25;
b=2.0;
c=15/29;
d=1.2;
tau1=12.7;
tau2=20.2;

% define point structure for a single steady state
stst.kind='stst';
stst.parameter=[a, b, c, d, tau1, tau2];
stst.x=[0;0]

%% Linear stability of initial equilibrium
method=df_mthod(funcs,'stst');
method.stability.minimal_real_part=-0.04;
[stst,success]=p_correc(funcs,stst,[],[],method.point)
% compute its stability:
stst.stability=p_stabil(funcs,stst,method.stability)
figure(1); clf;
p_splot(stst); % plot its stability:

%% initialize branch
% get an empty branch with ind_a as a free parameter:
ind_a=1;
branch1=df_brnch(funcs,ind_a,'stst')
branch1.parameter
branch1.parameter.min_bound
% set bounds for continuation parameter
branch1.parameter.min_bound(1,:)=[ind_a 0];
branch1.parameter.max_bound(1,:)=[ind_a .4];
branch1.parameter.max_step(1,:)=[ind_a 0.05];
% use stst as a first branch point:
branch1.point=stst;

%% Extend and continue branch of trivial equilibria
stst2=stst;
stst2.parameter(ind_a)=stst2.parameter(ind_a)-0.001;
[stst2,success]=p_correc(funcs,stst2,[],[],method.point)
% use as a second branch point:
branch1.point(2)=stst2;
branch1.method.continuation.plot=1;

figure(2); clf
% continue in one direction:
[branch1,s,f,r]=br_contn(funcs,branch1,300)
% turn the branch around:
branch1=br_rvers(branch1);
% continue in the other direction:
[branch1,s,f,r]=br_contn(funcs,branch1,300)

%% Stability of branch of equilibria
branch1.method.stability.minimal_real_part=-0.01;
branch1=br_stabl(funcs,branch1,0,0);

% obtain suitable scalar measures to plot stability along branch:
[xm,ym]=df_measr(1,branch1)
figure(3); clf;
br_plot(branch1,xm,ym,'b'); % plot stability along branch:
% ym.subfield='l0';
br_plot(branch1,xm,ym,'c');
plot([0 0.4],[0 0],'-.');
axis([0 0.4 -0.01 0.01]);
xlabel('$a$','Interpreter','LaTex');
ylabel('$\Re(\lambda)$','Interpreter','LaTex');
% plot stability versus point number:
figure(4); clf;
br_plot(branch1,[],ym,'b.');
br_plot(branch1,[],ym,'b-');
plot([5 45],[0 0],'-.');
% change markersize
hline = findobj(gcf, 'type', 'line');
set(hline,'MarkerSize',5,'linewidth',0.8);
xlabel('point number along branch');
ylabel('$\Re(\lambda)$','Interpreter','LaTex');
axis([5  45  -0.01   0.01]);
%% Locating the first Hopf point
ind_hopf=find(arrayfun(@(x)~isempty(x.stability.l0) ...
    && real(x.stability.l0(1))>0,branch1.point),1,'first');
hopf=p_tohopf(funcs,branch1.point(ind_hopf));

% get hopf calculation method parameters:
method=df_mthod(funcs,'hopf');
method.stability.minimal_real_part=-0.02;

% correct hopf
[hopf,success]=p_correc(funcs,hopf,ind_a,[],method.point)

% compute stability of hopf point
hopf.stability=p_stabil(funcs,hopf,method.stability);

% plot stability
figure(5); clf;
p_splot(hopf);

%% Compute the first Lyapunov coefficient
addpath('../../../ddebiftool_extra_nmfm') % load the normal form package
hopf=nmfm_hopf(funcs,hopf);
fprintf('First Lyapunov coefficient l1: %g\n',hopf.nmfm.L1);

%% plot dynamics near Hopf
figure(6)

% just before the Hopf bifurcation
history=[0.1 0];
par=hopf.parameter;
par(1)=par(1)-0.04;
tfinal=1800;
sol1 = dde23(@(t,y,Z) funcs.sys_rhs([y,Z],par),par(funcs.sys_tau()),...
    history,[0 tfinal]);

plot(sol1.x,sol1.y(2,:))

xlabel('$t$','Interpreter','LaTex');
ylabel('$x$','Interpreter','LaTex');

par=hopf.parameter;
par(1)=par(1)+0.04;
sol1 = dde23(@(t,y,Z) funcs.sys_rhs([y,Z],par),par(funcs.sys_tau()),...
    history,[0 tfinal]);
history=[1 0];
sol2 = dde23(@(t,y,Z) funcs.sys_rhs([y,Z],par),par(funcs.sys_tau()),...
    history,[0 tfinal]);

figure(7)
hold
plot(sol1.x,sol1.y(2,:))
plot(sol2.x,sol2.y(2,:))

xlabel('$t$','Interpreter','LaTex');
ylabel('$x$','Interpreter','LaTex');
