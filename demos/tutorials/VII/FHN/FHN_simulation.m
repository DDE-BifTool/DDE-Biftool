                                      %% Simulation with MATLAB dde23
%% Clean workspace and add DDE-BifTool scripts to the MATLAB search path
clear;      % clear variables
close all;	% close figures
ddebiftoolpath='../../../../';
addpath(strcat(ddebiftoolpath,'ddebiftool'),...
        strcat(ddebiftoolpath,'ddebiftool_extra_psol'),...
        strcat(ddebiftoolpath,'ddebiftool_extra_nmfm'),...
        strcat(ddebiftoolpath,'ddebiftool_utilities'));
load('FHN_results.mat')
%% Point in region I
beta0=1.8779;
alpha0=-1.1001;
% Point near the steady-state
x1=0;
x2=0.01;
% Integrate
tfinal=1000;
sol = dde23(@(t,y,Z) funcs.sys_rhs([y,Z],[beta0 alpha0 tau0]),...
    tau0,[x1 x2],[0 tfinal]);
t=linspace(0,tfinal,1000);
y=deval(sol,t);
% Plot
figure(1);clf;
title('Point in region I')
xlabel('$u_1$','Interpreter','LaTex')
ylabel('$u_2$','Interpreter','LaTex')
plot(y(1,:),y(2,:))
%% Point in region II
beta0=1.9130;
alpha0=-0.9008;
% Point near the steady-state
x2=0.1;
% Integrate
sol = dde23(@(t,y,Z) funcs.sys_rhs([y,Z],[beta0 alpha0 tau0]),...
    tau0,[x1 x2],[0 tfinal]);
t=linspace(0,tfinal,1000);
y=deval(sol,t);
% Plot
figure(2);clf;
title('Point in region II')
xlabel('$u_1$','Interpreter','LaTex')
ylabel('$u_2$','Interpreter','LaTex')
plot(y(1,:),y(2,:))
%% Point in region III
beta0=1.8890;
alpha0=-0.6081;
% Orbit converging to periodic orbit
x2=0.1;
% Integrate
tfinal=3000; % use lager time interval
sol = dde23(@(t,y,Z) funcs.sys_rhs([y,Z],[beta0 alpha0 tau0]),...
    tau0,[x1 x2],[0 tfinal]);
t=linspace(0,tfinal,4000);
y1=deval(sol,t);
% Plot
figure(3);clf;
title('Point in region III')
xlabel('$u_1$','Interpreter','LaTex')
ylabel('$u_2$','Interpreter','LaTex')
plot(y1(1,:),y1(2,:))
% Orbit converging to the stable steady-state
hold on;
% Integrate
x2=0.093;
sol = dde23(@(t,y,Z) funcs.sys_rhs([y,Z],[beta0 alpha0 tau0]),...
    tau0,[x1 x2],[0 tfinal]);
y2=deval(sol,t);
% Add to plot
plot(y2(1,:),y2(2,:),'Color',cm(2,:));
%% Time series of the previous solutions in region III
figure(4);clf;
title('Time series of solutions in region III')
plot(t,y1(1,:),t,y2(1,:))
xlabel('$t$','Interpreter','LaTex')
ylabel('$u_1$','Interpreter','LaTex')