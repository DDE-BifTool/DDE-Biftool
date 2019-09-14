%%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: tutorial_I.m 134 2016-09-12 11:10:44Z mmbosschaert $
%%

%% 2.1 Example
close all;
tspan=[0;4];
hspan=-1:0.01:0;
lags=1;

% solve the DDE
sol = dde23(@(t,y,Z)ddefun_ex1([y,Z]),lags,@history_ex1,tspan);

% plot solution
plot([hspan, sol.x],[history_ex1(hspan), sol.y])
hold

xlabel('$t$','Interpreter','LaTex')
ylabel('$x$','Interpreter','LaTex')

sol = dde23(@(t,y,Z)ddefun_ex1([y,Z]),lags,@(t) exp(t),tspan);
plot([hspan, sol.x],[exp(hspan), sol.y],':')

sol = dde23(@(t,y,Z)ddefun_ex1([y,Z]),lags,@(t) 1-t,tspan);
plot([hspan, sol.x],[1-hspan, sol.y],'--')

legend({'$\cos(t)$','$\exp(t)$','$1-t$'},'Interpreter','latex')

%% 3.1 Example with Hopf bifurcation
sys_rhs = @(xx,mu) -(pi*sqrt(3)/9+mu)*(xx(1,2)+xx(1,3))*(1-xx(1,1)^2);

% just beore the Hopf bifurcation
mu=-0.1;
history1=-0.99;
options = ddeset('MaxStep',0.1);
sol1 = dde23(@(t,y,Z)sys_rhs([y,Z],mu),[1 2],history1,[0 100],options);

% time series plot
figure
plot(sol1.x,sol1.y)
xlabel('$t$','Interpreter','LaTex')
ylabel('$x$','Interpreter','LaTex')

% plot in phase-plane
figure
plot(sol1.y,sol1.yp)
xlabel('$x$','Interpreter','LaTex')
ylabel('$\dot{x}$','Interpreter','LaTex')

% at critical value
mu=0;
history1=-0.99;
sol1 = dde23(@(t,y,Z)sys_rhs([y,Z],mu),[1 2],history1,[0 100],options)

% time series plot
figure
plot(sol1.x,sol1.y)
xlabel('$t$','Interpreter','LaTex')
ylabel('$x$','Interpreter','LaTex')

% plot in phase-plane
figure
plot(sol1.y,sol1.yp)
xlabel('$x$','Interpreter','LaTex')
ylabel('$\dot{x}$','Interpreter','LaTex')

% just after the Hopf bifurcation 
mu=0.1;
sol1 = dde23(@(t,y,Z)sys_rhs([y,Z],mu),[1 2],-0.99,[0 100],options)
sol2 = dde23(@(t,y,Z)sys_rhs([y,Z],mu),[1 2],0.01,[0 100],options)

% time series plot
figure
plot(sol1.x,sol1.y)
hold
plot(sol2.x,sol2.y)
xlabel('$t$','Interpreter','LaTex')
ylabel('$x$','Interpreter','LaTex')

% plot in phase-plane
figure
plot(sol1.y,sol1.yp)
hold
plot(sol2.y,sol2.yp)
xlabel('$x$','Interpreter','LaTex')
ylabel('$\dot{x}$','Interpreter','LaTex')

%% 3.2 Mathematical Model of Hematological Diseases
mackey_glass(7,[0;150],linspace(100,150,1000)); % (a) 1 
mackey_glass(7.75,[0;200],linspace(150,200,1000)); % (b) 2
mackey_glass(8.5,[0;250],linspace(200,250,1000)); % (c) 3 
mackey_glass(8.79,[0;400],linspace(300,400,1000)); % (d) 4
mackey_glass(9.65,[0;600],linspace(300,600,3000)); % (e) 5
mackey_glass(9.697148,[0;800],linspace(300,800,6000)); % (f) 6
mackey_glass(9.6975,[0;400],linspace(300,400,5000)); % (g) 7
mackey_glass(9.76,[0;400],linspace(300,400,3000)); % (h) 8
mackey_glass(10.0,[0;600],linspace(300,600,3000)); % (i) 9
mackey_glass(20.0,[0;400],linspace(300,400,1000)); % (j) 10

%% time series
tau=2;
beta=2;
gamma=1;
tspan=0:300;
t=linspace(200,300,1000);

h=figure;

n=7;
par=[beta gamma n];
sol = dde23(@(t,y,Z)ddefun_ex2([y,Z],par),tau,@(t,beta,gamma,n) 0.5,tspan);
y=deval(sol,t);
hold
plot(t,y)

n=7.75;
par=[beta gamma n];
sol = dde23(@(t,y,Z)ddefun_ex2([y,Z],par),tau,@(t,beta,gamma,n) 0.5,tspan);
y=deval(sol,t);
plot(t,y)

set(h,'Position',[91 1508 1920 508]);

xlabel('$t$','Interpreter','LaTex')
ylabel('$x$','Interpreter','LaTex')

%%
h2=figure;

n=8.5;
par=[beta gamma n];
sol = dde23(@(t,y,Z)ddefun_ex2([y,Z],par),tau,@(t,beta,gamma,n) 0.5,tspan);
t=linspace(200,300,1000);
y=deval(sol,t);

hold
plot(t,y)

n=8.79;
par=[beta gamma n];
sol = dde23(@(t,y,Z)ddefun_ex2([y,Z],par),tau,@(t,beta,gamma,n) 0.5,tspan);
y=deval(sol,t);
plot(t,y)

set(h2,'Position',[91 1508 1920 508]);

xlabel('$t$','Interpreter','LaTex');
ylabel('$x$','Interpreter','LaTex');

%% Example 3.3 Sinusoidal nonlinearity
%% Hopf bifurcation
figure
sys_rhs = @(xx,par) sin(xx(1,2));

% just before the Hopf bifurcation
par=pi/2-0.1;
history1=3.1;
options = ddeset('MaxStep',0.1);
sol1 = dde23(@(t,y,Z)sys_rhs([y,Z],par),par,history1,[0 200],options);

plot(sol1.x,sol1.y)
xlabel('$t$','Interpreter','LaTex')
ylabel('$x$','Interpreter','LaTex')

% just after the Hopf bifurcation
par=pi/2+0.1;
history1=3.1;
history2=4.5;
sol1 = dde23(@(t,y,Z)sys_rhs([y,Z],par),par,history1,[0 200],options);
sol2 = dde23(@(t,y,Z)sys_rhs([y,Z],par),par,history2,[0 200],options);

figure
plot(sol1.x,sol1.y)
hold
plot(sol2.x,sol2.y)
xlabel('$t$','Interpreter','LaTex')
ylabel('$x$','Interpreter','LaTex')

%% two coexisting cycles
tspan=[0 500];
par=4;
history1=2.988;
history2=2.99;
options = ddeset('MaxStep',0.1);
sol1 = dde23(@(t,y,Z)sys_rhs([y,Z],par),par,history1,tspan,options)
sol2 = dde23(@(t,y,Z)sys_rhs([y,Z],par),par,history2,tspan,options)

figure
% times series
plot(sol1.x,sol1.y)
hold
plot(sol2.x,sol2.y)
xlabel('$t$','Interpreter','LaTex');
ylabel('$x$','Interpreter','LaTex');

% (x,x')-phase space
figure
n1=4000;
plot(sol1.y(n1:end),sol1.yp(n1:end))
hold
plot(sol2.y(n1:end),sol2.yp(n1:end))
xlabel('$x(t)$','Interpreter','LaTex');
ylabel('$\dot{x}(t)$','Interpreter','LaTex');
axis([-1     7    -1.2     1.2])

% (x(t),x(t-\tau))-phase space
t=linspace(400,500,1000);
y1=deval(sol1,t);
y2=deval(sol2,t);
y1lag=deval(sol1,t-par);
y2lag=deval(sol2,t-par);

figure
plot(y1,y1lag)
hold
plot(y2,y2lag)
xlabel('$x(t)$','Interpreter','LaTex')
ylabel('$x(t-\tau)$','Interpreter','LaTex')

%% before period doubling
par1=4.82;
history1=2.99;
% options = ddeset('RelTol',1e-7);
sol1 = dde23(@(t,y,Z)sys_rhs([y,Z],par1),par1,history1,[0 3000])

t=linspace(2600,3000,1000);
y1=deval(sol1,t);
y1lag=deval(sol1,t-par1);

figure
plot(y1,y1lag)

xlabel('$x(t)$','Interpreter','LaTex')
ylabel('$x(t-\tau)$','Interpreter','LaTex')

% after period doubling
par2=4.87;
history2=2.99;
options = ddeset('RelTol',1e-05);
sol2 = dde23(@(t,y,Z)sys_rhs([y,Z],par2),par2,history2,[0 3000],options)

y2=deval(sol2,t);
y2lag=deval(sol2,t-par2);

figure
plot(y2,y2lag)
xlabel('$x(t)$','Interpreter','LaTex')
ylabel('$x(t-\tau)$','Interpreter','LaTex')
