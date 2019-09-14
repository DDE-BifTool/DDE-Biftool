function sol = mackey_glass(n,tspan,t)
%%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: mackey_glass.m 134 2016-09-12 11:10:44Z mmbosschaert $
%%

% close all;
% tspan=[0;150];

% set parameters
tau=2;
beta=2;
gamma=1;
% n=7;
par=[beta gamma n];

% options = ddeset('MaxStep',0.01);
options = ddeset('RelTol',1e-4);
sol = dde23(@(t,y,Z)ddefun_ex2([y,Z],par),tau,@(t,beta,gamma,n) 0.5,tspan,options);

figure;
y=deval(sol,t);
ylag=deval(sol,t-2);
plot(y,ylag);
xlabel('$x(t)$','Interpreter','LaTex');
ylabel('$x(t-\tau)$','Interpreter','LaTex');

end
