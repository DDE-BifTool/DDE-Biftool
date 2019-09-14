%% test some basic examples for fold and Hopf bifurcations
clear
addpath('../../ddebiftool',...
    '../../ddebiftool_extra_symbolic');
if dde_isoctave()
    pkg load symbolic
end
%% Variables for basic fold
ntau=2;
parnames={'p','tau1','tau2'};
x=sym('x',[2,ntau+1]);
syms(parnames{:});
par=dde_sym_from_cell(parnames);
%% Equations
f=[p+(tau1-1)*x(2,3)-x(2,2)^2*x(1,3);x(1,1)-x(2,1)];
%% Files
[fstr,derivs]=dde_sym2funcs(f,x,par,'filename','sym_basic_fold','directional_derivative',true);
%% Hopf with state-dependent delays
ntau=2;
parnames={'p','c'};
x=sym('x',[2,ntau+1]);
syms(parnames{:});
par=dde_sym_from_cell(parnames);
%% Equations
f=[p-x(2,2)*x(1,3);x(1,1)-x(2,1)];
