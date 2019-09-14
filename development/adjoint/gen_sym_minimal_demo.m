%% DDE-BIFTOOL Minimal demo for continuation of local bifurcations of periodic orbits
%% Create right-hand side from symbolic expressions
%
% <html>
% $Id$
% </html>
%
% This demo illustrates how to track local bifurcations of periodic orbits
% in DDEs with contant delay, using the extension |ddebiftool_extra_psol|
% and the auxiliary functions in |ddebiftool_utilities|.
% The example is the Duffing oscillator with delayed feedback discussed in
% the large-delay limit by Yanchuk & Perlikowski in (PRE79,0462211,2009):
%
% $$x''(t)+d*x'(t)+a*x(t)+x^3+b*[x(t)-x(t-\tau)]=0$$
%
% The parameters are $(\tau,a,b,d)$ (used in this order in the |parameter|
% vector).
%%
clear
addpath('../../ddebiftool',...
    '../../ddebiftool_extra_symbolic');
%% Set number of delays and parameter names
ntau=1;
parnames={'tau','a','b','d'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
%% Define system using symbolic algebra
% using arbitrary variable names
x=sym('x',[1,ntau+1]);
xp=sym('xp',[1,ntau+1]);
syms(parnames{:});
par=sym(parnames);
f=[xp(1);-d*xp(1)-a*x(1)-x(1)^3-b*(x(1)-x(2))];
%% Differentiate and generate code, exporting it to sym_HollingTanner.m
[fstr,derivs]=dde_sym2funcs(f,[x;xp],par,'filename','sym_minimal_demo');
%% Bifurcation analysis
% DDE-Biftool itself does not rely on the symbolic toolbox, but only on its
% generated code. See <minimal_demo.html> for the start of the demo.
