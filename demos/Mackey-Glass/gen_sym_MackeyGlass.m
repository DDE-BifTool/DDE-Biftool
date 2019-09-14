%% DDE-Biftool demo Mackey-Glass Equation 
%% Create right-hand side and derivatives from symbolic expressions
%
% The Mackey-Glass equation is given by
% 
% $$x'(t)=\beta \frac{x(t-\tau)}{1+x(t-\tau)^n}-\gamma x(t)$$
% 
% Parameters are (in this order) |beta|, |n|, |tau| (|gamma| is not part of
% parameter vector).
%
% <html>
% $Id$
% </html>
%%
clear
addpath('../../ddebiftool',...
    '../../ddebiftool_extra_symbolic');
%% Set number of delays and parameter names
ntau=1;
parnames={'beta','n','tau','gamma'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
%% Define system using symbolic algebra
% using arbitrary variable names
x=sym('x',[1,ntau+1]);
syms(parnames{:});
par=sym(parnames);
f=beta*x(2)/(1+x(2)^n)-gamma*x(1);
%% Differentiate and generate code, exporting it to sym_mackeyGlass.m
[fstr,derivs]=dde_sym2funcs(f,x,par,'filename','sym_MackeyGlass_mf','directional_derivative',false);
%% Bifurcation analysis
% DDE-Biftool itself does not rely on the symbolic toolbox, but only on its
% generated code. See <MackeyGlass_demo.html> for the demo.
