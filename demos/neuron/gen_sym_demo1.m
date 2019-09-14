%% DDE-BIFTOOL demo 1 - Neuron, 
%% Creation of right-hand side and derivatives using symbolic toolbox
%
% <html>
% $Id: gen_sym_demo1.m 367 2019-07-14 21:56:46Z jansieber $
% </html>
%
% This demo repeats the illustrative example using the convenience
% functions in |ddebiftool_utilities|. The system of
% delay differential equations is [Shay,99] again
% 
% $$\left\{\begin{array}{l}\dot{x_1}(t)=
% -\kappa x_1(t)+\beta \tanh(x_1(t-\tau_s))+a_{12}\tanh(x_2(t-\tau_2)) \\ 
% \dot{x_2}(t)=-\kappa x_2(t)+\beta \tanh(x_2(t-\tau_s))+a_{21}\tanh(x_1(t-\tau_1)).
% \end{array}\right.$$
%
% This system models two coupled neurons with time delayed connections.
% It has two components ($x_1$ and $x_2$), three delays 
% ($\tau_1$, $\tau_2$ and $\tau_s$), and four other parameters 
% ($\kappa$, $\beta$, $a_{12}$ and $a_{21}$).
%%
clear
addpath('../../ddebiftool',...
    '../../ddebiftool_extra_symbolic');
if dde_isoctave()
    pkg load symbolic
end
%% Set number of delays and parameter names
ntau=3;
parnames={'kappa','beta','a12','a21','tau1','tau2','taus'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
%% Define system using symbolic algebra
% using arbitrary variable names
x=sym('x',[2,ntau+1]);
syms(parnames{:});
par=sym(parnames);
f=[...
    -kappa*x(1,1)+beta*tanh(x(1,4))+a12*tanh(x(2,3));...
    -kappa*x(2,1)+beta*tanh(x(2,4))+a21*tanh(x(1,2))];
%% Differentiate and generate code, exporting it to sym_neuron
[fstr,derivs]=dde_sym2funcs(f,x,par,'filename','sym_neuron');
%% Bifurcation analysis
% DDE-Biftool itself does not rely on the symbolic toolbox, but only on its
% generated code. See <demo1_simple.html> for the demo.

