%% Create right-hand side and derivatives for demo Humphriesetal from symbolic expressions
%
% <html>
% $Id$
% </html>
%
% Humphries et al ( A. R. Humphries,  O. A. DeMasi,
% F. M. G. Magpantay,  F. Upham (2012), _Dynamics of a delay differential equation
% with multiple state-dependent delays_, Discrete and Continuous Dynamical
% Systems 32(8) pp. 2701-2727 <http://dx.doi.org/10.3934/dcds.2012.32.2701>)
%
% Consider a scalar DDE with two state-dependent delays:
% 
% $$x'(t)=-\gamma x(t)-\kappa_1 x(t-a_1-c x(t))-\kappa_2 x(t-a_2-c x(t))\mbox{.}$$
%%
clear
close all
addpath('../../ddebiftool/',...
    '../../ddebiftool_extra_symbolic/');
%% Define delays, parameter names and states
% |cind| is a cell array used to automatically number the parameter names.
% It is used to create the structure |ind|, which contains pointers into
% the parameter vector (so, |par(ind.beta)| is parameter beta).
ntaus=2;
parnames={'kappa1','kappa2','a1','a2','gamma','c'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
xx=sym('x',[1,ntaus+1]);
syms(parnames{:});
par=sym(parnames);
%% Right-hand side f and delay functions
f=-gamma*xx(1,1)-kappa1*xx(1,2)-kappa2*xx(1,3);
tau=[a1+c*xx(1,1);a2+c*xx(1,1)];
%% Generate code
[fstr,derivs]=dde_sym2funcs(f,xx,par,'sd_delay',tau,'filename','sym_humphriesetal',...
    'directional_derivative',true);
