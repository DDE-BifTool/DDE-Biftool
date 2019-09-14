%%  Active control system
% 
% Generate right-hand side and derivatives from symbolic expressions
%
%% Differential equations for active control system with delayed feedback
%
% From:
% Yuting Ding, Jun Cao, Weihua Jiang. Double Hopf bifurcation in active 
% control system with delayed feedback: application to glue dosing 
% processes for particleboard. In Nonlinear Dynamics. 83, 1567--1576 (2016)
% 
% $$x'=(x+m)(1-x-m)-\frac{xy}{ay+x}$$ 
%
% $$y'=\delta y\left(\beta-\frac{y(t-\tau)}{x(t-\tau)}\right)$$
%% Add paths and load sym package if GNU Octave is used
clear
ddebiftoolpath='../../../../';
addpath(strcat(ddebiftoolpath,'ddebiftool'),...
    strcat(ddebiftoolpath,'ddebiftool_extra_symbolic'));
if dde_isoctave()
    pkg load symbolic
end
%% Create parameter names as strings and define fixed parameters
% The demo has the parameters |zeta| and |tau_scaled|
parnames={'zeta','tau','tau_scaled'};
gu=sym(0.1,'r');
gv=sym(0.52,'r');
beta=sym(0.1,'r');
%% Create symbols for parameters, states and delays states
% The array |par| is the array of symbols in the same order as parnames.
% Due to the following two lines we may, for example, use either tau or
% par(1) to refer to the delay.
syms(parnames{:});       % create symbols for beta, alpha and tua
par=cell2sym(parnames);  % now beta is par(1) etc
%% Define system using symbolic algebra
syms x xtau y ytau % create symbols for x(t) x(t-tau), y(t), y(t-tau)
dx_dt=tau*y;
dy_dt=tau*(-x-gu*xtau-2*zeta*y-gv*ytau+beta*xtau^3);
%% Differentiate and generate code, exporting it to sym_acs_mf (multi-linear forms)
[fstr,derivs]=dde_sym2funcs(...
    [dx_dt;dy_dt],... % n x 1 array of derivative symbolic expressions
    [x,xtau;y,ytau],... % n x (ntau+1) array of symbols for states (current & delayed)
    par,... % 1 x np (or np x 1) array of symbols used for parameters
    'filename','sym_acs_mf',... % optional argument specifying output file
    'directional_derivative',false);
%% Differentiate and generate code, exporting it to sym_acs (directional_derivative)
[fstr,derivs]=dde_sym2funcs(...
    [dx_dt;dy_dt],... % n x 1 array of derivative symbolic expressions
    [x,xtau;y,ytau],... % n x (ntau+1) array of symbols for states (current & delayed)
    par,... % 1 x np (or np x 1) array of symbols used for parameters
    'filename','sym_acs',... % optional argument specifying output file
    'directional_derivative',true); 
%% The following trick creates an index struct |ind| of parameter names
% without manual counting. The indices can be accessed as |ind.tau|,
% |ind.a|, |ind.b|, |ind.d|. We save the structure and the names in a mat
% file
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
save('acs_parnames.mat','-v7','parnames','ind');