%%  System files for tutorial VI
% 
% Generate right-hand side and derivatives from symbolic expressions
%
%% Differential equations for neuron model
%
% This system represents a non-dimensionalized model of two interacting
% layers of neurons.
% 
% $$ \dot{x}_1(t) = - x_1(t) - a g(bx_1(t-\tau_1)) + cg (dx_2(t-\tau_2)), $$
%
% $$ \dot{x}_2(t) = - x_2(t) - a g(bx_2(t-\tau_1)) + cg (dx_1(t-\tau_2)), $$
%
% where $g:\bf{R} \rightarrow \bf{R}$ is of the sigmoidal form
%
% $$ g(z) = \left[\tanh(z-1) + \tanh(1)\right]\cosh(1)^2. $$
% 
% The variables $x_1(t)$ and $x_2(t)$ represent the population-averaged
% neural activity at time $t$ in layers one and two, respectively.
% The parameter $a > 0$ is a measure of the strength of inhibitory feedback,
% while $c > 0$ measures the strength of the excitatory effect of one layer
% on the other. The parameters $b > 0$ and $d > 0$ are saturation rates and
% the delays $\tau_{1,2}$ represent time lags in the inhibitory feedback
% loop and excitatory inter-layer connection. Note that the system is
% symmetric with respect to interchanging the labels $1$ and $2$, so
% equilibria are necessarily of the form $(x_0,x_0)$.
%% Add paths and load sym package if GNU Octave is used
clear
ddebiftoolpath='../../../';
addpath(strcat(ddebiftoolpath,'ddebiftool'),...
    strcat(ddebiftoolpath,'ddebiftool_extra_symbolic'));
if dde_isoctave()
    pkg load symbolic
end
%% Create parameter names as strings and define fixed parameters
% The demo has the parameters |zeta| and |tau_scaled|
parnames={'a','b','c','d','tau1','tau2'};
%% Create symbols for parameters, states and delays states
% The array |par| is the array of symbols in the same order as parnames.
% Due to the following two lines we may, for example, use either tau1 or
% par(5) to refer to the delay.
syms(parnames{:});       % create symbols for beta, alpha and tua
par=cell2sym(parnames);  % now beta is par(1) etc
%% Define system using symbolic algebra
syms x1 x1tau1 x1tau2 x2 x2tau1 x2tau2 % create symbols for x(t) x(t-tau1), y(t), y(t-tau1)
g=@(z) (tanh(z-1)+tanh(1))*cosh(1)^2;
dx1_dt=-x1-a*g(b*x1tau1)+c*g(d*x2tau2);
dx2_dt=-x2-a*g(b*x2tau1)+c*g(d*x1tau2);
%% Differentiate and generate code, exporting it to sym_acs_mf (multi-linear forms)
[fstr,derivs]=dde_sym2funcs(...
    [dx1_dt;dx2_dt],... % n x 1 array of derivative symbolic expressions
    [x1,x1tau1,x1tau2;x2,x2tau1,x2tau2],... % n x (ntau+1) array of symbols for states (current & delayed)
    par,... % 1 x np (or np x 1) array of symbols used for parameters
    'filename','sym_neuron_mf',... % optional argument specifying output file
    'directional_derivative',false);
%% Differentiate and generate code, exporting it to sym_acs (directional_derivative)
[fstr,derivs]=dde_sym2funcs(...
    [dx1_dt;dx2_dt],... % n x 1 array of derivative symbolic expressions
    [x1,x1tau1,x1tau2;x2,x2tau1,x2tau2],... % n x (ntau+1) array of symbols for states (current & delayed)
    par,... % 1 x np (or np x 1) array of symbols used for parameters
    'filename','sym_neuron',... % optional argument specifying output file
    'directional_derivative',true); 
%% The following trick creates an index struct |ind| of parameter names
% without manual counting. The indices can be accessed as |ind.a|,
% |ind.b|, |ind.c|, etc. We save the structure and the names in a mat
% file
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
save('neuron_parnames.mat','-v7','parnames','ind');