%% Zero-Hopf demo - generate right-hand side and derivatives from symbolic expressions
%
% <html>
% $Id$
% </html>
% 
% From:
% S. Ma and Z. Feng,
% Fold-Hopf bifurcation of the Rose-Hindmarsh model with time delay,
% Appl. Math. Comput., 219 (2013), pp. 10073--10081
% https://www.sciencedirect.com/science/article/abs/pii/S0096300313004141
%
%% Differential equations
%
% $$\dot{x}=y - ax^3 + bx^2(t-tau) - cz + I_{app},$$
%
% $$\dot{y}=c - dx^2-y,$$
%
% $$\dot{z]=r(S(x-\chi)-z).$$
%% Add paths and load sym package if GNU Octave is used
clear
ddebiftoolpath='../../../../';
addpath(strcat(ddebiftoolpath,'ddebiftool'),...
    strcat(ddebiftoolpath,'ddebiftool_extra_symbolic'));
if dde_isoctave()
    pkg load symbolic
end
%% Create parameter names as strings and define fixed parameters
% The demo has the parameters |beta|, |alpha| and |tau|
parnames={'Iapp','S','r','tau'};
a=sym(1.0,'r');
b=sym(3.0,'r');
c=sym(1.0,'r');
d=sym(5.0,'r');
chi=sym(-1.6,'r');
%% Create symbols for parameters, states and delays states
% The array |par| is the array of symbols in the same order as parnames.
% Due to the following two lines we may, for example, use either tau or
% par(1) to refer to the delay.
syms(parnames{:});       % create symbols for beta, alpha and tua
par=cell2sym(parnames);  % now beta is par(1) etc
%% Define system using symbolic algebra
syms x xtau y ytau z ztau % create symbols for u1(t) u1(t-tau), u2(t), u2(t-tau)
dx_dt=y-a*x^3+b*xtau^2-c*z+Iapp;
dy_dt=c-d*x^2-y;
dz_dt=r*(S*(x-chi)-z);
%% Differentiate and generate code, exporting it to sym_RH_mf (multi-linear forms)
[fstr,derivs]=dde_sym2funcs(...
    [dx_dt;dy_dt;dz_dt],... % n x 1 array of derivative symbolic expressions
    [x,xtau;y,ytau;z,ztau],... % n x (ntau+1) array of symbols for states (current & delayed)
    par,... % 1 x np (or np x 1) array of symbols used for parameters
    'filename','sym_RH_mf',... % optional argument specifying output file
    'directional_derivative',false); 
%% Differentiate and generate code, exporting it to sym_RH (directional derivatives)
[fstr,derivs]=dde_sym2funcs(...
    [dx_dt;dy_dt;dz_dt],... % n x 1 array of derivative symbolic expressions
    [x,xtau;y,ytau;z,ztau],... % n x (ntau+1) array of symbols for states (current & delayed)
    par,...            % 1 x np (or np x 1) array of symbols used for parameters
    'filename','sym_RH',...  % optional argument specifying output file
    'directional_derivative',true);
%% The following trick creates an index struct |ind| of parameter names
% without manual counting. The indices can be accessed as |ind.tau|,
% |ind.a|, |ind.b|, |ind.d|. We save the structure and the names in a mat
% file
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
save('RH_parnames.mat','-v7','parnames','ind');