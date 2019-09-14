%% Hopf-transcritical demo - generate right-hand side and derivatives from symbolic expressions
%
% <html>
% $Id$
% </html>
% 
% From:
%  J. Bramburger, B. Dionne, and V. G. LeBlanc, 
% Zero-Hopf bifurcation in the Van der Pol oscillator with delayed position 
% and velocity feedback, Nonlinear Dynam., 78 (2014), pp. 2959–2973,
% http://dx.doi.org/10.1007/s11071-014-1638-0.
%
%% Differential equations
%
% $$\dot{x}=(\tau_0+\mu_2) y,$$
%
% $$\dot{y}=(\tau_0+\mu_2)[-x-\epsilon(x^2-1)y+(1+\mu_1)x(t-1)-0.2y(t-1)].$$
% $$ \qquad\qquad-0.2x^2(t-1)-0.2x(t-1)y(t-1)-0.2y^2(t-1)+0.5x^3(t-1)$$
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
parnames={'mu1','mu2','tau'};
epsilon=sym(0.3,'r');
a=sym(-0.2,'r');
omega0=sqrt(2-epsilon^2+a^2);
tau0=1/omega0*acos((1-(1+epsilon*a)*omega0^2)/(a^2*omega0^2+1));
%% Create symbols for parameters, states and delays states
% The array |par| is the array of symbols in the same order as parnames.
% Due to the following two lines we may, for example, use either tau or
% par(1) to refer to the delay.
syms(parnames{:});       % create symbols for beta, alpha and tua
par=cell2sym(parnames);  % now beta is par(1) etc
%% Define system using symbolic algebra
syms x xtau y ytau % create symbols for u1(t) u1(t-1), u2(t), u2(t-1)
dx_dt=(tau0+mu2)*y;
dy_dt=(tau0+mu2)*(...
    (1+mu1)*xtau-0.2*ytau-0.2*xtau^2-...
    0.2*xtau*ytau-0.2*ytau^2 ...
    -epsilon*(x^2-1)*y-x+0.5*x^3);
%% Differentiate and generate code, exporting it to sym_FHN_mf (multi-linear forms)
[fstr,derivs]=dde_sym2funcs(...
    [dx_dt;dy_dt],... % n x 1 array of derivative symbolic expressions
    [x,xtau;y,ytau],... % n x (ntau+1) array of symbols for states (current & delayed)
    par,... % 1 x np (or np x 1) array of symbols used for parameters
    'filename','sym_vdpo_mf',... % optional argument specifying output file
    'directional_derivative',false); 
%% Differentiate and generate code, exporting it to sym_FHN_mf (directional derivatives)
[fstr,derivs]=dde_sym2funcs(...
    [dx_dt;dy_dt],... % n x 1 array of derivative symbolic expressions
    [x,xtau;y,ytau],... % n x (ntau+1) array of symbols for states (current & delayed)
    par,...            % 1 x np (or np x 1) array of symbols used for parameters
    'filename','sym_vdpo',...  % optional argument specifying output file
    'directional_derivative',true);
%% The following trick creates an index struct |ind| of parameter names
% without manual counting. The indices can be accessed as |ind.tau|,
% |ind.a|, |ind.b|, |ind.d|. We save the structure and the names in a mat
% file
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
save('vdpo_parnames.mat','-v7','parnames','ind');