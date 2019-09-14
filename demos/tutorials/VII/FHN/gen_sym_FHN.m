%% Generaluzed-Hopf demo - generate right-hand side and derivatives from symbolic expressions
%
% <html>
% $Id$
% </html>
% 
% From:
% Zhen, Bin (PRC-TONG-AEM); Xu, Jian [Xu, Jian2] (PRC-TONG-AEM)
% Bautin bifurcation analysis for synchronous solution of a coupled FHN
% neural system with delay.
% Commun. Nonlinear Sci. Numer. Simul. 15 (2010), no. 2, 442–458.
%
%% Differential equations
%
% $$u_1'=-\frac{u_1^3+(c+\alpha)u_1^2+d u_1-u_2+2\beta*tanh(u_1(t-\tau))$$
%
% $$u_2'=\epsilon(u_1-b u_2)$$
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
parnames={'beta','alpha','tau'};
b=sym(0.9,'r');
epsilon=sym(0.08,'r');
c=sym(2.0528,'r');
d=sym(-3.2135,'r');
%% Create symbols for parameters, states and delays states
% The array |par| is the array of symbols in the same order as parnames.
% Due to the following two lines we may, for example, use either beta or
% par(1) to refer to the delay.
syms(parnames{:});       % create symbols for beta, alpha and tua
par=cell2sym(parnames);  % now beta is par(1) etc
%% Define system using symbolic algebra
syms u1 u1t u2 u2t % create symbols for u1(t) u1(t-tau), u2(t), u2(t-tau)
du1_dt=-u1^3/3+(c+alpha)*u1^2+d*u1-u2+2*beta*tanh(u1t);
du2_dt=epsilon*(u1-b*u2);
%% Differentiate and generate code, exporting it to sym_FHN_mf (multi-linear forms)
[fstr,derivs]=dde_sym2funcs(...
    [du1_dt;du2_dt],... % n x 1 array of derivative symbolic expressions
    [u1,u1t;u2,u2t],... % n x (ntau+1) array of symbols for states (current & delayed)
    par,... % 1 x np (or np x 1) array of symbols used for parameters
    'filename','sym_FHN_mf',... % optional argument specifying output file
    'directional_derivative',false); 
%% Differentiate and generate code, exporting it to sym_FHN (directional derivatives)
[fstr,derivs]=dde_sym2funcs(...
    [du1_dt;du2_dt],... % n x 1 array of derivative symbolic expressions
    [u1,u1t;u2,u2t],... % n x (ntau+1) array of symbols for states (current & delayed)
    par,...            % 1 x np (or np x 1) array of symbols used for parameters
    'filename','sym_FHN',...  % optional argument specifying output file
    'directional_derivative',true);
%% The following trick creates an index struct |ind| of parameter names
% without manual counting. The indices can be accessed as |ind.beta|,
% |ind.alpha| and |ind.tau|. We save the structure and the names in a mat
% file
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
save('FHN_parnames.mat','-v7','parnames','ind');