%% DDE-Biftool demo -- Connecting orbits
%% Create a right-hand side from symbolic expressions
%
% <html>
% $Id: gen_sym_hom_demo.m 375 2019-09-14 14:26:50Z jansieber $
% </html>
%
% This file demonstrates how one can use the symbolic toolbox
% to automatically create the right-hand side and all of its derivatives
% yusing the symbolic toolbox, using the routines in
% |ddebiftool_extra_symbolic|. 
%
% The right-hand side is for the demo |hom_demo|, showing the comptuation
% of connecting orbits.
% The model is given as a two-dimensional system of differential equations
%
% $$
% \begin{array}[t]{rcl}
% \dot{x}_{1}(t)&=&-x_{1}(t)+q_{11}\frac{1}{1+\mathrm{e}^{-4x_{1}(t-\tau)}}
%   -q_{12}x_2(t-\tau) +e_1\\
% \dot{x}_{2}(t)&=&-x_2(t)+q_{21}\frac{1}{1+e^{-4x_{1}(t-\tau)}}+e_2
% \end{array}
% $$
%
% The parameters of the system are $[q_{11},q_{12}, q_{21}, e_1, e_2,\tau]$ (in
% this order in the |parameter| vector). Our main bifurcation parameters
% are $q_{12}$ (|q12|) and $e_1$ (|e1|).
%
%% Add paths and load sym package if octave is used
clear
addpath('../../ddebiftool',...
    '../../ddebiftool_extra_symbolic');
if dde_isoctave()
    pkg load symbolic
end
%% Set number of delays and create parameter names as strings
parnames={'q11','q12','q21','e1','e2','tau'};
syms(parnames{:});       % create symbols for parameters
par=cell2sym(parnames);  % now q11 is par(1) etc
x=sym('x',[2,2]);        % create symbols for x1(t) x1(t-tau), x2(t), x2(t-tau)
g=1/(1+exp(-4*x(1,2)));
f=[...
    -x(1,1)+q11*g-q12*x(2,2)+e1;...
    -x(2,1)+q21*g+e2];
%% Differentiate and generate code, exporting it to sym_minimal_demo
[fstr,derivs]=dde_sym2funcs(...
    f,...              % n x 1 array of derivative symbolic expressions
    x,...              % n x (ntau+1) array of symbols for states (current & delayed)
    par,...            % 1 x np (or np x 1) array of symbols used for parameters
    'filename','sym_hom_demo','directional_derivative',true); % optional argument specifying output file


