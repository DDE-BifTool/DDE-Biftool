%% DDE-BIFTOOL Minimal demo for continuation of local bifurcations of periodic orbits
%% Create a right-hand side from symbolic expressions
%
% <html>
% $Id: gen_sym_minimal_demo.m 363 2019-07-14 17:20:58Z jansieber $
% </html>
%
% This file demonstrates how one can use the symbolic toolbox
% to automatically create the right-hand side and all of its derivatives
% using the symbolic toolbox, using the routines in
% |ddebiftool_extra_symbolic|. The example is the Duffing oscillator with
% delayed feedback discussed in the large-delay limit by Yanchuk &
% Perlikowski in (PRE79,0462211,2009):
%
% $$x''(t)+d x'(t)+a x(t)+x^3+b [x(t)-x(t-\tau)]=0$$
%
% The parameters are $(\tau,a,b,d)$ (used in this order in the |parameter|
% vector).
%% Add paths and load sym package if octave is used
clear
addpath('../../ddebiftool',...
    '../../ddebiftool_extra_symbolic');
if dde_isoctave()
    pkg load symbolic
end
%% Set number of delays and create parameter names as strings
% The demo has the parameters |tau|, |a|, |b| and |d|.
parnames={'tau','a','b','d'};
%% Create symbols for parameters, states and delays states
% The array |par| is the array of symbols in the same order as parnames.
% Due to the following two lines we may, for example, use either tau or
% par(1) to refer to the delay.
syms(parnames{:});       % create symbols for tau, a, b, d
par=cell2sym(parnames);  % now tau is par(1) etc
syms x xt xp xpt         % create symbols for x(t) x(t-tau), x'(t), x'(t-tau)
%% Define system using symbolic algebra
% using the variables created above. The left-hand side will be a symbol,
% with the value given by the expression on the right-hand side if the
% right-hand side contains symbols.
dx_dt = xp;
dxp_dt=-d*xp-a*x-x^3-b*(x-xt);
%% Differentiate and generate code, exporting it to sym_minimal_demo
[fstr,derivs]=dde_sym2funcs(...
    [dx_dt;dxp_dt],... % n x 1 array of derivative symbolic expressions
    [x,xt; xp,xpt],... % n x (ntau+1) array of symbols for states (current & delayed)
    par,...            % 1 x np (or np x 1) array of symbols used for parameters
    'filename','sym_minimal_demo',... % optional argument specifying output file
    'directional_derivative',true);
%% Detailed description of |dde_sym2funcs|
% The function |dde_sym2funcs| has the header
%
%   function [funcstr,derivatives]=dde_sym2funcs(fs,xxs,ps,...) 
%
% It converts the input fs (an array of symbolc expressions) into a
% right-hand side function and its derivatives that can be used for
% DDE-Biftool. The wrapper |set_symfuncs| around |set_funcs| can be used to
% create a structure containing right-hand sides and derivatives.
%
%% Inputs:
%
% * |fs|:  n x 1 array of symbolic expressions, defining the right-hand side
% * |xxs|: n x (ntau+1) array of symbols for states, |xxs(j,k)| is the
% $x_j(t-\tau_{k-1})$, and |xxs(j,1)|=$x_j(t)$.
% * |par|: 1 x np (or np x 1) array of symbols for parameters, |par(j)| is
% the |j|the parameter.
%
%% Common optional name-value input pairs
%
% * |'sd_delay'| (default |sym([]))|): ntau x 1 (or 1 x ntau) array tau of
% symbolic expressions depending on symbols in |xxs| and |par|. Entry
% |tau(k)| is delay number |k| ($\tau_{k}(x,p)$) and may depend on
% |xxs(:,nu)| for |nu=1..k|.
% * |'write'| (default |true|): write output to file?
% * |'filename'| (default |'sys'|): results are written to function file
% with this name.
% * |'maxorder'| (default 5): maximal order of derivatives to be computed.
%
%% Outputs (often not needed)
%
% * |fstr|: string containing the code for writing the function.
% * |derivatives|: symbolic expressions for all dervatives (a structure,
% containing the fields |df|, |xx|, |parameter|, |dx| and |dp|. |df|
% contains the expressions, |xx|, |parameter|, |dx| and |dp| the symbols
% used in these expressions. 
%
% The routine will create a vectorized matlab
% function |y=f(xx,p)| from the expressions in |fs|, and its directional
% derivatives up to a maximum order (default 5). For state-dependent delays
% (optional argument |sd_delay| non-empty), it also computes the values and
% directional derivatives of all delays.
%% Warning
% The function dde_sym2funcs will write to temporary files, since
% matlabFunction can only write to files.
% 
%% The following trick creates an index struct |ind| of parameter names
% without manual counting. The indices can be accessed as |ind.tau|,
% |ind.a|, |ind.b|, |ind.d|. We save the structure and the names in a mat
% file
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
save('minimal_demo_parnames.mat','-v7','parnames','ind');
%% Bifurcation analysis
% DDE-Biftool itself does not rely on the symbolic toolbox, but only on its
% generated code. See <minimal_demo.html> for the start of the demo.
