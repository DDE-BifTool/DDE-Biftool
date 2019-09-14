%% test branch point detection
%
% $Id$
%%
%% Add paths and load sym package if octave is used
clear
addpath('../../ddebiftool',...
    '../../ddebiftool_extra_symbolic');
if dde_isoctave()
    pkg load symbolic
end
%% Set number of delays and create parameter names as strings
% The demo has the parameters |tau|, |a|, |b| and |d|.
parnames={'p','q','tau'};
%% Create symbols for parameters, states and delays states
% The array |par| is the array of symbols in the same order as parnames.
% Due to the following two lines we may, for example, use either tau or
% par(1) to refer to the delay.
syms(parnames{:});       % create symbols for tau, a, b, d
par=cell2sym(parnames);  % now tau is par(1) etc
syms x y xt yt           % create symbols for x(t) x(t-tau), y(t), y(t-tau)
%% Define system using symbolic algebra
% using the variables created above. The left-hand side will be a symbol,
% with the value given by the expression on the right-hand side if the
% right-hand side contains symbols.
dx_dt = p*x+x^2*y;
dy_dt= q-yt-y^3;
%% Differentiate and generate code, exporting it to sym_minimal_demo
[fstr,derivs]=dde_sym2funcs(...
    [dx_dt;dy_dt],... % n x 1 array of derivative symbolic expressions
    [x,xt; y,yt],...  % n x (ntau+1) array of symbols for states (current & delayed)
    par,...           % 1 x np (or np x 1) array of symbols used for parameters
    'filename','sym_bp_demo','directional_derivative',true); % optional argument specifying output file
