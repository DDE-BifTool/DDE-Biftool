%% Add paths and load sym package if GNU Octave is used
clear
ddebiftoolpath='../../';
addpath(strcat(ddebiftoolpath,'ddebiftool'),...
        strcat(ddebiftoolpath,'ddebiftool_extra_symbolic'));
%% Create parameter names as strings and define fixed parameters
% The demo has the parameters |tau|
parnames={'tau'};
%% Create symbols for parameters, states and delays states
% The array |par| is the array of symbols in the same order as parnames.
% Due to the following two lines we may, for example, use either tau or
% par(1) to refer to the delay.
syms(parnames{:});       % create symbol for tua
par=cell2sym(parnames);  % now a is par(1) etc
%% Define system using symbolic algebra
syms x1 x2 x1tau x2tau  % create symbols for x(t) x(t-tau), y(t), y(t-tau), ect.
%% Define the system
A0 = [-5 1; 2 6];
A1 = [-2 1; 4 -1];
dx_dt = A0*[x1; x2] + A1*[x1tau; x2tau];
%% Differentiate and generate code, exporting it to sym_strangeExample_mf (multi-linear forms)
[fstr,derivs]=dde_sym2funcs(...
    [dx_dt(1); dx_dt(2)],... % n x 1 array of derivative symbolic expressions
    [x1,x1tau;x2,x2tau],... % n x (ntau+1) array of symbols for states (current & delayed)
    par,... % 1 x np (or np x 1) array of symbols used for parameters
    'filename','sym_GyanSwarupNag_mf'... % optional argument specifying output file
);