%% Clean workspace and add DDE-BifTool scripts to search path
clear;      % clear variables
close all;	% close figures
ddebiftoolpath='../../';
addpath(strcat(ddebiftoolpath,'ddebiftool'),...
    strcat(ddebiftoolpath,'ddebiftool_extra_psol'),...
    strcat(ddebiftoolpath,'ddebiftool_extra_symbolic'),...
    strcat(ddebiftoolpath,'ddebiftool_extra_nmfm'),...
    strcat(ddebiftoolpath,'ddebiftool_utilities'));
%% Set parameter names
parnames={'tau'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
%% Set the funcs structure using sysbolic generated derivatives
% We load the precalculated multilinear forms. These have been generated
% with the file gen_sym_bandlimitedfeedback.m.
funcs=set_symfuncs(@sym_GyanSwarupNag_mf, 'sys_tau',@() ind.tau);
%% Construct steady-state point
% construct steady-state point
stst=dde_stst_create('x', [0;0]);
tau = 1;
stst.parameter(ind.tau) = tau;
% Calculate stability
method_stst=df_mthod(funcs,'stst');
method_stst.stability.minimal_real_part=-10;
stst.stability=p_stabil(funcs, stst, method_stst.stability);
plot(stst.stability.l1, '*')
stst.stability.l1(1:3)