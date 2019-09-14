%% test NDDE modifications - generation of r.h.s.
%
% This file demonstrates how one can use the symbolic toolbox to
% automatically create the right-hand side and its first directional
% derivatives using the symbolic toolbox, using the routines in
% |ddebiftool_extra_symbolic|. The example is the laser model by Hessel et
% al.
clear
base=[pwd(),'/../../'];
addpath([base,'ddebiftool'],[base,'ddebiftool_extra_symbolic']);
if dde_isoctave()
    pkg load symbolic
end
%% Set number of delays and create parameter names as strings
parnames={'Rc','R','Is_Rc_by_V0','tau_by_C'};
%% Create symbols for parameters, states and delayed states
% The array |par| is the array of symbols in the same order as parnames.
% Due to the following two lines we may, for example, use either tau or
% par(1) to refer to the delay.
syms(parnames{:});       % create symbols for parameters
par=cell2sym(parnames);  % now tau is par(1) etc
syms y v ytau vtau       % create symbols varaiables
u=y-vtau;
dy_dt=(-(1/Rc-1/R)*u+(1/Rc+1/R)*vtau);
v_a=v-u+Is_Rc_by_V0*(exp(v+u)-1);
%% Differentiate and generate code, exporting it to sym_LangKobayashi
[fstr,derivs]=dde_sym2funcs(...
    [dy_dt;v_a],...     % n x 1 array of derivative symbolic expressions
    [y,ytau;v,vtau],... % n x (ntau+1) array of symbols for states (current & delayed)
    par,...            % 1 x np (or np x 1) array of symbols used for parameters
    'filename','sym_tline',... % optional argument specifying output file
    'directional_derivative',true);
