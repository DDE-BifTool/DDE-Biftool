%% test NDDE modifications - generation of r.h.s.
%
% This file demonstrates how one can use the symbolic toolbox to
% automatically create the right-hand side and its first directional
% derivatives using the symbolic toolbox, using the routines in
% |ddebiftool_extra_symbolic|. The example is the laser model by Hessel et
% al.
%
% $Id: gen_sym_ddae.m 372 2019-08-27 22:13:34Z jansieber $
%%
clear
base=[pwd(),'/../../'];
addpath([base,'ddebiftool'],[base,'ddebiftool_extra_symbolic']);
if dde_isoctave()
    pkg load symbolic
end
%% Set number of delays and create parameter names as strings
parnames={'a','b','h','tau'};
%% Create symbols for parameters, states and delayed states
% The array |par| is the array of symbols in the same order as parnames.
% Due to the following two lines we may, for example, use either tau or
% par(1) to refer to the delay.
syms(parnames{:});       % create symbols for parameters
par=cell2sym(parnames);  % now tau is par(1) etc
for i=1:length(parnames)
    assume(par(i),'real');
end
syms Er Ei Ertau Eitau Yr Yi Yrtau Yitau real % create symbols for E(t) E(t-tau), n(t)
E=Er+1i*Ei;
Etau=Ertau+1i*Eitau;
Y=Yr+1i*Yi;
Ytau=Yrtau+1i*Yitau;
dE_dt=(a-(1+1i)*conj(E)*E)*E+h*Y;
Y_a=b*(Etau-Ytau)-Y;
dEr_dt=real(dE_dt);
dEi_dt=imag(dE_dt);
Yr_a=real(Y_a);
Yi_a=imag(Y_a);
%% Differentiate and generate code, exporting it to sym_LangKobayashi
[fstr,derivs]=dde_sym2funcs(...
    [dEr_dt;dEi_dt;Yr_a;Yi_a],... % n x 1 array of derivative symbolic expressions
    [Er,Ertau; Ei,Eitau; Yr,Yrtau; Yi,Yitau],... % n x (ntau+1) array of symbols for states (current & delayed)
    par,...            % 1 x np (or np x 1) array of symbols used for parameters
    'filename','sym_ddae_demo',... % optional argument specifying output file
    'directional_derivative',true,...
    'maxorder',2);
