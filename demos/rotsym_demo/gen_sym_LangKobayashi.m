%% DDE-BIFTOOL Demo for continuation of local bifurcations with rotational symmetry
%% Create a right-hand side from symbolic expressions
%
% <html>
% $Id: gen_sym_LangKobayashi.m 354 2019-06-30 23:15:16Z jansieber $
% </html>
%
% This file demonstrates how one can use the symbolic toolbox to
% automatically create the right-hand side and its first directional
% derivatives using the symbolic toolbox, using the routines in
% |ddebiftool_extra_symbolic|. The example is the Lang-Kobayashi equation.
clear
base=[pwd(),'/../../'];
addpath([base,'ddebiftool'],[base,'ddebiftool_extra_symbolic']);
if dde_isoctave()
    pkg load symbolic
end
%% Set number of delays and create parameter names as strings
parnames={...
    'pump',...    % injection current
    'eta',...     % feedback strength
    'phi',...     % feedback phase
    'tau',...     %  delay
    'alpha',...   % alpha factor
    'epsilon'};     % rotation velocity
%% Create symbols for parameters, states and delayed states
% The array |par| is the array of symbols in the same order as parnames.
% Due to the following two lines we may, for example, use either tau or
% par(1) to refer to the delay.
syms(parnames{:});       % create symbols for parameters
par=cell2sym(parnames);  % now tau is par(1) etc
for i=1:length(parnames)
    assume(par(i),'real');
end
syms Er Ei Ertau Eitau n ntau real % create symbols for E(t) E(t-tau), n(t)
E=Er+1i*Ei;
Etau=Ertau+1i*Eitau;
dE_dt=(1+1i*alpha)*n*E+eta*exp(1i*phi)*Etau;
dn_dt=epsilon*(pump-n-(2*n+1)*conj(E)*E);
dEr_dt=real(dE_dt);
dEi_dt=imag(dE_dt);
dn_dt=simplify(dn_dt);
%% Differentiate and generate code, exporting it to sym_LangKobayashi
[fstr,derivs]=dde_sym2funcs(...
    [dEr_dt;dEi_dt;dn_dt],... % n x 1 array of derivative symbolic expressions
    [Er,Ertau; Ei,Eitau; n,ntau],... % n x (ntau+1) array of symbols for states (current & delayed)
    par,...            % 1 x np (or np x 1) array of symbols used for parameters
    'filename','sym_LangKobayashi',... % optional argument specifying output file
    'directional_derivative',true,...
    'maxorder',2);
