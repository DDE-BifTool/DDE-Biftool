%% DDE-BIFTOOL demo for nested state-dependent delays
%% Creation of right-hand side and derivatives using symbolic toolbox
%
% <html>
% $Id: gen_sym_nested_demo.m 362 2019-07-14 15:49:40Z jansieber $
% </html>
%
% This demo shows that DDE-Biftool can treat DDEs with nested
% state-dependent delays. The DDE is
% 
% $$\dot x(t)=-x(t-\tau_k)$$
%
% where $\tau_0=0$ and 
%
% $$\tau_{j+1}=\tau_j+x(t-\tau_j)+p x(t-\tau_j)^2$$
%
% The system has two parameters, $\tau_1$ and $p$. A Hopf bifurcation
% occurs at %\tau_1=\pi/2$, which is sub- or supercritical, depending on
% $p$ and $k$.
%%
clear
addpath('../../ddebiftool',...
    '../../ddebiftool_extra_symbolic');
if dde_isoctave()
    pkg load symbolic
end
%% Set number of delays and parameter names
ntau=3; % plays the role of k
parnames={'tau1','p'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
%% Define system using symbolic algebra
% using arbitrary variable names
x=sym('x',[1,ntau+1]);
syms(parnames{:});
par=sym(parnames);
f=-x(ntau+1);
tau=sym('tau',[ntau+1,1]);
tau(1)=tau1;
for i=1:ntau
    tau(i+1)=tau1+x(i)+p*x(i)^2;
end
%% Differentiate and generate code, exporting it to sym_nested
[fstr,derivs]=dde_sym2funcs(f,x,par,'sd_delay',tau(2:ntau+1),...
    'directional_derivative',true,'filename','sym_nested');
%% Bifurcation analysis
% DDE-Biftool itself does not rely on the symbolic toolbox, but only on its
% generated code. See <demo1_simple.html> for the demo.

