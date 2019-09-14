%% DDE-BIFTOOL  state-dependent delays sd-demo
%% Generation of right-hand side and its derivatives using symbolic toolbox
% 
% <html>
% $Id: gen_sym_sd_demo.m 362 2019-07-14 15:49:40Z jansieber $
% </html>
% 
%% Differential equations
% The differential equations for this example are
%
% $$\begin{array}{l}
% \frac{\mathrm{d}}{\mathrm{d} t}x_1(t)=\frac{1}{p_1+x_2(t)}\left(1-p_2x_1(t)x_1(t-\tau_3)
% x_3(t-\tau_3)+p_3x_1(t-\tau_1)x_2(t-\tau_2)\right),\\
% \frac{\mathrm{d}}{\mathrm{d} t}x_2(t)=\frac{p_4 x_1(t)}{p_1+x_2(t)}+
%         p_5\tanh(x_2(t-\tau_5))-1,\\
% \frac{\mathrm{d}}{\mathrm{d} t}x_3(t)=p_6(x_2(t)-x_3(t))-p_7(x_1(t-\tau_6)-x_2(t-\tau_4))e^{-p_8 \tau_5},\\
% \frac{\mathrm{d}}{\mathrm{d} t}x_4(t)=x_1(t-\tau_4)e^{-p_1 \tau_5} -0.1,\\
% \frac{\mathrm{d}}{\mathrm{d} t}x_5(t)=3(x_1(t-\tau_2)-x_5(t))-p_9,
% \end{array} $$
%
% where $\tau_1$ and $\tau_2$ are constant delays and
%
% $$
% \begin{array}{l}
% \tau_3=2+p_5\tau_1x_2(t)x_2(t-\tau_1),\\
% \tau_4=1-\frac{1}{1+x_1(t)x_2(t-\tau_2)},\\
% \tau_5=x_4(t),\\
% \tau_6=x_5(t).
% \end{array}
% $$
%
% This system has five components $(x_1,\ldots,x_5)$, six delays
% $(\tau_1,\ldots,\tau_6)$ and eleven parameters $(p_1,\ldots,p_{11})$,
% where $p_{10}=\tau_1$ and $p_{11}=\tau_2$. 
%
%% Add paths for symbolic routines
clear
addpath('../../ddebiftool',...
    '../../ddebiftool_extra_symbolic');
if dde_isoctave()
    pkg load symbolic
end
%% Set number of delays and parameter names
ntau=6;
parnames=[strcat('p',num2cell('1':'9')),{'tau1','tau2'}];
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
%% Introduce symbols for x variables, delays and parameters
x=sym('x',[5,ntau+1]);
syms(parnames{:});
par=sym(parnames);
%% Right-hand side
f=[...
    (1/(p1+x(2,1)))*(1-p2*x(1,1)*x(1,4)*x(3,4)+p3*x(1,2)*x(2,3));...
    p4*x(1,1)/(p1+x(2,1))+p5*tanh(x(2,6))-1;...
    p6*(x(2,1)-x(3,1))-p7*(x(1,7)-x(2,5))*exp(-p8*x(4,1));...
    x(1,5)*exp(-p1*x(4,1))-0.1;...
    3*(x(1,3)-x(5,1))-p9];
%% Delay
delays=[...
    tau1;...
    tau2;...
    2+p5*tau1*x(2,1)*x(2,2);...
    1-1/(1+x(2,3)*x(1,1));...
    x(4,1);...
    x(5,1)];
%% Differentiate and generate code, exporting it to sym_sd_demo
[fstr,derivs]=dde_sym2funcs(f,x,par,'sd_delay',delays,...
    'directional_derivative',true,'filename','sym_sd_demo');
%% Bifurcation analysis
% DDE-Biftool itself does not rely on the symbolic toolbox, but only on its
% generated code. See <sd_demo.html> for the demo.
