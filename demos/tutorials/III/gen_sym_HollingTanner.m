%% Holling-Tanner model - Generate right-hand side and derivatives from symbolic expressions
% This script demonstrates that variable and parameter names in the
% creation of the right-hand side do not have to be xx or par.
%
%% Differential equations for predator-prey model
%
% From:
% Liu, X., Liu, Y., and Wang, J. (2013). Bogdanov-Takens bifurcation of a 
% delayed ratio-dependent Holling-Tanner predator prey system. 
% In Abstract and Applied Analysis, volume 2013
% 
% $$x'=(x+m)(1-x-m)-\frac{xy}{ay+x}$$ 
%
% $$y'=\delta y\left(\beta-\frac{y(t-\tau)}{x(t-\tau)}\right)$$
% 
% <html>
% $Id: gen_sym_HollingTanner.m 178 2017-03-14 00:00:01Z jansieber $
% </html>
%%
clear
addpath('../../../ddebiftool',...
    '../../../ddebiftool_extra_symbolic');
%% Set number of delays and parameter names
ntau=1;
parnames={'beta','tau','a','m','h','delta'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
%% Define system using symbolic algebra
% using arbitrary variable names
x=sym('x',[1,ntau+1]);
y=sym('y',[1,ntau+1]);
syms(parnames{:});
par=sym(parnames);
f=[(x(1)+m)*(1-x(1)-m)-...
    x(1)*y(1)/(a*y(1)+x(1))-h;...
    delta*y(1)*(beta-y(2)/x(2))];
%% Differentiate and generate code, exporting it to sym_HollingTanner.m
[fstr,derivs]=dde_sym2funcs(f,[x;y],par,'filename','sym_HollingTanner',...
    'directional_derivative',true);
%% Differentiate and generate code, exporting it to sym_HollingTanner.m
[fstr,derivs]=dde_sym2funcs(f,[x;y],par,'filename','sym_HollingTanner_mf',...
    'directional_derivative',false);
%% Bifurcation analysis
% DDE-Biftool itself does not rely on the symbolic toolbox, but only on its
% generated code. See <HollingTanner_demo.html> for the demo.
