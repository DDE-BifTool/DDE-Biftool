clear
addpath('../../ddebiftool',...
    '../../ddebiftool_extra_psol',...
    '../../ddebiftool_extra_nmfm',...
    '../../ddebiftool_utilities');
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
%% Differentiate and generate code, exporting it to sym_Holling.m
[fstr,derivs]=dde_sym2funcs(f,[x;y],par,'filename','sym_Holling');
