%% generate rhs for position control
clear
addpath('../../ddebiftool/',...
    '../../ddebiftool_extra_nmfm/',...
    '../../ddebiftool_extra_psol/',...
    '../../ddebiftool_utilities/',....
    '../../ddebiftool_extra_symbolic/');
%% Set number of delays and parameter names
ntau=3;
parnames={'tau0','s0','K','c'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
x=sym('x',[1,ntau+1]);
s=sym('s',[1,ntau+1]);
syms(parnames{:});
par=sym(parnames);
f=[-K*c/2*(s(2)-s0);...
    (-K*c/2*(s(4)+s(2)-2*s0)-c*s(1)+x(3)+x(1))/(c+K*c/2*(s(4)-s0))];
tau=[tau0;s(1);tau0+s(1)];
%% Generate code
[fstr,derivs]=dde_sym2funcs(f,[x;s],par,'sd_delay',tau,'filename','sym_poscontrol');


