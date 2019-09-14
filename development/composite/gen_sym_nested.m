clear
addpath('../../ddebiftool/',...
    '../../ddebiftool_extra_nmfm/',...
    '../../ddebiftool_extra_rotsym/',...
    '../../ddebiftool_utilities/',....
    '../../ddebiftool_extra_symbolic/');
ntaus=1;
parnames={'p'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
xx=sym('x',[1,ntaus+1]);
syms(parnames{:});
par=sym(parnames);
%% Right-hand side f and delay functions
f=-xx(1,2);
tau=p+xx(1,1);
%% Generate code
[fstr,derivs]=dde_sym2funcs(f,xx,par,'sd_delay',tau,'filename','sym_nested');
