ddebiftoolpath='../../';
addpath(strcat(ddebiftoolpath,'ddebiftool'),...
    strcat(ddebiftoolpath,'ddebiftool_extra_psol'),...
    strcat(ddebiftoolpath,'ddebiftool_extra_symbolic'),...
    strcat(ddebiftoolpath,'ddebiftool_extra_nmfm'),...
    strcat(ddebiftoolpath,'ddebiftool_utilities'));
format compact
format short g
%% define parameter names and values
if dde_isoctave()
    pkg load symbolic
end
ntau=1;
parnames={'Iapp','S','tau','a','b','c','d','chi','r'};
syms(parnames{:});
par=dde_sym_from_cell(parnames);
syms x y z xt yt zt
rose_hindmarsh_sys=[...
    y-a*x^3+b*xt^2-c*z+Iapp;...
    c-d*x^2-y;...
    r*(S*(x-chi)-z)];
%% Differentiate and generate code, exporting it to sym_RH
[fstr,derivs]=dde_sym2funcs(rose_hindmarsh_sys,...
    [x,xt;y,yt;z,zt],par,'filename','sym_RH','keeptemp',true,'diretional_derivative',true);
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
save('RH_symdefs.mat','parnames','ind');
assert(exist('sym_RH.m','file')==2);
