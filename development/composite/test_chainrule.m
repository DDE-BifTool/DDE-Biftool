%% test chain rule formula
clear
addpath('../../ddebiftool_extra_nmfm/',...
    '../../ddebiftool/');
%% define test points and functions
x0=[2;1];
dx=[1;2];
f=@(y)[sin(0.5*y(1));cos(y(2));cos(y(1)*y(2)/4)];
g=@(x)[cos(x(1)*x(2)/2);sin(x(1))];
%% use numerical finite differences
gh=@(h)g(x0+h*dx);
fhy=@(h,dy)f(g(x0)+h*dy);
fgh=@(h)f(gh(h));
args={'isvectorized',false};
dg=@(k)nmfm_dfdx_scalar(gh,k,args{:});
dfg=@(k)nmfm_dfdx_scalar(fgh,k,args{:});
dfhy=@(k,dy)nmfm_dfdx_scalar(@(h)fhy(h,dy),k,args{:});
%% test derivatives
mxorder=5;
cf=dde_chainrule_combinatorics(mxorder);
for i=mxorder:-1:1
    cfp{i}=dde_pol_from_chain(cf(i));
end
nf=length(f(g(x0)));
for i=mxorder:-1:1
    a_dfg(:,i)=dfg(i);
    a_cdfg(:,i)=nmfm_chainrule(dfhy,dg,i,cfp{i},nf);
end
