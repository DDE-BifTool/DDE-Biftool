%% Demo for state-dependent delays --- example from Humpries etal (DCDS-A 2012)
%
% <html>
% $Id: humphriesetal_demo.m 176 2017-03-13 00:25:33Z jansieber $
% </html>
%
% Humphries et al ( A. R. Humphries,  O. A. DeMasi,
% F. M. G. Magpantay,  F. Upham (2012), _Dynamics of a delay differential equation
% with multiple state-dependent delays_, Discrete and Continuous Dynamical
% Systems 32(8) pp. 2701-2727 <http://dx.doi.org/10.3934/dcds.2012.32.2701>)
%
% consider a scalar DDE with two state-dependent delays:
% 
% $$x'(t)=-\gamma x(t)-\kappa_1 x(t-a_1-c x(t))-\kappa_2 x(t-a_2-c x(t))\mbox{.}$$
%
% The system exhibits double Hopf bifurcations, folds of periodic orbits
% and torus bifurcations, which can be computed and demonstrated using
% DDE-Biftool.
%%
%% Add path to DDE-Biftool
clear
close all
addpath('../../ddebiftool/',...
    '../../ddebiftool_extra_psol/',...
    '../../ddebiftool_extra_nmfm/',...
    '../../ddebiftool_utilities/');
%% Function definitions (right-hand side and delays)
% Define the right-hand side and the state-dependent delays. We also provide their
% derivatives (in separate functions) to speed up execution.
% 
%% Order of the parameters: 
% |p(1:2)=kappa1:2, p(3:4)=a1:2, p(5)=gamma, p(6)=c|
parnames={'kappa1','kappa2','a1','a2','gamma','c'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
%% Options for defining the problem
% below are several ways in which one can define the problem. One can
% choose one of them by setting |funcs=...| at the end.
%%  (Alternative 1) Minimal - Right-hand side and delays
% All derivatives are approximated by finite differences
ntaus=2;
rhs=@(x,p)-p(ind.gamma)*x(1,1,:)-p(ind.kappa1)*x(1,2,:)-p(ind.kappa2)*x(1,3,:);
sys_ntau=@()ntaus;
tau=@(nr,x,p)p(ind.a1-1+nr)+p(ind.c)*x(1,1,:);
fnum=set_funcs('sys_rhs',rhs,'sys_ntau',sys_ntau,'sys_tau',tau,'x_vectorized',true);
%% (Alternative 2) traditional provision of first- and 2nd-order derivatives
% This is sufficient for all standard bifurcation analysis. Normal form
% computations require nested higher-order derivatives.
drhs=@(x,p,nx,np,v)sys_deri_humphries_etal(rhs,x,p,nx,np,v);
dtau=@(nr,x,p,nx,np)sys_dtau_humphries_etal(tau,nr,x,p,nx,np);
ftraditional=set_funcs('sys_rhs',rhs,'sys_ntau',sys_ntau,'sys_tau',tau,...
    'sys_deri',drhs,'sys_dtau',dtau,'x_vectorized',true);
%% (Alternative 3) provide directional derivatives of rhs and tau
% Only 1st- and 2nd-order is needed (not used for normal form computations).
rhs_dir{1}=@(x,p,dx,dp)...
    -dp(ind.gamma)*x(1,1,:)-p(ind.gamma)*dx(1,1,:)...
    -dp(ind.kappa1)*x(1,2,:)-p(ind.kappa1)*dx(1,2,:)...
    -dp(ind.kappa2)*x(1,3,:)-p(ind.kappa2)*dx(1,3,:);
rhs_dir{2}=@(x,p,dx,dp)...
    -2*dp(ind.gamma)*dx(1,1,:)-2*dp(ind.kappa1)*dx(1,2,:)-2*dp(ind.kappa2)*dx(1,3,:);
tau_dir{1}=@(nr,x,p,dx,dp)dp(ind.a1-1+nr)+dp(ind.c)*x(1,1,:)+p(ind.c)*dx(1,1,:);
tau_dir{2}=@(nr,x,p,dx,dp)2*dp(ind.c)*dx(1,1,:);
fdirectional=set_funcs('sys_rhs',rhs,'sys_ntau',sys_ntau,'sys_tau',tau,...
    'sys_dirderi',rhs_dir,'sys_dirdtau',tau_dir,'x_vectorized',true);
%% (Alternative 5) create all derivatives with symbolic toolbox
% The code generation script is <gen_sym_humphriesetal.html>. The function
% |set_symfuncs| is a wrapper around |set_funcs|. The call does not need to
% specify the number of delays (as this has been given in
% gen_sym_humphriesetal.html>).
fsymbolic=set_symfuncs(@sym_humphriesetal);
%% First step: continuation of equilibria and their bifurcations and normal forms
% see <humphriesetal_equilibria.html>.
save('humphriesetal_demo_funcs_results.mat');
