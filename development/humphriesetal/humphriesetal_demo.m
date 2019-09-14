%% Demo for state-dependent delays --- example from Humpries etal (DCDS-A 2012)
%
% <html>
% $Id: humphriesetal_demo.m 49 2014-05-02 17:43:07Z jan.sieber $
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
% Order of the parameters: 
% |p(1:2)=kappa1:2, p(3:4)=a1:2, p(5)=gamma, p(6)=c|
rhs=@(x,p)-p(5)*x(1,1,:)-p(1)*x(1,2,:)-p(2)*x(1,3,:);
sys_ntau=@()2;
tau=@(nr,x,p)p(2+nr)+p(6)*x(1,1,:);
drhs=@(x,p,nx,np,v)sys_deri_humphries_etal(rhs,x,p,nx,np,v);
dtau=@(nr,x,p,nx,np)sys_dtau_humphries_etal(tau,nr,x,p,nx,np);
%% Set the following switch to false to test numerical finite differences 
% for normal form computations.
symbolic_derivatives=true;
if symbolic_derivatives 
    dirderi=@humphriesetal_dirderi; 
else
    dirderi=[]; %#ok<UNRCH>
end
funcs=set_funcs('sys_rhs',rhs,'sys_ntau',sys_ntau,'sys_tau',tau,...
    'sys_deri',drhs,'sys_dtau',dtau,'sys_dirderi',dirderi,'x_vectorized',true);
par_ini=[0,2.3,1.3,6,4.75,1];
indkappa1=1;
indkappa2=2;
%% First step: continuation of equilibria and their bifurcations and normal forms
% see <humphriesetal_equilibria.html>.
save('humphriesetal_demo_funcs_results.mat');
