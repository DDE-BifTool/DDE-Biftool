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
addpath('new_nmfm/');
rmpath('../../ddebiftool_extra_nmfm/');
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
funcs=set_funcs('sys_rhs',rhs,'sys_ntau',sys_ntau,'sys_tau',tau,...
    'sys_deri',drhs,'sys_dtau',dtau,'x_vectorized',true);
par_ini=[0, 2.3, 1.3, 6, 4.75, 1];
indkappa1=1;
indkappa2=2;
%% Trivial equilibrium branch
% The equilibrium $x=0$ changes its stability in Hopf bifurcations
[eqbr,suc]=SetupStst(funcs,'contpar',indkappa1,'x',0,'parameter',par_ini,...
    'max_bound',[indkappa1,12],'max_step',[indkappa1,0.1]);
if ~suc
    error('equilibrium not found');
end
clf
eqbr=br_contn(funcs,eqbr,100);
% stability
[eqnunst,dom,triv_defect,eqbr.point]=...
GetStability(eqbr,'funcs',funcs,'points',2:length(eqbr.point));
%% Periodic orbits branching off from 1st Hopf bifurcation
% We continue the periodic orbits in |parameter(1)| (|kappa1|), and calculate
% their stability.
indhopf=find(eqnunst>0,1,'first');
[hopf,suc]=SetupHopf(funcs,eqbr,indhopf,'contpar',[indkappa1,indkappa2],...
    'corpar',indkappa1,'dir',indkappa2,'step',1e-2);
if ~suc
    error('initialization of periodic orbits failed');
end
%%
hopf.parameter.max_step=[indkappa1,0.1; indkappa2,0.1];
hopf.parameter.max_bound=[indkappa1,8;indkappa2,5];
hopf=br_contn(funcs,hopf,200);
disp('calculate stability')
[hopfnunst,dom,triv_defect,hopf.point]=...
    GetStability(hopf,'exclude_trivial',true,'funcs',funcs); %#ok<*ASGLU>
fprintf('maximum error of trivial eigenvalue: %g\n',max(abs(triv_defect)));
%%
%load('humphries1dbif.mat');
[L1,L1low]=HopfLyapunovCoefficients(funcs,hopf);
k1=arrayfun(@(x)x.parameter(indkappa1),hopf.point);
%% Find Hopf hopf interaction
indhopf2=find(hopfnunst==2,1,'first')-1;
[hoho,hoho_low]=HopfHopfNormalform(funcs,hopf,indhopf2+(0:1));
%%
save('humphries1dbif.mat');
