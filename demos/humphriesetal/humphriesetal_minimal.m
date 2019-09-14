%% Equilibrium bifurcations and normal forms for example from Humpries etal (DCDS-A 2012)
%
% <html>
% $Id: humphriesetal_minimal.m 374 2019-09-14 14:02:58Z jansieber $
% </html>
%
% Humphries et al ( A. R. Humphries,  O. A. DeMasi,
% F. M. G. Magpantay,  F. Upham (2012), _Dynamics of a delay differential equation
% with multiple state-dependent delays_, Discrete and Continuous Dynamical
% Systems 32(8) pp. 2701-2727 <http://dx.doi.org/10.3934/dcds.2012.32.2701>)
%
% This is the fast and convenient version using normal form branching
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
addpath('../../ddebiftool/',...
    '../../ddebiftool_extra_psol/',...
    '../../ddebiftool_extra_nmfm/',...
    '../../ddebiftool_utilities/');
%% Function definitions (right-hand side and delays)
% Define the right-hand side and the state-dependent delays. We also provide their
% derivatives (in separate functions) to speed up execution.
funcs=set_symfuncs(@sym_humphriesetal);
%% Order of the parameters: 
% |p(1:2)=kappa1:2, p(3:4)=a1:2, p(5)=gamma, p(6)=c|
parnames={'kappa1','kappa2','a1','a2','gamma','c'};
cind=[parnames;num2cell(1:length(parnames))];
ip=struct(cind{:});
par_ini([ip.kappa1,ip.kappa2,ip.a1,ip.a2,ip.gamma,ip.c])=...
    [      0,        2.3      1.3,   6,    4.75,  1];
%% Set continuation bounds during two-parameter continuation
brbounds={...
    'max_bound',[ip.kappa1,13; ip.kappa2,5],...
    'min_bound',[ip.kappa1,0; ip.kappa2,0],...
    'max_step',[ip.kappa1,0.1; ip.kappa2,0.1]};
%% Trivial equilibrium branch
% The equilibrium $x=0$ changes its stability in Hopf bifurcations
[eqbr,suc]=SetupStst(funcs,'contpar',ip.kappa1,'x',0,'parameter',par_ini,...
    'max_bound',[ip.kappa1,14],'max_step',[ip.kappa1,0.1]);
if ~suc
    error('equilibrium not found');
end
figure(1);clf;ax1=gca;xlabel(ax1,'kappa1');ylabel(ax1,'x');
eqbr=br_contn(funcs,eqbr,200);
%% Stability, bifurcations and their normal forms
[eqbr,~,ind_hopf,eqc1biftype]=LocateSpecialPoints(funcs,eqbr);
%%
figure(2);clf;ax2=gca;
Plot2dBranch(eqbr,'y',@(p)p.x(1),'ax',ax2);
set(ax2,'ylim',[-0.1,1]);xlabel(ax2,'kappa1');ylabel(ax2,'max x');
%% First Hopf bifurcation in two parameters
% Continuation parameters are $\kappa_1$ and $\kappa_2$.
nhopf=length(ind_hopf);
figure(3);clf;ax3=gca;xlabel(ax3,'kappa1');ylabel(ax3,'kappa2');
nhopfpoints=100;
fprintf('Computation of 1st Hopf bifurcation found\n');
[hbranch0,suc]=SetupHopf(funcs,eqbr,ind_hopf(1),...
    'contpar',[ip.kappa1,ip.kappa2],'dir',ip.kappa1,'step',0.1,brbounds{:});
hbranch0=br_contn(funcs,hbranch0,nhopfpoints,'plotaxis',ax3);
hbranch0=br_rvers(hbranch0);
hbranch0=br_contn(funcs,hbranch0,nhopfpoints,'plotaxis',ax3);
[hbranch0,~,indhoho0,eqc2biftype0]=LocateSpecialPoints(funcs,hbranch0);
%% Compute all branches coming out of Hopf-Hopf interaction
t_args={'submesh','cheb','degree',5,'intervals',50,...
    'collocation_parameters','cheb','print_residual_info',1,...
    'stop_1_2',true,'minimal_angle',0.5};
nbr2_steps=200;
C1info=BranchFromCodimension2(funcs,hbranch0,t_args{:},brbounds{:});
clear C1branches testfuncs;
figure(4);clf;ax4=gca;xlabel(ax4,'kappa1');ylabel(ax4,'kappa2');hold(ax4,'on');
Plot2dBranch(hbranch0,'ax',ax4);
for i=length(C1info):-1:1
    fprintf('Continuing %s\n',[C1info(i).C1type]);
    C1branches(i)=br_contn(C1info(i).funcs,C1info(i).branch,nbr2_steps,'plotaxis',ax3);
    if C1info(i).reverse
        fprintf('Reversing %s\n',[C1info(i).C1type]);
        C1branches(i)=br_rvers(C1branches(i));
        C1branches(i)=br_contn(C1info(i).funcs,C1branches(i),nbr2_steps,'plotaxis',ax3);
    end
    [C1_nunst{i},C1dom{i},C1defect{i},C1branches(i).point]=...
        GetStability(C1branches(i),'funcs',C1info(i).funcs,'exclude_trivial',true);
    [C1branches(i),testfuncs{i},indc2{i},ind2ctype{i}]=...
        LocateSpecialPoints(C1info(i).funcs,C1branches(i));
    Plot2dBranch(C1branches(i),'funcs',C1info(i).funcs,'ax',ax4);
end
%% The first torus bifurcation crosses a fold
% We detect it by changing stability and continue the fold of periodic
% orbits.
ind_pfold=find(diff(C1_nunst{1})==1,1,'first');
[pffuncs,pfbr,suc]=SetupPOfold(C1info(1).funcs,C1branches(1),ind_pfold,...
    'contpar',[ip.kappa1,ip.kappa2],'dir',ip.kappa1,'step',0.1,brbounds{:})
figure(3);
pfbr=br_contn(pffuncs,pfbr,100,'plotaxis',ax3);
pfbr=br_rvers(pfbr);
pfbr=br_contn(pffuncs,pfbr,400,'plotaxis',ax3);
[pf_nunst,pf_dom,pf_def,pfbr.point]=GetStability(pfbr,'funcs',pffuncs,...
    'exclude_trivial',true);
%%
Plot2dBranch(pfbr,'funcs',pffuncs,'ax',ax4);