%% Test state-dependent delay equations with many levels of nesting
%
% <html>
% $Id: nested_demo.m 161 2017-02-23 21:58:11Z jansieber $
% </html>
%
% The equation is 
% 
% $$x'(t)=-x(t+\theta_k)$$
%
% where $\theta_0(t)=0$ and
% 
% $$ \theta_{j+1}(t)=-p_1-x(t+\theta_j(t))-p_2 x(t+\theta_j(t))^2 $$
% 
% resulting in a nesting level of delays of order |k|. The parameter |p(1)|
% controls the delay at the Hopf bifurcation, |p(2)| controls the
% criticiality of the Hopf bifurcation (at least for |k=3| this is
% demonstrated below).
%
%%
clear
addpath('../../ddebiftool/',...
    '../../ddebiftool_extra_psol/',...
    '../../ddebiftool_extra_nmfm/',...
    '../../ddebiftool_utilities/');
format compact
%#ok<*NOPTS>
%% Equation definition
% Change |ntau| to change the level of nesting.
parnames={'tau1','p'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
dim=1;
ntau=3;
par0([ind.tau1,ind.p])=[0,0];
%rhs=@(x,p)-x(1,ntau+1,:);
%sys_ntau=@()ntau;
%tau=@(nr,x,p)p(1,ind.tau1,:)+x(1,nr,:)+p(1,ind.p,:)*x(1,nr,:).^2;
%funcs=set_funcs('sys_rhs',rhs,'sys_ntau',sys_ntau,'sys_tau',tau,'x_vectorized',true);
funcs=set_symfuncs(@sym_nested);
pbounds={...
    'max_bound',[ind.tau1,2;     ind.p,2],...
    'min_bound',[ind.tau1,-1e-3; ind.p,-2],...
    'max_step', [ind.tau1,0.1;   ind.p,0.1]};
%% Branch of trivial Equilibria
% x=0 is always the only equlibrium.
[eqbr,suc]=SetupStst(funcs,'contpar',ind.tau1,'x',zeros(dim,1),'parameter',par0,...
    pbounds{:});
if ~suc
    error('equilibrium not found');
end
figure(1);clf;ax1=gca;
eqbr=br_contn(funcs,eqbr,100);
%% Stability of equilibria
% The family of trivial equlibria loses stability in a Hopf bifurcation at
% |p(1)=pi/2|.
[eqbr,eqtests,indhopf,eqbiftype]=LocateSpecialPoints(funcs,eqbr,'debug',true)
%% Branch off at Hopf bifurcation
% At the Hopf bifurcation a family of periodic orbits branches off. Its
% stability close to the Hopf bifurcation depends on the level of nesting |ntau|.
[per,suc]=SetupPsol(funcs,eqbr,indhopf,'degree',3,'intervals',30,...
    'print_residual_info',1,pbounds{:});
if ~suc
    error('initialization of periodic orbits failed');
end
%% Periodic orbits continued in |p(1)|
% The family of periodic orbits folds back and forth for |ntau=3|.
per.parameter.max_step=[0,0.01];
per=br_contn(funcs,per,30,'plotaxis',ax1);
per.parameter.max_step=[];
per=br_contn(funcs,per,40,'plotaxis',ax1);
%% Stability of periodic orbits
% The only source of instability is the fold such that periodic orbits are
% either stable or order-1 unstable.
[pernunst,dom,triv_defect,per.point]=...
    GetStability(per,'exclude_trivial',true,'funcs',funcs); %#ok<ASGLU>
fprintf('maximum error of trivial Floquet multiplier: %g\n',max(abs(triv_defect)));
%% Profiles of periodic orbits
ppars=arrayfun(@(x)x.parameter(ind.tau1),per.point);
pmeshes=cell2mat(arrayfun(@(x)x.mesh(:),per.point,'uniformoutput',false));
pprofs=cell2mat(arrayfun(@(x)x.profile(1,:)',per.point,'uniformoutput',false));
figure(2);clf
plot(pmeshes,pprofs,'.-');
xlabel('t/T');
ylabel('x');
title('profiles');
%% Plot 1d bifurcation diagram
figure(1);clf;ax1=gca;hold(ax1,'on');
Plot2dBranch(eqbr,'ax',ax1);
Plot2dBranch(per,'ax',ax1);
Plot2dBranch(per,'ax',ax1,'y',@(p)min(p.profile));
xlabel(parnames{ind.tau1});
ylabel('min x, max x');
%% Save data
save(sprintf('sd_basic_per%d.mat',ntau));
%% Continue Hopf bifurcation and check criticality along the branch
hbranch=SetupHopf(funcs,eqbr,indhopf,'contpar',[ind.tau1,ind.p],'dir',ind.p,...
    'step',0.1,pbounds{:});
figure(3);clf;ax3=gca;
hbranch=br_contn(funcs,hbranch,30,'plotaxis',ax3);
hbranch=br_rvers(hbranch);
hbranch=br_contn(funcs,hbranch,30,'plotaxis',ax3);
[hbranch,hopftests,hc2_ind,hc2_types]=LocateSpecialPoints(funcs,hbranch) 
%% Plot criticality of Hopf bifurcation along Hopf line
figure(4);clf;
plot(arrayfun(@(x)x.parameter(2),hbranch.point),hopftests.genh,'.-');
grid on
xlabel('p2');
ylabel('L1');
%% Fold of periodic orbits
% A two-parameter continuation of the folds of periodic orbits shows that
% the fold is not present for large |p(2)|. We start from the generalized
% Hopf point detected along the Hopf line.
disp('PO Fold initialization')
C1info=BranchFromCodimension2(funcs,hbranch,pbounds{:},'print_residual_info',1);
figure(3);
pfuncs=C1info.funcs;
pbr=br_contn(C1info.funcs,C1info.branch,50,'plotaxis',ax3);
%% Stability of periodic orbits at fold
% We compute the stability of the periodic orbits at the fold excluding the
% two Floquet multipliers closest to unity in |pfnunst|. All orbits are
% stable transversal to the fold direction.
[pfnunst,dom,triv_defect,pbr.point]=GetStability(pbr,'exclude_trivial',true,'funcs',pfuncs);
celldisp({'transversally unstable Flqouet multipliers at fold=',pfnunst.'});
%% Plot two-parameter bifurcation diagram in parameter plane
cla(ax3);hold(ax3,'on');
lg2=Plot2dBranch(hbranch,'ax',ax3);
lg2=Plot2dBranch(pbr,'ax',ax3,'oldlegend',lg2,'funcs',pfuncs);
xlabel(ax3,parnames{ind.tau1});
ylabel(ax3,parnames{ind.p});
%% Plot two-parameter bifurcation diagram in 3d
% axes: |p(1)|, |p(2)| and |max(x)| and |min(x)|. The Hopf bifurcation is
% independent of |p(2)|
pfs=pfuncs.get_comp(pbr,'solution');
pfpars=cell2mat(arrayfun(@(x)x.parameter(1:2)',pfs,'uniformoutput',false));
pfmeshes=cell2mat(arrayfun(@(x)x.mesh(:),pfs,'uniformoutput',false));
pfprofs=cell2mat(arrayfun(@(x)x.profile(1,:)',pfs,'uniformoutput',false));
hpars=cell2mat(arrayfun(@(x)x.parameter(1:2)',hbranch.point,'uniformoutput',false));
figure(5);clf
pl3=plot3(pfpars(ind.tau1,:),pfpars(ind.p,:),max(pfprofs),'k.-',...
    pfpars(ind.tau1,:),pfpars(ind.p,:),min(pfprofs),'k.-',...
    hpars(ind.tau1,:),hpars(ind.p,:),0*hpars(ind.tau1,:),'.-');
grid on
legend(pl3([1,3]),{'Fold of periodic orbits','Hopf bifurcation'},'location','south')
xlabel(parnames{ind.tau1});
ylabel(parnames{ind.p});
zlabel('max(x)-min(x)');
%% Save all data
save('NestedPOfold.mat');
