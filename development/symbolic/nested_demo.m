%% Test state-dependent delay equations with many levels of nesting
%
% <html>
% $Id: nested_demo.m 174 2017-03-10 21:59:17Z jansieber $
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
close all
addpath('../../ddebiftool/',...
    '../../ddebiftool_extra_psol/',...
    '../../ddebiftool_extra_nmfm/',...
    '../../ddebiftool_utilities/');
format compact
%#ok<*NOPTS>
%% Equation definition
% Change |ntau| to change the level of nesting.
dim=1;
ntau=3;
rhs=@(x,p)-x(1,ntau+1,:);
sys_ntau=@()ntau;
tau=@(nr,x,p)p(1)+x(1,nr,:)+p(2)*x(1,nr,:).^2;
drhsdir{1}=@(x,p,dx,dp)-dx(1,ntau+1,:);
drhsdir(2:5)={@(x,p,dx,dp)0*x(:,1,:)};
dtaudir{1}=@(nr,x,p,dx,dp)dp(1)+dx(1,nr,:)+dp(2)*x(1,nr,:).^2+p(2)*2*x(1,nr,:).*dx(1,nr,:);

%% Different levels of problem definition
%% DefinitinoSymbolic 
% One may use the symbolic toolbox to define the right-hand side and delays
% (see <generate_sym_nested.html>). For this case the wrapper
% |set_symfuncs| convertes the symbolic toolbox output into functions
% usable by DDE-Biftool.
fsym=set_symfuncs(sprintf('sym_nested_ntau_%d',ntau));
%% I maple is available
args={'sys_dirderi',str2func(sprintf('nested_dirderi%d',ntau))};
fmaple=set_funcs('sys_rhs',rhs,'sys_ntau',sys_ntau,'sys_tau',tau,args{:},'x_vectorized',true);
fnum=set_funcs('sys_rhs',rhs,'sys_ntau',sys_ntau,'sys_tau',tau,'x_vectorized',true);
funcs=fsym;
fdir=set_funcs('sys_rhs',rhs,'sys_ntau',sys_ntau,'sys_tau',tau,...
    'sys_dirderi',drhsdir,'sys_dirdtau',dtaudir,'x_vectorized',true);
pbounds={'max_bound',[1,2; 2,2],'min_bound',[1,-1e-3; 2,-2],'max_step',[1,0.1; 2,0.1]};
%% Branch of trivial Equilibria
% x=0 is always the only equlibrium.
[eqbr,suc]=SetupStst(funcs,'contpar',1,'x',zeros(dim,1),'parameter',[0,0],...
    pbounds{:});
if ~suc
    error('equilibrium not found');
end
figure(1);clf;ax1=gca;
eqbr=br_contn(funcs,eqbr,100);
%% Stability of equilibria
% The family of trivial equlibria loses stability in a Hopf bifurcation at
% |p(1)=pi/2|.
[eqbr,eqtests,indhopf,eqbiftype]=LocateSpecialPoints(funcs,eqbr)
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
per=br_contn(funcs,per,20,'plotaxis',ax1);
%% Stability of periodic orbits
% The only source of instability is the fold such that periodic orbits are
% either stable or order-1 unstable.
[pernunst,dom,triv_defect,per.point]=...
    GetStability(per,'exclude_trivial',true,'funcs',funcs); %#ok<ASGLU>
fprintf('maximum error of trivial Floquet multiplier: %g\n',max(abs(triv_defect)));
%% Profiles of periodic orbits
ppars=arrayfun(@(x)x.parameter(1),per.point);
pmeshes=cell2mat(arrayfun(@(x)x.mesh(:),per.point,'uniformoutput',false));
pprofs=cell2mat(arrayfun(@(x)x.profile(1,:)',per.point,'uniformoutput',false));
figure(2);clf
plot(pmeshes,pprofs,'.-');
xlabel('t/T');
ylabel('x');
title('profiles');
%% Plot 1d bifurcation diagram
cla(ax1);hold(ax1,'on');
lg1=Plot2dBranch(eqbr,'ax',ax1);
lg1=Plot2dBranch(per,'ax',ax1,'oldlegend',lg1);
lg1=Plot2dBranch(per,'ax',ax1,'oldlegend',lg1,'y',@(p)min(p.profile));
xlabel('p(1)');
ylabel('min x, max x');
%% Save data
save(sprintf('sd_basic_per%d.mat',ntau));
%% Continue Hopf bifurcation and check criticality along the branch
% Note that one should specify |'correc'| false, since otherwise, the
% initial correction will attempt to correct using contpar(2) and parameter
% 2 has no influence on the Hopf bifurcation.
hbranch=SetupHopf(funcs,eqbr,indhopf,'contpar',[1,2],'correc',false,'dir',2,...
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
% the fold is not present for large |p(2)|. We start from the fold
% (last point where the number of Floquet multipliers changes by one).
disp('Fold initialization')
pf_ind0=find(abs(diff(pernunst))==1,1,'last');
[pfuncs,pbr,suc]=SetupPOfold(funcs,per,pf_ind0,'contpar',[1,2],'dir',1,...
    'step',1e-3,pbounds{:},'print_residual_info',1);
if ~suc
    error('initialization of folds of periodic orbits failed');
else
    disp('PO folds initialized');
end
%% Fold of periodic orbits
% continued in two parameters |p(1:2)|.
figure(3);
pbr=br_contn(pfuncs,pbr,20,'plotaxis',ax3);
pbr=br_rvers(pbr);
pbr=br_contn(pfuncs,pbr,20,'plotaxis',ax3);
%% Stability of periodic orbits at fold
% We compute the stability of the periodic orbits at the fold excluding the
% two Floquet multipliers closest to unity in |pfnunst|. All orbits are
% stable transversal to the fold direction.
[pfnunst,dom,triv_defect,pbr.point]=GetStability(pbr,'exclude_trivial',true,'funcs',pfuncs);
display(pfnunst.','transversally unstable Floquet multipliers at fold=');
%% Plot two-parameter bifurcation diagram in parameter plane
cla(ax3);hold(ax3,'on');
lg2=Plot2dBranch(hbranch,'ax',ax3);
lg2=Plot2dBranch(pbr,'ax',ax3,'oldlegend',lg2,'funcs',pfuncs);
xlabel(ax3,'p(1)');
ylabel(ax3,'p(2)');
%% Plot two-parameter bifurcation diagram in 3d
% axes: |p(1)|, |p(2)| and |max(x)| and |min(x)|. The Hopf bifurcation is
% independent of |p(2)|
pfs=pfuncs.get_comp(pbr,'solution');
pfpars=cell2mat(arrayfun(@(x)x.parameter(1:2)',pfs,'uniformoutput',false));
pfmeshes=cell2mat(arrayfun(@(x)x.mesh(:),pfs,'uniformoutput',false));
pfprofs=cell2mat(arrayfun(@(x)x.profile(1,:)',pfs,'uniformoutput',false));
hpars=per.point(1).parameter;
figure(5);clf
pl3=plot3(pfpars(1,:),pfpars(2,:),max(pfprofs),'k.-',...
    pfpars(1,:),pfpars(2,:),min(pfprofs),'k.-',...
    hpars(1)+0*pfpars(1,:),pfpars(2,:),0*pfpars(1,:),'.-');
grid on
legend(pl3([1,3]),{'Fold of periodic orbits','Hopf bifurcation'},'location','south')
xlabel('p(1)');
ylabel('p(2)');
zlabel('max(x)-min(x)');
%% Save all data
save('NestedPOfold.mat');
