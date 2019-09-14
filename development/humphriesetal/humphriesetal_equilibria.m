%% Equilibrium bifurcations and normal forms for example from Humpries etal (DCDS-A 2012)
%
% <html>
% $Id$
% </html>
%% Load definition of right-hand side
clear
load('humphriesetal_demo_funcs_results.mat');
%% Set continuation bounds during two-parameter continuation
eqbounds={'max_bound',[indkappa1,14; indkappa2,5],...
        'min_bound',[indkappa1,0; indkappa2,0],...
        'max_step',[indkappa1,0.1; indkappa2,0.1]};
%% Trivial equilibrium branch
% The equilibrium $x=0$ changes its stability in Hopf bifurcations
[eqbr,suc]=SetupStst(funcs,'contpar',indkappa1,'x',0,'parameter',par_ini,...
    'max_bound',[indkappa1,14],'max_step',[indkappa1,0.1]);
if ~suc
    error('equilibrium not found');
end
figure(1);clf
eqbr=br_contn(funcs,eqbr,200);
%% Stability, bifurcations and their normal forms
[eqbr,~,ind_hopf,eqc1biftype]=LocateSpecialPoints(funcs,eqbr);
figure(3);clf;ax3=gca;
lg1=Plot2dBranch(eqbr,'y',@(p)p.x(1),'ax',ax3);
set(ax3,'ylim',[-0.1,1]);
xlabel(ax3,'kappa1');
ylabel('max x');
%% Hopf bifurcation in two parameters
% Continuation parameters are $\kappa_1$ and $\kappa_2$.
nhopf=length(ind_hopf);
figure(2);clf;ax2=gca;
xlabel('kappa1');
ylabel('kappa2');
nhopfpoints=100;
fprintf('Computation of Hopf bifurcations found\n');
%% Loop over all found Hopf branches
for i=nhopf:-1:1
    %% continue branch
    [hbranch{i},suc]=SetupHopf(funcs,eqbr,ind_hopf(i),...
        'contpar',[indkappa1,indkappa2],'dir',indkappa1,'step',0.1,eqbounds{:});
    if ~suc
        error('Hopf initialization failed');
    end
    hbranch{i}=br_contn(funcs,hbranch{i},nhopfpoints,'plotaxis',ax2);
    hbranch{i}=br_rvers(hbranch{i});
    hbranch{i}=br_contn(funcs,hbranch{i},nhopfpoints,'plotaxis',ax2);
    %% compute stability, criticality and codim2 bifurcations along Hopf curves
    tic;
    [hbranch{i},~,indhoho{i},eqc2biftype{i}]=LocateSpecialPoints(funcs,hbranch{i},'stabilityfield','l1');
    ti=toc;
    %% extract criticality along Hopf curve in nf
    hpars{i}=cell2mat(arrayfun(@(x)x.parameter([indkappa1,indkappa2])',...
        hbranch{i}.point,'uniformoutput',false));
    nf{i}=arrayfun(@(x)x.nmfm.L1,hbranch{i}.point);
    fprintf('Stability changes: %d, %d subcritical, %d supercritical\ntiming=%g\n\n',...
        length(indhoho{i}),sum(nf{i}>0),sum(nf{i}<0),ti);
end
%% One additional Hopf bifurcation crosses the above Hopf curves
[hbranch{nhopf+1},suc]=SetupHopf(funcs,hbranch{1},indhoho{1}(1),...
    'contpar',[indkappa1,indkappa2],'dir',indkappa2,'step',0.1,eqbounds{:});
hbranch{nhopf+1}=br_contn(funcs,hbranch{nhopf+1},nhopfpoints,'plotaxis',ax2);
hbranch{nhopf+1}=br_rvers(hbranch{nhopf+1});
hbranch{nhopf+1}=br_contn(funcs,hbranch{nhopf+1},nhopfpoints,'plotaxis',ax2);
%% Compute stability, criticality and codim2 bifurcations along last Hopf curve
% several of Hopf-Hopf interactions should be identical to those found on
% the other Hopf branches, since these are the intersections. They may
% possibly have permuted coefficients.
[hbranch{nhopf+1},~,indhoho{nhopf+1},eqc2biftype{nhopf+1}]=...
    LocateSpecialPoints(funcs,hbranch{nhopf+1},'stabilityfield','l1');
hpars{nhopf+1}=cell2mat(arrayfun(@(x)x.parameter([indkappa1,indkappa2])',...
    hbranch{nhopf+1}.point,'uniformoutput',false));
nf{nhopf+1}=arrayfun(@(x)x.nmfm.L1,hbranch{nhopf+1}.point);
fprintf('Stability changes: %d, %d subcritical, %d supercritical\ntiming=%g\n\n',...
    length(indhoho{nhopf+1}),sum(nf{nhopf+1}>0),sum(nf{nhopf+1}<0),ti);
%% Bifucation diagram for equilibria
figure(4);clf;ax4=gca;hold(ax4,'on');
args={'ax',ax4,'stability',0.75};
lg2=Plot2dBranch(hbranch{nhopf+1},args{:});
for i=1:nhopf
    lg2=Plot2dBranch(hbranch{i},args{:},'oldlegend',lg2);
end
grid on
legend(lg2{1},lg2{2},'location','southeast');
xlabel('kappa1');
ylabel('kappa2');
%%
save(sprintf('humphries_equilibria_symbolic=%d.mat',symbolic_derivatives));
