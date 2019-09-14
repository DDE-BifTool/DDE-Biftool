%% Minimal demo - Equilibria, Hopf bifurcations, periodic orbits
% This part creates the computations that are possible with standard
% DDE-Biftool. This demo requires <minimal_demo.html> to
% have run beforehand.
%
% <html>
% $Id: minimal_demo_stst_psol.m 178 2017-03-14 00:00:01Z jansieber $
% </html>
%
%% Definition of initial (trivial) equilibrium
% The system has only the trivial equilibrium for positive values of the
% parameter $a$. Still, we initialize a branch using the convenience
% functions 
% 
% function [br,suc]=SetupStst(funcs,varargin)
%
% Apart from the first argument (|funcs|) all arguments
% of |gen_stst| are name-value pairs. Important arguments:
%
% * |'x'|: initial guess for equilibrium state of system
% * |'parameter'|: parameter values of initial equilibrium
% * |'contpar'|: index of the continuation parameter
% * |'step'|: size of initial step along continuation parameter (can be
%    negative)
% * |'dir'|: index of parameter, which will be varied in initial step
% (default is |contpar(1)|)
% Output |br| is a branch of |stst| (equilibria) containing the first two
% points, |suc| is a flag indicating success. Other name-value inputs pairs
% get passed on to the substructures of the output.
% 
%#ok<*NOPTS>
triv_eqs=SetupStst(funcs,'x',[0;0],'parameter',[0.2,0.5,0.6,0.05],...
    'contpar',ind.tau,'max_step',[ind.tau,0.3],'max_bound',[ind.tau,20],...
    'newheuristics_tests',0,'minimal_real_part',-1);
%% Stability of trivial equilibria
% We continue the trivial equilibrium in |tau| and compute its stability
% (which changes).  The convenience function |GetStability| recomputes the
% eigenvalues if not yet present and returns as its first output
% |nunst_eqs| the number of unstable eigenvalues for bifurcation detection.
% Its first argument is the |branch| structure for which stability
% information is required. The optional argument |'funcs'| is required if
% the stability had not been computed beforehand.
clear
close all
addpath('../../ddebiftool',...
    '../../ddebiftool_extra_psol',...
    '../../ddebiftool_extra_nmfm/',...
    '../../ddebiftool_utilities');
format compact
format short g
%% Set number of delays and parameter names
ntau=1;
parnames={'tau','a','b','d'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
funcs=set_symfuncs(@sym_minimal_demo,'sys_tau',@()ind.tau);
%% Choose one of |fnum|, |fmaple| or |fsymbolic|
disp('Trivial equilibria');
figure(1);clf;ax1=gca;
triv_eqs=SetupStst(funcs,'x',[0;0],'parameter',[0.2,0.5,0.6,0.05],...
    'contpar',ind.tau,'max_step',[ind.tau,0.3],'max_bound',[ind.tau,20],...
    'newheuristics_tests',0,'minimal_real_part',-1);
triv_eqs=br_contn(funcs,triv_eqs,60,'plotaxis',ax1);
nunst_eqs=GetStability(triv_eqs,'funcs',funcs);
%% Find Hopf bifurcations and their criticality
% The convenience function LocateSpecialPoints detects bifurcations along
% the branch and comptues their criticality. This function is only
% implemented for equilibrium type branches (including Hopf or fold
% branches).
[triv_eqs,~,ind_hopf,stst_bifs]=LocateSpecialPoints(funcs,triv_eqs) 
%% Branch off at 2nd Hopf bifurcation to find periodic orbits
% We find the first point at which the number of unstable eigenvalues is
% greater than 2. The convenience function |SetupPsol| performs the
% initial corrections for the first two points along the branch of periodic
% orbits.
%
% function [per,suc]=SetupPsol(funcs,branch,ind,varargin)
%
% The first input |funcs| is the system definition structure, the second
% input |branch| is the branch from which one wants to branch off, and
% |ind| is the index of the point along this branch. The point can be either
% of type |stst| with stability information or aof type |hopf|. The other
% arguments are name-value pairs. Important inputs:
%
% * |'contpar'| (default |[]|): index of continuation parameter
% * |'radius'| (default |0.01|): amplitude of first periodic orbit on branch
% * |'degree'| (default |3|): degree of collocation polynomials used for
%     periodic orbits
% * |'intervals'| (default |20|): number of collocation intervals
% * |'hopfcorrection'| (default |true|): call p_correc to find the Hopf
%     point accurately
%
% Output |per| is a branch of periodic orbits with the first two points,
% |suc| is a flag indicating success. Other name-value input pairs get
% passed on to the substructures of |per|.
%
% Note how the |'max_step'| field of the |per.parameter| permits now to use
% index 0 to indicate a restriction on the overall secant length (measured
% with p_norm).
%
% Subsequently the newly created branch of periodic orbits is continued and
% its stability is determined. The 3rd output of GetStability shows the
% distance of the Floquet multiplier closest to 1 to its theoretical value
% of 1. The first ouyput gives the number of unstabble Flouqet multipliers
% ignoreing the trivial one (enforced by the optional input
% |'exclude_trivial'|).

disp('Branch off at 2nd Hopf bifurcation');
fprintf('Initial correction of periodic orbits at Hopf:\n');
[per_orb,suc]=SetupPsol(funcs,triv_eqs,ind_hopf(2),'intervals',20,'degree',4,...
    'max_bound',[ind.tau,20],'max_step',[0,2;ind.tau,0.5]);
hold on
per_orb=br_contn(funcs,per_orb,300,'plotaxis',ax1);
per_orb.point=per_orb.point(arrayfun(@(x)norm(x.profile)>1e-10,per_orb.point));
[nunst_per,~,err_per,per_orb.point]=...
    GetStability(per_orb,'funcs',funcs,'exclude_trivial',true);
%% 1d bifurcation diagram
% We can manually extract all information from the branches and plot and
% color the branches according to their stability. Colors indicate number of
% unstable Floquet multipliers.
colors='kbcrm';
xp1=arrayfun(@(x)x.parameter(ind.tau),per_orb.point);
yp1=arrayfun(@(x)max(x.profile(1,:)),per_orb.point);
pl1={};
for i=0:max(nunst_per)
    pl1=[pl1,{xp1(nunst_per==i),yp1(nunst_per==i),[colors(i+1),'o']}]; %#ok<AGROW>
end
figure(2);
plot(pl1{:},'linewidth',2);
grid on
xlabel('tau');
ylabel('max(x)');
legend({'0','1','2','3','4'})
title('1d bif diagram of p.o''s in tau, color=stability');
%% Automatic plot of single-parameter bifurcation diagram
% The convenience function |Plot2dBranch| makes a few sensible default
% guesses to produce a reasonable plot for each branch. Its output |lg| is
% legend text that can be fed into subsequent calls to produce a joint
% legend.
figure(4);clf;ax4=gca;hold(ax4,'on');
lg1=Plot2dBranch(triv_eqs,'ax',ax4);
lg1=Plot2dBranch(per_orb,'ax',ax4,'oldlegend',lg1);
grid(ax4,'on');
xlabel(ax4,'tau');
ylabel(ax4,'max(x)');
%% Continue 2nd Hopf bifurcation
% Similar to |SetupPsol| the convenience function |SetupHopf| creates
% the initial Hopf branch. Its first arguments are |funcs|, the branch
% along which the Hopf bifurcation was detected (here |triv_eqs|), and the
% index of the point near which the Hopf bifurcation was detected.
%
% function [hbranch,suc]=SetupHopf(funcs,branch,ind,varargin)
%
% Important parameters:
%
% * |'contpar'|: continuation parameters (vector of length >=2), 
% * |'dir'|: index of parameter, which is varied at initial step. The
% default is [], which means that only one point on the branch is computed.
% This is useful if one wants to correct only a single Hopf point.
% * |'step'|: initial step along branch (default |1e-3|)
% * |'excudefreqs'|: list of frequencies that should be excluded (default
% []). The initial guess for the Hopf frequency is the complex conjugate
% pair closest to the imaginary axis, after one takes away a pair of
% eigenvalues for each frequency listed in |excludefreqs|.
%
% All other name-value pairs can be used to replace fields in the
% structures of the Hopf branch. Otherwise, the output |branch2| inherits
% all values from the input |branch|. The subsequent continuation computes
% toward smaller |b| for 80 steps.
hopfopts={'step',-1e-3,'max_step',[ind.tau,0.2;ind.b,0.01],'max_bound',[ind.b,0.6]};
hopf=SetupHopf(funcs,triv_eqs,ind_hopf(2),'contpar',[ind.b,ind.tau],'dir',ind.b,...
    hopfopts{:});
figure(1);clf;ax1=gca;
hopf=br_contn(funcs,hopf,200,'plotaxis',ax1);
%% Continue 1st Hopf bifurcation
% We also find the first hopf bifurction of the branch |triv_eqs| and
% continue that, too.
hopf1=SetupHopf(funcs,triv_eqs,ind_hopf(1),'contpar',[ind.b,ind.tau],'dir',ind.b,...
    hopfopts{:});
hopf1=br_contn(funcs,hopf1,200,'plotaxis',ax1);
%% Criticality along Hopf bifurcation, codimension 2 points and their normal forms
% The function LocateSpecialPoints also works along Hopf branches. Its task
% can be controlled and performed manually (see
% <minimal_demo_extra_nmfm.html>). The field |hopftests.genh| contains the
% criticality of the Hopf bifurcation
[hopf,hopftests,hc2_indices,hc2_types]=LocateSpecialPoints(funcs,hopf)
[hopf1,hopf1tests,h1c2_indices,h1c2_types]=LocateSpecialPoints(funcs,hopf1)
%% Automated plot
% We may use |Plot2dBranch| to get a rough bifurcation diagram.
figure(3);clf;ax3=gca;hold(ax3,'on');
lg2=Plot2dBranch(hopf,'ax',ax3);
lg2=Plot2dBranch(hopf1,'ax',ax3,'oldlegend',lg2);
grid(ax3,'on');
xlabel(ax3,'b');
ylabel(ax3,'tau');

%% Save and continue
% For continuation of folds and torus bifurcations of periodic orbits, see
% <minimal_demo_extra_psol.html>. Final results in
% <minimal_demo_plot_2dbif.html>.
save('minimal_demo_stst_psol_results.mat')