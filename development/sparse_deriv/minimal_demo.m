%% DDE-BIFTOOL Minimal demo for continuation of local bifurcations of periodic orbits
%
% <html>
% $Id: minimal_demo.m 375 2019-09-14 14:26:50Z jansieber $
% </html>
%
% This demo illustrates how to track local bifurcations of periodic orbits
% in DDEs with contant delay, using the extension |ddebiftool_extra_psol|
% and the auxiliary functions in |ddebiftool_utilities|.
% The example is the Duffing oscillator with delayed feedback discussed in
% the large-delay limit by Yanchuk & Perlikowski in (PRE79,0462211,2009):
%
% $$x''(t)+d*x'(t)+a*x(t)+x^3+b*[x(t)-x(t-\tau)]=0$$
%
% The parameters are $(\tau,a,b,d)$ (used in this order in the |parameter|
% vector).

%% Define path and system
% First we load the folder of DDE-Biftool and its extensions into the
% Matlab path, and define the right-hand side and the delays. We create the
% structure containing the user-defined functions using |set_funcs|. We
% define the right-hand side such that it ca nbe called in vectorized form.
% 
clear
addpath('../../ddebiftool',...
    '../../ddebiftool_extra_psol',...
    '../../ddebiftool_extra_nmfm/',...
    '../../ddebiftool_utilities');
format compact
format short g
%#ok<*ASGLU>
%% Set number of delays and parameter names
% This could also be loaded from mat file |'minimal_demo_parnames.mat'|.
parnames={'tau','a','b','d'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
%% Definition of right-hand side and delay parameters
% The simplest way to define the right-hand side is by creating a
% right-hand-side function (see manual) of the format |y=rhs(xx,p)| with
% |[n,ntaup1,nvec]=size(xx)|, and a 1 x np parameter array, where |size(y)|
% is |n x nvec|. The dimension |n| is the system size, |xx(:,k,:)| is the
% value of |xx| at time |t| minus |(k-1)|th delay (|xx(:,1,:)| is the value
% of |xx| at the current time. Enabling the function to be called with a 3d
% array argument for |xx| and returning a 2d array y, permits simultaneous
% (vectorized) calls of the function in periodic-orbit continuations for
% speed-up. The function can be created in a separate m file or as an
% anonymous function handle.
%
rhs=@(xx,p)[...
     xx(2,1,:);...
    -p(ind.d)*xx(2,1,:)-p(ind.a)*xx(1,1,:)-xx(1,1,:).^3-p(ind.b)*(xx(1,1,:)-xx(1,2,:))];
%% Call to set_funcs(...)
% For DDEs with parametric delays, the r.h.s. is passed on as named input
% |'sys_rhs'|, and the function indicating which parameters are delays is
% passed on as |'sys_tau'|. Set the optional input |'x_vectorized'| in
% |set_funcs(...)| to |true| to indicate vectorization of the r.h.s.
% function.
%funcs=set_funcs('sys_rhs',rhs,'sys_tau',@()ind.tau,'x_vectorized',true); %#ok<NASGU>
%% Definition of r.h.s. and derivatives using symbolic toolbox
% The function sym_minimal_demo has been created with the symbolc toolbox
% in <gen_sym_minimal_demo.html>. The wrapper |set_symfuncs| converts this
% function into a structure as needed for further DDE-Biftool computations.
% (|x_vectorized| and |p_vectorized| are on by default in this case.)
% Uncomment the statement below to use the simple definition of |rhs|
% above.
funcs=set_symfuncs(@sym_minimal_demo,'sys_tau',@()ind.tau);
%% Set bifurcation parameter range and step size bounds
% These name-value pairs get passed on to SetupStst and will replace fields
% in the method fields of the branch controlling the computations. Note
% that the "parameter" index 0 in |'max_step'| refers to the maximum norm
% of the secant predictor.
br_args={'max_bound',[ind.tau,20;ind.b,0.6],...
    'min_bound',[ind.tau 0; ind.b,0],...
    'max_step',[0,2;ind.tau,0.2;ind.b,0.01]};
%% Definition of initial (trivial) equilibrium
% The system has only the trivial equilibrium for positive values of the
% parameter $a$. Still, we initialize a branch using the convenience
% functions 
% 
% function [br,suc]=SetupStst(funcs,varargin)
%
% Apart from the first argument (|funcs|) all arguments
% of |setupStst| are name-value pairs. Important arguments:
%
% * |'x'| (mandatory): initial guess for equilibrium state of system
% * |'parameter'| (mandatory): parameter values of initial equilibrium
% * |'contpar'| (mandatory): index of the continuation parameter
% * |'step'|: size of initial step along continuation parameter (can be
%    negative)
% * |'dir'|: index of parameter, which will be varied in initial step
% (default is |contpar(1)|) Output |br| is a branch of |stst| (equilibria)
% containing the first two points, |suc| is a flag indicating success.
% Other name-value inputs pairs get passed on to the substructures of the
% output. For example, |'minimal_real_part'| is a field controlling the
% computation of linear stability.
% 
%#ok<*NOPTS>
par(ind.tau)=0;
par(ind.a)=0.5;
par(ind.b)=0.6;
par(ind.d)=0.05;
[triv_eqs,suc]=SetupStst(funcs,'x',[0;0],'parameter',par,...
    'contpar',ind.tau,br_args{:},'minimal_real_part',-1)
%% Location and stability of trivial equilibria
% We continue the trivial equilibrium in |tau| and compute its stability
% (which changes).  The convenience function |GetStability| recomputes the
% eigenvalues if not yet present and returns as its first output
% |nunst_eqs| the number of unstable eigenvalues for bifurcation detection.
% Its first argument is the |branch| structure for which stability
% information is required. The optional argument |'funcs'| is required if
% the stability had not been computed beforehand.
disp('Trivial equilibria');
figure(1);clf;ax1=gca;
triv_eqs=br_contn(funcs,triv_eqs,60,'plotaxis',ax1);
[nunst_eqs,~,~,triv_eqs.point]=GetStability(triv_eqs,'funcs',funcs);
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
[per_orb,suc]=SetupPsol(funcs,triv_eqs,ind_hopf(2),'intervals',30,'degree',4,...
    'max_bound',[ind.tau,20],'max_step',[0,4;ind.tau,0.5]);
hold on

per_orb=br_contn(funcs,per_orb,300,'plotaxis',ax1);
amplitudes=arrayfun(@(p)norm(diff(p.profile,[],2)),per_orb.point);
per_orb.point=per_orb.point(amplitudes>1e-10); % remove singular "periodic orbits"=equilibria
[nunst_per,~,err_per,per_orb.point]=...
    GetStability(per_orb,'funcs',funcs,'exclude_trivial',true);
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
nbr2_steps=200;
hopf=br_contn(funcs,hopf,nbr2_steps,'plotaxis',ax1);
%% Criticality along Hopf bifurcation, codimension 2 points and their normal forms
% The function LocateSpecialPoints also works along Hopf branches. Its task
% can be controlled and performed manually (see
% <minimal_demo_extra_nmfm.html>). The field |hopftests.genh| contains the
% criticality of the Hopf bifurcation.
[hopf,hopftests,hc2_indices,hc2_types]=LocateSpecialPoints(funcs,hopf)
%% Plot Hopf bifurcation curve in two parameters
% This uses again the convenience function |Plot2dBranch|. Note that, if
% stability of points is not yet computed |Plot2dBranch| also needs the
% input (name-value pair) |'funcs'| with the right-hand side definition.
figure(2);clf;ax2=gca;hold(ax2,'on');
lg2=Plot2dBranch(hopf,'ax',ax2,'parameter',[ind.tau,ind.b]);
set(ax2,'box','on');
set(findobj(figure(2), 'Type', 'Legend'),'location','eastoutside');
xlabel('tau');
ylabel('b');
title('Two-parameter bifurcation diagram in tau and b');
%% Find branches of Codimension 1 bifurcations emerging from special points
% This helper function performs initialization only. By default, it finds
% the first 3 points along each branch.  The see
% <../../../ddebiftool_utilities/html/BranchFromCodimension2.html> for
% description of optional inputs and fields of output.
C1info=BranchFromCodimension2(funcs,hopf,br_args{:},'print_residual_info',1);%,'print_residual_info',1);
%% Continue all newly found branches
% in both directions, if C1info(...).reverse is true.
clear C1branches;
for i=length(C1info):-1:1
    fprintf('Continuing %s\n',[C1info(i).C1type]);
    C1branches(i)=br_contn(C1info(i).funcs,C1info(i).branch,nbr2_steps,'plotaxis',ax1);
    if C1info(i).reverse
        fprintf('Reversing %s\n',[C1info(i).C1type]);
        C1branches(i)=br_rvers(C1branches(i));
        C1branches(i)=br_contn(C1info(i).funcs,C1branches(i),nbr2_steps,'plotaxis',ax1);
    end
end
%% Compute stability and (where possible) detect special points along curve
for i=length(C1branches):-1:1
    [nunst_br{i},dom_ev{i},triv_err{i},C1branches(i).point]=...
        GetStability(C1branches(i),'funcs',C1info(i).funcs,'exclude_trivial',true);
    C1branches(i)=LocateSpecialPoints(C1info(i).funcs,C1branches(i));
end
%% plot all branches
for i=1:length(C1branches)
    lg2=Plot2dBranch(C1branches(i),'parameter',[ind.tau,ind.b],'ax',ax2,...
        'oldlegend',lg2,'funcs',C1info(i).funcs);
end
