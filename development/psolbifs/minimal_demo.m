%% DDE-BIFTOOL Minimal demo for continuation of local bifurcations of periodic orbits
%
% <html>
% $Id: minimal_demo.m 362 2019-07-14 15:49:40Z jansieber $
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
base=[pwd(),'/../../'];
addpath([base,'ddebiftool/'],...
    [base,'ddebiftool_extra_psol/'],...
    [base,'ddebiftool_utilities/'],...
    [base,'ddebiftool_extra_nmfm']);
format compact
format short g
load('minimal_demo_parnames.mat');
load('minimal_demo_stst_psol_results.mat');
%#ok<*ASGLU>
%%
funcs=set_symfuncs(@sym_minimal_demo,'sys_tau',@()ind.tau);
%% Set bifurcation parameter range and step size bounds
% These name-value pairs get passed on to SetupStst and will replace fields
% in the method fields of the branch controlling the computations. Note
% that the "parameter" index 0 in |'max_step'| refers to the maximum norm
% of the secant predictor.
br_args={'max_bound',[ind.tau,20;ind.b,0.6],...
    'min_bound',[ind.tau 0; ind.b,0],...
    'max_step',[0,2;ind.tau,0.2;ind.b,0.01]};
%% Try to use psolfold for continuation
ipf=find(nunst_per==0,1,'first');
pf0=per_orb.point(ipf);
pfoldini=dde_psolfold_init(funcs,pf0,[],per_orb.method);
method=dde_psolfold_method(per_orb.method);
method.point.print_residual_info=1;

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
